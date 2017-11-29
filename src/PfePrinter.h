/*
 * ParOSol: a parallel FE solver for trabecular bone modeling
 * Copyright (C) 2011, Cyril Flaig
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef PFEPRINTER_H
#define PFEPRINTER_H

#include "Config.h"
#include "OctreeGrid.h"
#include "Postprocessing.h"
#include "GWriter.hpp"

#include <eigen3/Eigen/Core>

//! A class to write the result in the format that ParFE uses.

/*! PfePrinter print the results in the same format as ParFE. In this way,
 *  ParOSol can be used as a fast mesher for ParFE. 
 */


template <class T>
class PfePrinter {
    public:
        PfePrinter(std::string filename, OctreeGrid<T> &grid):_MyPID(grid.GetPID()), _Size(grid.GetNrCPU()),_filename(filename), _grid(grid)
    {
        Writer = new HDF5_GWriter(filename, MPI_COMM_WORLD);
    }
        int _MyPID;
        int _Size;
        ~PfePrinter()
        {
            Writer->Close();
            delete Writer;
        }

        void OctKey_to_Coord(long key, int &x, int &y, int &z)
        {
            x = y = z = 0;
            for(int i = 0; i < 16; i++) {
                x += (key & 1) << i;
                key = key >> 1;
                y += (key & 1) << i;
                key = key >> 1;
                z += (key & 1) << i;
                key = key >> 1;
            }
        }


        void PrintGrid() {
            PrintCoord("Coordinates");
            PrintElems("Elements");
            PrintMaterialIDs("Material IDs");
            PrintBC();
        }

        //! Print all Boundary Conditions in to the fixed node Data Set.
        void PrintBC() {
            Writer->Select("/Boundary conditions");
            int my_n_fixed_nodes = _grid.bc.FixedNodes_Ind.size();
            int n_fixed_nodes;
            MPI_Allreduce(&my_n_fixed_nodes, &n_fixed_nodes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            int *sense_fixed = new int[my_n_fixed_nodes];
            double *value_fixed = new double[my_n_fixed_nodes];
            int *node_numbers_fixed = new int[my_n_fixed_nodes];
            int myNFixedNodes = 0;
            double bc_val, abs_bc_val;

            t_octree_key offset = _grid.GetNodeOffset();
            for(int i =0;i<my_n_fixed_nodes;i++) {
                bc_val = _grid.bc.FixedNodes[i];
                abs_bc_val = fabs(bc_val);
                if (abs_bc_val > 1e-14)	{
                    node_numbers_fixed[myNFixedNodes] = _grid.bc.FixedNodes_Ind[i]/3 +offset+1;
                    sense_fixed[myNFixedNodes] = _grid.bc.FixedNodes_Ind[i]%3+1;
                    value_fixed[myNFixedNodes] = bc_val;
                    myNFixedNodes++;
                }
                // TODO: needed?
                /*
                   if (value[i] == 0)
                   value[i] = 1e-16;
                   */
            }

            int my_fixed_nodes_offset, N_actual_fixed_nodes;
            MPI_Scan(&myNFixedNodes, &my_fixed_nodes_offset,1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&myNFixedNodes, &N_actual_fixed_nodes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            my_fixed_nodes_offset -= myNFixedNodes; //include it selfs

            Writer->Write("Fixed nodes size", (int)N_actual_fixed_nodes);
            if(N_actual_fixed_nodes > 0)	{
                Writer->Write("Fixed nodes", "Node number", node_numbers_fixed,"Sense", sense_fixed, "Value", value_fixed, N_actual_fixed_nodes, myNFixedNodes, 1,1,1, my_fixed_nodes_offset );
            }
            Writer->Write("Loaded nodes size", (int)0);
            Writer->Write("Restrained nodes size", (int)0);

            delete[] sense_fixed;
            delete[] value_fixed;
            delete[] node_numbers_fixed;
        }

        void PrintMaterialIDs(std::string dset) {
            Writer->Select("/Mesh");
            Eigen::VectorXi MatIDs(_grid.GetNrElem());
            t_index i =0;
            for(_grid.initIterateOverElements(); _grid.TestIterateOverElements(); _grid.IncIterateOverElements())	{
                MatIDs[i] = _grid.GetElementWeight();
                i++;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            Writer->Write(dset, MatIDs.data(), _grid.GetNrElemGlobal(),_grid.GetNrElem(), 1, _grid.GetElemOffset());
        }


        void PrintCoord(std::string dset) {
            Writer->Select("/Mesh");


            t_octree_key k;
            std::vector<OctreeNode> &grid = _grid.GetOctGrid();
            std::vector<OctreeNode>::iterator iter; 
            int x =0,y=0,z=0;
            double res[3];
            _grid.GetRes(res);
            Eigen::VectorXd coord(_grid.GetNrPrivateNodes()*3);
            long i=0;
            for(iter = grid.begin(); iter != _grid._GridIteratorEnd; ++iter ) {
                k = iter->key;
                OctKey_to_Coord(k, x, y, z);
                coord[i++] = x*res[0];
                coord[i++] = y*res[1];
                coord[i++] = z*res[2];
            }
            Writer->Write(dset, coord.data(), _grid.GetNrNodesGlobal(),_grid.GetNrPrivateNodes(), 3, _grid.GetNodeOffset());
        }

        void PrintElems(std::string dset) {
            Writer->Select("/Mesh");

            //print element to node
            int *elems = new int[_grid.GetNrElem()*8];

            t_index local_nodes[8];

            //Quick an dirty hack:
            //compute the offset with double
            Eigen::VectorXd ind(_grid.GetNrDofs());
            _grid.Recv_import_Ghost(ind);
            ind.setZero(_grid.GetNrDofs());
            t_octree_key offset = _grid.GetNodeOffset();
            for (unsigned int i = 0; i < _grid.GetNrNodes(); i++) {
                ind[3*i] = i + offset+1;
            }
            _grid.Send_import_Ghost(ind);
            _grid.Wait_import_Ghost();

            t_index e =0;
            for(_grid.initIterateOverElements(); _grid.TestIterateOverElements(); _grid.IncIterateOverElements()){
                _grid.SearchIndexes(local_nodes); 
                for(int i =0; i <8; i++)
                    elems[e*8+i] = (int) ind[3*local_nodes[i]];
                e++;
            }
            Writer->Write(dset, elems, _grid.GetNrElemGlobal(),_grid.GetNrElem(), 8, _grid.GetElemOffset());
            delete[] elems;
        }

        void PrintAll(Eigen::VectorXd &x, Eigen::VectorXd &force, Eigen::VectorXd &res) {

            PrintGrid();

            PostProcess<OctreeGrid<T> > post(_grid);
            Eigen::VectorXd m, s, eff;
            Eigen::VectorXd stresses, strains;
            post.ComputeStressAndStrain(x,m,s,eff,stresses, strains);
            MPI_Barrier(MPI_COMM_WORLD);

            Writer->Select("/Solution");
            Writer->Write("Nodal displacements", x.data(), _grid.GetNrNodesGlobal(),_grid.GetNrPrivateNodes(), 3, _grid.GetNodeOffset());
            Writer->Write("VonMises", m.data(), _grid.GetNrElemGlobal(),_grid.GetNrElem(), 1, _grid.GetElemOffset());
            Writer->Write("SED", s.data(), _grid.GetNrElemGlobal(),_grid.GetNrElem(), 1, _grid.GetElemOffset());
            Writer->Write("EFF", eff.data(), _grid.GetNrElemGlobal(),_grid.GetNrElem(), 1, _grid.GetElemOffset());
            Writer->Write("Element stress", stresses.data(), _grid.GetNrElemGlobal(),_grid.GetNrElem(), 6, _grid.GetElemOffset());
            Writer->Write("Element strain", strains.data(), _grid.GetNrElemGlobal(),_grid.GetNrElem(), 6, _grid.GetElemOffset());
            Writer->Write("Nodal forces", force.data(), _grid.GetNrNodesGlobal(),_grid.GetNrPrivateNodes(), 3, _grid.GetNodeOffset());
        }

        std::string _filename;

        OctreeGrid<T> &_grid;

        HDF5_GWriter *Writer;

};
#endif /* PFEPRINTER */
