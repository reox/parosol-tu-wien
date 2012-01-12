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


#ifndef HDF5PRINTER_H
#define HDF5PRINTER_H

#include "Config.h"
#include "OctreeGrid.h"
#include "Postprocessing.h"
#include "GWriter.hpp"

#include <eigen2/Eigen/Core>
USING_PART_OF_NAMESPACE_EIGEN

//! This class prints the grid and result into a HDF5 file.

/*! The printer writes the whole into a single HDF5-file. This file can be read by ParaView.
 *  To read it with ParaView a XML-file is used. This file can be generate with a helper script.
 */

template <class T>
class HDF5Printer {
	public:
		HDF5Printer(std::string filename, OctreeGrid<T> &grid):_MyPID(grid.GetPID()), _Size(grid.GetNrCPU()),_filename(filename), _grid(grid)
		{
		  Writer = new HDF5_GWriter(filename, MPI_COMM_WORLD);
		}
		int _MyPID;
		int _Size;
		~HDF5Printer()
		{
		  Writer->Close();
		  delete Writer;
		}
		
		void OctKey_to_Coord(long key, int &x, int &y, int &z)
		{
			x = y = z = 0;
			for(int i = 0; i < 20; i++) {
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
            PrintEmoduli();
		}
		
		void PrintCoord(std::string dset) {
		  Writer->Select("/Mesh");
		  
		  
			t_octree_key k;
			std::vector<OctreeNode> &grid = _grid.GetOctGrid();
			std::vector<OctreeNode>::iterator iter; 
			int x =0,y=0,z=0;
			double res[3];
			_grid.GetRes(res);
			VectorXf coord(_grid.GetNrPrivateNodes()*3);
			long i=0;
			T keys;
			t_octree_key tmp;
			for(iter = grid.begin(); iter != _grid._GridIteratorEnd; ++iter ) {
				k = iter->key;
				OctKey_to_Coord(k, x, y, z);
				tmp = keys(x,y,z);
				if (tmp != k) {
				  PCOUT(_MyPID, "ERROR not the same key\n")
				  break;
				}
				coord[i++] = x*res[0];
				coord[i++] = y*res[1];
				coord[i++] = z*res[2];
			}
			Writer->Write(dset, coord.data(), _grid.GetNrNodesGlobal(),_grid.GetNrPrivateNodes(), 3, _grid.GetNodeOffset());
		}

		void PrintElems(std::string dset) {
			Writer->Select("/Mesh");

			//print element to node
			t_octree_key *elems = new t_octree_key[_grid.GetNrElem()*8];
			
			t_index local_nodes[8];
			
			//Quick an dirty hack:
			//compute the offset with double
			VectorXd ind(_grid.GetNrDofs());
			_grid.Recv_import_Ghost(ind);
			ind.setZero(_grid.GetNrDofs());
			t_octree_key offset = _grid.GetNodeOffset();
			MPI_Barrier(MPI_COMM_WORLD);
			for (t_index i = 0; i < _grid.GetNrNodes(); i++) {
				ind[3*i] = i + offset;
			}
			MPI_Barrier(MPI_COMM_WORLD);
			_grid.Send_import_Ghost(ind);
			_grid.Wait_import_Ghost();
			MPI_Barrier(MPI_COMM_WORLD);

			t_index e =0;
			for(_grid.initIterateOverElements(); _grid.TestIterateOverElements(); _grid.IncIterateOverElements()){
				_grid.SearchIndexes(local_nodes); 
				for(int i =0; i <8; i++)
					elems[e*8+i] = (long) ind[3*local_nodes[i]];
				e++;
			}
			MPI_Barrier(MPI_COMM_WORLD);

			Writer->Write(dset, elems, _grid.GetNrElemGlobal(),_grid.GetNrElem(), 8, _grid.GetElemOffset());
			delete[] elems;
		}
		
        void PrintEmoduli() {
            Writer->Select("/Mesh");
            VectorXd emoduli(_grid.GetNrElem());
			t_index i =0;
			for(_grid.initIterateOverElements(); _grid.TestIterateOverElements(); _grid.IncIterateOverElements()){
				emoduli[i]=_grid.GetElementWeight()*1000; 
				i++;
			}
            MPI_Barrier(MPI_COMM_WORLD);
            Writer->Write("Emoduli", emoduli.data(), _grid.GetNrElemGlobal(),_grid.GetNrElem(), 1, _grid.GetElemOffset());
        }

        //x displacement, res residuum
		void PrintAll(VectorXd &x, VectorXd &res) {
			PrintGrid();

			PostProcess<OctreeGrid<T> > post(_grid);
			VectorXd m, s, eff;
			post.ComputeStressAndStrain(x,m,s,eff);
			MPI_Barrier(MPI_COMM_WORLD);
			
			Writer->Select("/Solution");
			Writer->Write("disp", x.data(), _grid.GetNrNodesGlobal(),_grid.GetNrPrivateNodes(), 3, _grid.GetNodeOffset());
			MPI_Barrier(MPI_COMM_WORLD);
			Writer->Write("VonMises", m.data(), _grid.GetNrElemGlobal(),_grid.GetNrElem(), 1, _grid.GetElemOffset());
			MPI_Barrier(MPI_COMM_WORLD);
			Writer->Write("SED", s.data(), _grid.GetNrElemGlobal(),_grid.GetNrElem(), 1, _grid.GetElemOffset());
			Writer->Write("EFF", eff.data(), _grid.GetNrElemGlobal(),_grid.GetNrElem(), 1, _grid.GetElemOffset());
		}

		void PrintPartition(std::string dset) {
		  VectorXi part(_grid.GetNrElem());
		  part.setConstant(_grid.GetNrElem(), _MyPID);
		  MPI_Barrier(MPI_COMM_WORLD);
		  Writer->Write(dset, part.data(), _grid.GetNrElemGlobal(),_grid.GetNrElem(), 1, _grid.GetElemOffset());
		}




		

		std::string _filename;

		OctreeGrid<T> &_grid;
		
		HDF5_GWriter *Writer;

};
#endif /* HDF5PRINTER_H */
