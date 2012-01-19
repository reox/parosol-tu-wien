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


#ifndef HDF5PARFEPRINTER_H
#define HDF5PARFEPRINTER_H

#include "Config.h"
#include "OctreeGrid.h"
#include "Postprocessing.h"
#include "GWriter.hpp"

#include <eigen2/Eigen/Core>
USING_PART_OF_NAMESPACE_EIGEN

//! A class to write the result in the format that ParFE uses.

/*! HDF5ParfePrinter print the results in the same format as ParFE. In this way,
 *  ParOSol can be used as a fast mesher for ParFE. 
*/


template <class T>
class HDF5ParfePrinter {
	public:
		HDF5ParfePrinter(std::string filename, OctreeGrid<T> &grid):_MyPID(grid.GetPID()), _Size(grid.GetNrCPU()),_filename(filename), _grid(grid)
		{
		  Writer = new HDF5_GWriter(filename, MPI_COMM_WORLD);
		}
		int _MyPID;
		int _Size;
		~HDF5ParfePrinter()
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
			PrintBC();
			PrintParameters();
		}
		
		//!Prints the Parameter that is needed by ParFe
		void PrintParameters() {
		  int n_elements = _grid.GetNrElemGlobal();
		  int n_nodes = _grid.GetNrNodesGlobal();
		  Writer->Select("/Parameters");
		    Writer->Write("element_type", "'hexahedron'");
		    Writer->Write("nr_elements", n_elements);
		    Writer->Write("nr_nodes", n_nodes);
		    Writer->Write("nr_integration_points", (int)8);
		    Writer->Write("nr_dofs_per_node", (int)3);
		    Writer->Write("nr_nodes_per_element", (int)8);
		    Writer->Write("size_of_stress_strain_matrix", (int) 6);
		    Writer->Write("nr_dimensions", (int)3);
		    Writer->Write("aa", (double)0.2);
		    Writer->Write("bb", (double)0.2);
		    Writer->Write("cc", (double)0.2);
		   Writer->Write("nr_material_properties", (int) 2);
		   Writer->Write("nr_material_types", (int)1);

		    double mat_props[2] = {1000,0.3};
		    Writer->Write("materials", mat_props, 1, 1, 2, 0);
		    Writer->Write("iteration_limit", (int)2000);
		    Writer->Write("tolerance", (double) 1e-5);
		    
		}

		//! Print all Boundary Conditions in to the fixed node Data Set.
		void PrintBC() {
		  Writer->Select("/Boundary_conditions");
		  int my_n_fixed_nodes = _grid.bc.FixedNodes_Ind.size();
		  int n_fixed_nodes;
		  MPI_Allreduce(&my_n_fixed_nodes, &n_fixed_nodes, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		  int *sense = new int[my_n_fixed_nodes];
		  double *value = new double[my_n_fixed_nodes];
		  int *node_numbers = new int[my_n_fixed_nodes];
		  
		  t_octree_key offset = _grid.GetNodeOffset();
		  for(int i =0;i<my_n_fixed_nodes;i++) {
		    node_numbers[i] = _grid.bc.FixedNodes_Ind[i]/3 +offset+1;
		    sense[i] = _grid.bc.FixedNodes_Ind[i]%3+1;
		    value[i] = _grid.bc.FixedNodes[i];
		    if (value[i] == 0)
		      value[i] = 1e-16;
		  }
		  
		  int my_fixed_nodes_offset;
		  MPI_Scan(&my_n_fixed_nodes, &my_fixed_nodes_offset,1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
		  my_fixed_nodes_offset -= my_n_fixed_nodes; //include it selfs
		  Writer->Write("fixed_nodes_size", (int)n_fixed_nodes);
		  Writer->Write("fixed_nodes", "Node_number", node_numbers,"Sense", sense, "Value", value, n_fixed_nodes, my_n_fixed_nodes, 1,1,1, my_fixed_nodes_offset );
		  Writer->Write("loaded_nodes_size", (int)0);
		  Writer->Write("restrained_nodes_size", (int)0);
		  delete[] sense;
		  delete[] value;
		  delete[] node_numbers;
		}
		
		  
		void PrintCoord(std::string dset) {
		  Writer->Select("/Mesh");
		  
		  
			t_octree_key k;
			std::vector<OctreeNode> &grid = _grid.GetOctGrid();
			std::vector<OctreeNode>::iterator iter; 
			int x =0,y=0,z=0;
			double res[3];
			_grid.GetRes(res);
			VectorXd coord(_grid.GetNrPrivateNodes()*3);
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
			VectorXd ind(_grid.GetNrDofs());
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
		
        void PrintAll(VectorXd &x, VectorXd &force, VectorXd &res) {

			PrintGrid();

			PostProcess<OctreeGrid<T> > post(_grid);
			VectorXd m, s, eff;
			post.ComputeStressAndStrain(x,m,s,eff);
			MPI_Barrier(MPI_COMM_WORLD);
			
			Writer->Select("/Solution");
			Writer->Write("disp", x.data(), _grid.GetNrNodesGlobal(),_grid.GetNrPrivateNodes(), 3, _grid.GetNodeOffset());
			Writer->Write("VonMises", m.data(), _grid.GetNrElemGlobal(),_grid.GetNrElem(), 1, _grid.GetElemOffset());
			Writer->Write("SED", s.data(), _grid.GetNrElemGlobal(),_grid.GetNrElem(), 1, _grid.GetElemOffset());
		}

		std::string _filename;

		OctreeGrid<T> &_grid;
		
		HDF5_GWriter *Writer;

};
#endif /* HDF5PARFEPRINTER */
