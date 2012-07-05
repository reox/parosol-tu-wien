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

#ifndef VTKPRINTER_H
#define VTKPRINTER_H

#include "Config.h"
#include "OctreeGrid.h"
#include "Postprocessing.h"
#include <fstream>

#define CPU0SAVES(msg) \
			if (_MyPID ==0) { \
				std::ofstream foutcc(file.c_str(), std::ios_base::app);\
				foutcc << msg; \
				foutcc.close();\
			}
/*! Basic VTKPrinter to visualize the mesh in the VTK format.
 *
 * It is used for debung reason. It is not optimized for huge number of cores. 
 */
template <class T>
class VTKPrinter {
	public:
	  	/** 
		 * @brief Constructor
		 * 
		 * @param filename output file
		 * @param grid the finite element grid
		 */
		VTKPrinter(std::string filename, OctreeGrid<T> &grid):_MyPID(grid.GetPID()), _Size(grid.GetNrCPU()),_filename(filename), _grid(grid)
		{
		}

	  	/** 
		 * @brief Writes the Mesh into the VTU file.
		 */
		void PrintGrid() {
			std::string suffix = ".vtu";
			std::string file = _filename + suffix;
			PrintVTKHeader(file);
			PrintCoord(file);
			PrintElems(file);
			PrintVTKFooter(file);
		}
		
		
	  	/** 
		 * @brief Writes the Mesh, solution, Partiotion, and the VonMises stress into the VTU file.
		 */

		void PrintAll(Eigen::VectorXd &x, Eigen::VectorXd &res) {

			PostProcess<OctreeGrid<T> > post(_grid);
            Eigen::VectorXd m, s;
			post.ComputeStressAndStrain(x,m,s);

			std::string suffix = ".vtu";
			std::string file = _filename + suffix;
			PrintVTKHeader(file);
			PrintCoord(file);
			PrintElems(file);

			CPU0SAVES( "<CellData>\n")
			PrintPartition(file);
			PrintCellVector("VonMises",file,  m);
			//PrintCellVector("SED",file, s);
			CPU0SAVES("</CellData>\n");
		
			CPU0SAVES( "<PointData>\n")

			PrintPointVector("solution", file, x);
			//PrintPointVector("residual", file, res);

			CPU0SAVES("</PointData>\n");

			PrintVTKFooter(file);
		}
private:
		 /** 
		 * @brief Print the partions
		 * 
		 * @param file filename
		 */
		void PrintPartition(std::string file) {
			CPU0SAVES("<DataArray type=\"Int32\" Name=\"partition\" NumberOfComponents=\"1\" format=\"ascii\">")
			int nr_elem = _grid.GetNrElem();
			Lock();
			std::ofstream fout(file.c_str(), std::ios_base::app);
			for(int i =0; i < nr_elem; i++) {
				fout << _MyPID << " "; 
			}
			fout.close();
			UnLock();

			CPU0SAVES("</DataArray>\n")
		}

		 /** 
		 * @brief Print coordinates
		 * 
		 * @param file filename
		 */
		void PrintCoord(std::string file) {
			CPU0SAVES("<Points>\n<DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n")
			
			t_octree_key k;
			std::vector<OctreeNode> &grid = _grid.GetOctGrid();
			std::vector<OctreeNode>::iterator iter; 
			int x =0,y=0,z=0;
			double res[3];
			_grid.GetRes(res);

			Lock();
			std::ofstream fout(file.c_str(), std::ios_base::app);
			for(iter = grid.begin(); iter != _grid._GridIteratorEnd; ++iter ) {
				k = iter->key;
				OctKey_to_Coord(k, x, y, z);
				fout << x*res[0] << " " << y*res[1] << " " << z*res[2] << std::endl;
			}
			fout.close();
			UnLock();
			CPU0SAVES("</DataArray>\n</Points>\n")
		}

		 /** 
		 * @brief Print the Elements
		 * 
		 * @param file filename
		 */
		void PrintElems(std::string file) {

			CPU0SAVES("<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">")
			//print element to node
			t_index local_nodes[8];

			int offset =0;
			int tmpoffset =0;
			Lock(offset);
			tmpoffset = offset+_grid.GetNrPrivateNodes();
			UnLock(tmpoffset);

            Eigen::VectorXd ind(_grid.GetNrDofs());
			_grid.Recv_import_Ghost(ind);
			ind.setZero(_grid.GetNrDofs());
			for (unsigned int i = 0; i < _grid.GetNrNodes(); i++) {
				ind[3*i] = i + offset;
			}
			_grid.Send_import_Ghost(ind);
			_grid.Wait_import_Ghost();

			Lock();

			std::ofstream fout(file.c_str(), std::ios_base::app);
			for(_grid.initIterateOverElements(); _grid.TestIterateOverElements(); _grid.IncIterateOverElements()){
				_grid.SearchIndexes(local_nodes); 
				for(int i =0; i <8; i++)
					fout <<(long) ind[3*local_nodes[i]] << " ";
				fout << std::endl;
			}
			fout.close();
			offset += _grid.GetNrPrivateNodes();
			UnLock();
			CPU0SAVES("</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n")

			//print types
			int nr_elem = _grid.GetNrElem();
			Lock();
			fout.open(file.c_str(), std::ios_base::app);
			for(int i=0; i < nr_elem; i++)
				fout << "12 ";
			fout.close();
			UnLock();
			CPU0SAVES(std::endl << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n")

			//print offsets
			int count =8;
			Lock(count);
			fout.open(file.c_str(), std::ios_base::app);
			for(int i=0; i < nr_elem; i++) {
				fout << count << " ";
				count +=8;
			}
			fout.close();
			UnLock(count);

			CPU0SAVES(std::endl << "</DataArray>\n</Cells>\n")
		}
		
		int _MyPID;
		int _Size;
		
		/**
		 * Lock until CPU mypid-1 is finished
		 */
		void Lock() {
			int i =0;
			Lock(i);
		}

		/**
		 * Lock until CPU mypid-1 is finished and receive a message
		 * @param offset message from CPU mypid -1
		 */
		void Lock(int &offset)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Status status;
			if ( 0 != _MyPID) {
				 MPI_Recv(&offset, 1, MPI_INT, _MyPID -1,0, MPI_COMM_WORLD, &status);
			}
		}

		/**
		 * Give the lock free
		 */
		void UnLock()
		{
			int i=0;
			UnLock(i);
		}

		/**
		 * Give the lock free and send a message
		 * 
		 * @param offset message to next cpu
		 */
		void UnLock(int offset) {
			if (_MyPID != _Size -1) {
				MPI_Send(&offset, 1, MPI_INT, _MyPID+1, 0, MPI_COMM_WORLD);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}  
  
  		/**
		 * Converts the key to the three coordinates.
		 * 
		 * @param [in] key
		 * @param [out] x x coodrinate
		 * @param [out] y y coodrinate
		 * @param [out] z z coodrinate
		 */
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
		void PrintVTKHeader(std::string file) {
			t_octree_key nr_elem = _grid.GetNrElemGlobal();
			t_octree_key nr_nodes = _grid.GetNrNodesGlobal();


			std::string head1 = "<?xml version=\"1.0\"?><VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"><UnstructuredGrid>\n<Piece NumberOfPoints=\"";
			std::string head2 = "\" NumberOfCells=\"";
			std::string head3 = "\">";
			CPU0SAVES( head1 << nr_nodes << head2 << nr_elem << head3)
		}

		 /** 
		 * @brief Print a 3 dofs per Vertex
		 * 
		 * @param file filename
		 * @param x Vector that holds the dofs
		 * @param s name of the section
		 */
		void PrintPointVector(std::string s, std::string file, Eigen::VectorXd &x) {
			CPU0SAVES("<DataArray type=\"Float32\" Name=\"" << s << "\" NumberOfComponents=\"3\" format=\"ascii\">")

			int nr_nodes = _grid.GetNrPrivateNodes();
			Lock();
			std::ofstream fout(file.c_str(), std::ios_base::app);
			for(int i =0; i < nr_nodes; i++) {
				fout << x[3*i] << " " << x[3*i+1] << " " << x[3*i+2] << std::endl;
			}
			fout.close();
			UnLock();
			CPU0SAVES("</DataArray>" << std::endl)
		}

		 /** 
		 * @brief Print a 1 information per Cell
		 * 
		 * @param file filename
		 * @param x Vector that holds the dofs
		 * @param s name of the section
		 */
		void PrintCellVector(std::string s, std::string file, Eigen::VectorXd &x) {
			CPU0SAVES("<DataArray type=\"Float32\" Name=\"" << s << "\" NumberOfComponents=\"1\" format=\"ascii\">")

			int nr_elem = _grid.GetNrElem();
			Lock();
			std::ofstream fout(file.c_str(), std::ios_base::app);
			for(int i =0; i < nr_elem; i++) {
				fout << x[i] << " ";
			}
			fout.close();
			UnLock();
			CPU0SAVES("</DataArray>" << std::endl)
		}


		void PrintVTKFooter(std::string file) {
			CPU0SAVES("</Piece>\n</UnstructuredGrid>\n</VTKFile>")
		}
				

		std::string _filename;

		OctreeGrid<T> &_grid;

};
#endif /* TESTTOOLS_H */
