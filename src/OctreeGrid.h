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


#ifndef OCTREEGRID_H
#define OCTREEGRID_H

#include <ostream>

#include <fstream>    //ugly
#include <vector>
#include <set>
#include <map>


#include <utility>
#include <algorithm>

#include <eigen2/Eigen/Core>

#include "ImageReader.h"
#include "BaseGrid.h"
#include "OctreeNode.h"
#include "BoundaryCondition.h"
#include "Config.h"

#include <iterator>

using namespace Eigen;

//Activate for debugging
#define DOUT(msg)
//#define DOUT(msg) 
//	COUT << msg;

//! A class to store a OctreeGridFE mesh.
typedef std::vector<OctreeNode>::iterator t_octree_iterator;
typedef std::vector<OctreeNode> t_octree_vec;

/*! This class stores the grid. It also provides an iterator over the grid. The template parameter
 * defines the type of node in the grid. For a implementation see OctreeNode.
 */
template <typename T>
class OctreeGrid : public BaseGrid
{
	public:
		/** @name Constructors/destructors*/
		//@{
		//!OctreeGrid Constructor
		OctreeGrid() 
		{
			init_Members();
		}

		/** 
		 * @brief Constructor that read in the grid from an ImageeReader
		 * 
		 * @param ir ImageReader
		 * @param elas Array with emoduli
		 */
		OctreeGrid(ImageReader& ir);

		//! OctreeGrid Destructor
		virtual ~OctreeGrid()
		{
			delete[] _sending_status;
			delete[] _sending_request;
			delete[] _remote_request;
			delete[] _remote_status;
		}


		/** @name Some Access function */
		//@{

		t_index GetNrDofs()
		{
			return _nr_nodes*3;
		}

		t_index GetNrNodes()
		{
			return _nr_nodes;
		}

		t_index GetNrPrivateNodes()
		{
			return _GridIteratorEnd - _OctreeGrid.begin();
		}

		t_index GetNrElem()
		{
			return _nr_elem;
		}

		t_octree_key  GetNrElemGlobal()
		{
			return _nr_elem_global;
		}

		t_octree_key  GetNrNodesGlobal()
		{
			return _nr_nodes_global;
		}
		t_octree_key GetNodeOffset()
		{
		  return _my_nodes_offset;
		}
		t_octree_key GetElemOffset()
		{
		  return _my_elem_offset;
		}

		std::vector<OctreeNode>& GetOctGrid()
		{
			return _OctreeGrid;
		}

		void GetLocalDim(t_coord dim[3]) const
		{
			for(int i=0; i< 3; i++)
				dim[i] = this->ldim[i];
		}

		void GetGlobalDim(t_coord dim[3]) const
		{
			for(int i=0; i< 3; i++)
				dim[i] = this->gdim[i];
		}

		void GetRes(double dim[3]) const
		{
			for(int i=0; i< 3; i++)
				dim[i] = this->res[i];
		}

		int GetPID() {
			return _My_PID;
		}
		int GetNrCPU() {
			return _Nr_CPU;
		}
		//@}

		/** @name Iterator on the Grid **/
		//@{

		/** 
		 * @brief Initialize the Iterator on the Grid
		 */
		void initIterateOverElements()
		{
			_GridIterator = _OctreeGrid.begin();
			while((_GridIterator != _OctreeGrid.end()) && (_GridIterator->w <=0))
				++_GridIterator;
		}

		/** 
		 * @brief Tests the termination criteria of the iteration on the grid. 
		 * 
		 * @return True if the end is not reached
		 */
		bool TestIterateOverElements()
		{
			//return _GridIterator != _GridIteratorEnd; 
			return _GridIterator < _GridIteratorEnd; 
		}


		/** 
		 * @brief Increments the Iterator
		 */
		void IncIterateOverElements()
		{
			++_GridIterator;
			while ( (_GridIterator < _GridIteratorEnd) && (_GridIterator->w <=0) ) 
			{
				++_GridIterator;
			}
		}


		/** 
		 * @brief Get the weight of the actual element
		 * 
		 * @return Weight of the element
		 */
		double GetElementWeight()
		{
			return _GridIterator->w;
		}
		//@}

		/** @name Grid and BoundaryCondition methods **/
		//@{

		/** 
		 * @brief Search the 8 index of the neighbours of the actual node. 
         * It stores the indices in the member array called neighbours
		 */
		void SearchIndexes0()
		{
			int ind = (t_index) (_GridIterator - _OctreeGrid.begin());
			t_octree_iterator match = _GridIterator;
			t_octree_key tmp = _GridIterator->key;
			t_octree_key keys[8] = {};
			keys[0] = tmp;
			_neighbours[0] = ind;

			t_octree_key diffx = KeyGen.incX(tmp) ^ tmp;
			t_octree_key diffy = KeyGen.incY(tmp) ^ tmp;
			t_octree_key diffz = KeyGen.incZ(tmp) ^ tmp;

			keys[1] = keys[0]^diffx;
			keys[2] = keys[1]^diffy;
			keys[3] = keys[0]^diffy;
			keys[4] = keys[0]^diffz;
			keys[5] = keys[1]^diffz;
			keys[6] = keys[2]^diffz;
			keys[7] = keys[3]^diffz;
			//indexeoffsets for the points
			//
			//    7------6
			//   /|     /|
			//  4------5 | 
			//  | |    | |
			//  | 3----|-2
			//  |/     |/
			//  0------1

			for(int i =1; i<8;i++) {
				match = _GridIterator+1;
				ind = GetIndexOfKey(match,keys[i] );
				ind = GetIndexOfKey(keys[i] );
				if (ind < 0) {
					COUT << " ERRRROOOOOORRRRRRRRR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! on PID: " << _My_PID << std::endl;
					COUT << " Could not find key: " << keys[i] << std::endl;
					COUT << " main key is: " << keys[0] << " weight is " << _OctreeGrid[_neighbours[0]].w << std::endl;
				}
				_neighbours[i] = ind;
			}
			return;
		}


		/** 
		 * @brief Get the global keys of the nodes of the element
		 * 
		 * @param keys  keys of the nodes
		 * @param key0 global key of the element
		 */
		inline	void ComputeKeysofNodes(t_octree_key keys[8], t_octree_key key0) {
			//indexeoffsets for the points
			//
			//    7------6
			//   /|     /|
			//  4------5 | 
			//  | |    | |
			//  | 3----|-2
			//  |/     |/
			//  0------1
			keys[0] = key0;
			t_octree_key diffx = KeyGen.incX(key0) ^ key0;

			t_octree_key diffy = KeyGen.incY(key0) ^ key0;
			t_octree_key diffz = KeyGen.incZ(key0) ^ key0;
			keys[1] = keys[0]^diffx;
			keys[2] = keys[1]^diffy;
			keys[3] = keys[0]^diffy;
			keys[4] = keys[0]^diffz;
			keys[5] = keys[1]^diffz;
			keys[6] = keys[2]^diffz;
			keys[7] = keys[3]^diffz;
		}


		void SearchIndexes(t_index *neigh) {
			SearchIndexes();
			for(int i=0; i<8; i++)
				neigh[i] = _neighbours[i];
		}

		/** 
		 * @brief Search the 8 index of the neighbours and stores in the neighbour array member.
		 * This is the optimized version.
		 */
		void SearchIndexes()
		{
			t_index ind = (t_index) (_GridIterator - _OctreeGrid.begin());
			t_octree_iterator tmpx, tmpy, tmpz, tmpxy, tmp;
			t_octree_key key0 = _GridIterator->key;
			t_octree_key keys[8];
			_neighbours[0] = ind; 

			ComputeKeysofNodes(keys, key0);



			tmpx = _GridIterator+1;
			DOUT("key 1\n")
			ind = GetIndexOfKeyA(tmpx,keys[1] );
			_neighbours[1] = ind; 

			tmpy = _GridIterator+1;
			DOUT("key 3\n")
			ind = GetIndexOfKeyA(tmpy, keys[3] );
			_neighbours[3] = ind; 

			tmpz = _GridIterator+1;
			DOUT("key 4\n")
			ind = GetIndexOfKeyA(tmpz, keys[4] );
			_neighbours[4] = ind; 

			if (keys[3] > keys[1])
				tmpxy = tmpy+1;
			else
				tmpxy = tmpx+1;
			DOUT("key 2\n")
			ind = GetIndexOfKeyA(tmpxy, keys[2] );
			_neighbours[2] = ind; 




			if (keys[4] > keys[1])
				tmp = tmpz+1;
			tmp = tmpx+1;

			DOUT("key 5\n")
				ind = GetIndexOfKeyA(tmp, keys[5] );
			_neighbours[5] = ind; 

			if (keys[5] > keys[2])
				tmp = tmp+1;
			else
				tmp = tmpxy+1;

			DOUT("key 6\n")
				ind = GetIndexOfKeyA(tmp, keys[6] );
			_neighbours[6] = ind; 

			if (keys[4] > keys[3])
				tmp = tmpz+1;
			else
				tmp = tmpy+1;
			DOUT("key 7\n");
			ind = GetIndexOfKeyA(tmp, keys[7] );
			_neighbours[7] = ind; 

		}

		/** 
		 * @brief Copies the 8x3 nodal values of the actual element in an fixed size vector
		 * The iterator must be fixed on the actual element.
		 * 
		 * @param x Source vector
		 * @param pref fixed size vector
		 */
		void GetNodalDisplacementsOfElement(VectorXd &x, Matrix<double, 24, 1> &pref)
		{
			SearchIndexes();

			for(int i =0; i<8;i++)  {
				Copy3(x, _neighbours[i], pref, i);
			}

			return;
		}



		/** 
		 * @brief Updates the 8x3 nodal values of the actual element in the vector.
		 * The iterator must be fixed on the actual element.
		 *
		 * @param x Vector that gets NrCPUed
		 * @param pref fixed size vector with 24 values that are summed into the vector
		 * @param factor x+= factor * pref
		 */
		void SumInToNodalDisplacementsOfElement(VectorXd &x, Matrix<double, 24, 1> &pref, double factor)
		{
			for(int i =0; i<8;i++) {
				for(t_index j=0; j < 3; ++j) 
					x[3*_neighbours[i] +j] +=factor * pref[3*i+j];
			}
		}

		// with returns -1 if is not in the grid 

		/** 
		 * @brief Search the Index of the key in the Octree with binary search.
		 * It ist the slowest variant.
		 * 
		 * @param key Global key of the node
		 * 
		 * @return if key is found the index is returned else -1
		 */
		inline int GetIndexOfKey(t_octree_key key)
		{
			t_octree_iterator tmp = lower_bound(_OctreeGrid.begin(), _OctreeGrid.end(), key);
			if (tmp != _OctreeGrid.end() && !(key < *tmp))
				return tmp - _OctreeGrid.begin();
			else
				return -1;
		}

		/** 
		 * @brief Search the Index of the key in the Octree with binary search
         * in the array above match. The key must exist in the vector.
		 * 
		 * @param match Starting iterator for the search. After the call the iterator points to the found object.
		 * @param key Global key
		 * 
		 * @return returns the index of the searches key
		 */
		inline t_index GetIndexOfKey(t_octree_iterator &match, t_octree_key key)
		{
			t_octree_iterator tmp = match;
			match = lower_bound(match, _OctreeGrid.end(), key);
			return (t_index) (match - _OctreeGrid.begin());
		}

		/** 
		 * @brief Search the Index of the key in the Octree with extend/binary
         * search in the array above the  match iterator. The key must exist in the vector.
		 * 
		 * @param match Starting iterator for the search. After the call the iterator points to the found object.
		 * @param key Global key
		 * 
		 * @return returns the index of the searches key
		 */
		inline t_index GetIndexOfKeyA(t_octree_iterator &match, t_octree_key key)
		{
			if (key == *match)
				return (t_index) (match - _OctreeGrid.begin());

			++match;

			t_octree_iterator tmp = match;
			t_octree_iterator end = _OctreeGrid.end();
			int i = 1, j = 8;
			int count = _OctreeGrid.end() - match;
			//cout << "Intervall " << end-match << endl;
			while (count > i) {
				if (key < *(match+i)) {
					end = match+i;
					break;
				}
				match = match+i;
				count = count -i;

				i=i *j;
			} 

			tmp = lower_bound(match, end, key);
			//assert(tmp->key== key);
			match = tmp;
			return (t_index) (match - _OctreeGrid.begin());
		}


		/*
		 * @brief Computes the index of the Boundary Condition.
		 * 
		 * @param[out] list index of the boundary condition nodes
		 */
		void GetIndexOfBC(std::vector<t_boundary_node> &list) {
			std::vector<  t_boundary_node >::iterator iter;
			t_index index;
			t_octree_key key;

			for (iter = _BC_list.begin(); iter !=_BC_list.end(); ++iter) {
				key = iter->first;
				t_octree_iterator tmp = lower_bound(_OctreeGrid.begin(), _GridIteratorEnd, key); //hmhmhm
				if (tmp != _GridIteratorEnd && !(key < *tmp))
					index = tmp - _OctreeGrid.begin();
				else {
				 index = -1; COUT << "Error! index of BC not found!\n";
				}
				list.push_back(t_boundary_node(index, iter->second));
			}
		}

		void GenerateBC() {
			std::vector<t_boundary_node> list;
			Distribute_the_BC();
			GetIndexOfBC(list);
			CheckIndexOfBC(list);
			bc.GenerateBC(list);
			Delete_BC_vector();
			return;
		}

		void Distribute_the_BC()
		{
			std::vector<int> send_count(_Nr_CPU);
			std::vector<int> send_pos(_Nr_CPU);
			std::vector< std::vector< t_boundary_node >::iterator> send_pos_iters(_Nr_CPU+1);
			std::vector<int> recv_count(_Nr_CPU);
			std::vector<int> recv_pos(_Nr_CPU); 

			std::vector< t_boundary_node > new_bc_list;


			std::vector<  t_boundary_node >::iterator iter = _BC_list.begin();
			send_pos_iters[0] = _BC_list.begin();
			send_pos_iters[_Nr_CPU] = _BC_list.end();

			for(int i = 0; i < _Nr_CPU; i++)
			{
				while( (iter != _BC_list.end()) && (iter->first < _smallest_elements[i+1]))
					iter++;
				send_pos_iters[i+1] = iter;
			}

			for(int i = 0; i < _Nr_CPU; i++) {
				send_count[i] = send_pos_iters[i+1] - send_pos_iters[i];
			}
			for(int i = 0; i < _Nr_CPU; i++) {
				send_count[i] *= sizeof(t_boundary_node);
			}

			
			MPI_Alltoall (&send_count[0], 1, MPI_INT, &recv_count[0], 1, MPI_INT, MPI_COMM_WORLD);

			send_pos[0] = 0;
			recv_pos[0] = 0;
			for(int i=1; i < _Nr_CPU; i++) {
				send_pos[i] = (send_pos[i-1] + send_count[i-1]);
				recv_pos[i] = recv_pos[i-1] + recv_count[i-1];
			}
			
			int newsize = (recv_pos[_Nr_CPU-1] + recv_count[_Nr_CPU-1])/sizeof(t_boundary_node);
			new_bc_list.resize(newsize);

			void * send_addr;
			void * recv_addr;
			if (_BC_list.size() > 0)
				send_addr = &_BC_list[0];
			else
				send_addr = NULL;

			if (new_bc_list.size() > 0)
				recv_addr = &new_bc_list[0];
			else
				recv_addr = NULL;

			MPI_Alltoallv(send_addr, &send_count[0], &send_pos[0], MPI_CHAR,
						  recv_addr, &recv_count[0], &recv_pos[0], MPI_CHAR,
						  MPI_COMM_WORLD);
			std::set< t_boundary_node > bc_set;
			for(unsigned r =0; r < new_bc_list.size(); r++)
				bc_set.insert(new_bc_list[r]);
			Delete_BC_vector(new_bc_list);
			Delete_BC_vector();
			_BC_list.reserve(bc_set.size());
			_BC_list.insert(_BC_list.begin(), bc_set.begin(), bc_set.end());
			

		}

		/** 
		 * @brief Deletes the temporary storage of the boundary condition
		 */
		void Delete_BC_vector() {
			std::vector< t_boundary_node > tmp;
			_BC_list.swap(tmp);
		}
		void Delete_BC_vector(std::vector< t_boundary_node > &x) {
			std::vector< t_boundary_node> tmp;
			x.swap(tmp);
		}

		//!@}

		friend std::ostream& operator<<(std::ostream& stream, const OctreeGrid &grid) {
			grid.print(stream);
			return stream;
		}

		//! @name Debug function
		//!@{
		void CheckGrid(const char* msg)
		{
			t_octree_iterator iter;
			for(iter = _OctreeGrid.begin(); iter < _OctreeGrid.end(); ++iter) {
				if (iter->key == 0) {
					std::cout << "Error after " << msg << "  on PID: " << _My_PID << std::endl;
					std::cout.flush();
				}
			}

		}
		void PrintVector(VectorXd &x) {
			unsigned numberelem = _GridIteratorEnd - _OctreeGrid.begin();
			numberelem *=3;
			MPI_Barrier(MPI_COMM_WORLD);

			for(int j =0; j < _Nr_CPU; j++) {
				if (_My_PID == j) {
					COUT << "CPU "<< j <<" holds " << numberelem << " elements" << std::endl; 
				}
				MPI_Barrier(MPI_COMM_WORLD);
			}

			for(int j =0; j < _Nr_CPU; j++) {
				if (_My_PID == j) {
					for(unsigned i=0; i < numberelem; i++)
							COUT << x[i] << std::endl; 
					COUT << std::endl; 
						}
				MPI_Barrier(MPI_COMM_WORLD);
			}
		}

		void Print();  

		void PrintGrid(const char *filename)
		{
			int i=0;
			int PID = _My_PID;
			int SIZE = _Nr_CPU;
			MPI_Status status;
			MPI_Barrier(MPI_COMM_WORLD);

			if ( 0 != PID) { 
				MPI_Recv(&i, 1, MPI_INT, PID -1,0, MPI_COMM_WORLD, &status); 
			} 
			std::fstream fb(filename, std::fstream::out| std::fstream::app);
			t_octree_iterator Iter;
			fb << "On Proc:" << PID << std::endl;
			for(Iter = _OctreeGrid.begin(); Iter != _GridIteratorEnd; Iter++)
			{
					fb << Iter->key << " " << Iter->w << std::endl;
			}
			fb.flush();
			fb.close();

			if (PID != SIZE -1) { 
				MPI_Send(&i, 1, MPI_INT, PID+1, 0, MPI_COMM_WORLD);  
			} 
			MPI_Barrier(MPI_COMM_WORLD); 
		}

		void PrintBC(const char *filename)
		{
			int i=0;
			int PID = _My_PID;
			int SIZE = _Nr_CPU;
			MPI_Status status;
			MPI_Barrier(MPI_COMM_WORLD);

			if ( 0 != PID) { 
				MPI_Recv(&i, 1, MPI_INT, PID -1,0, MPI_COMM_WORLD, &status); 
			} 
			std::fstream fb(filename, std::fstream::out| std::fstream::app);
			fb << "On Proc:" << PID << std::endl;
			std::vector<t_index>::iterator bc_it;
			int j =0;
			for(bc_it = bc.FixedNodes_Ind.begin(); bc_it != bc.FixedNodes_Ind.end(); ++bc_it)
			{
				//if (Iter->w > 0) {
					fb << _OctreeGrid[*bc_it/3] << " " << bc.FixedNodes[j++] << std::endl;
				//}
			}
			fb.flush();
			fb.close();

			if (PID != SIZE -1) { 
				MPI_Send(&i, 1, MPI_INT, PID+1, 0, MPI_COMM_WORLD);  
			} 
			MPI_Barrier(MPI_COMM_WORLD); 
		}




		bool CheckIndexOfBC(std::vector<t_boundary_node> &list) {
			std::vector< t_boundary_node>::iterator BC_iter = _BC_list.begin();
			std::vector<t_boundary_node>::iterator list_iter = list.begin();
			bool ret = true;

			for(;list_iter != list.end(); ++list_iter, ++BC_iter) {
				if (_OctreeGrid[list_iter->first] != BC_iter->first) {
					COUT << "BC key: " << BC_iter->first << " OCT Key " << _OctreeGrid[list_iter->first] << std::endl;
					ret = false;
				}
			}
			return ret;
		}
		//!@}
		
		BoundaryCondition bc;

		/** @name Method for communication **/
		//@{

		/** 
		 * @brief Prints the amount of communication in a matrix form.
		 */
		void Print_Communication_Matrix(std::ostream &stream) {
			int nr_cpu = _remote_elements_pid.size()+_sending_elements_pid.size();
			int total;
			MPI_Reduce(&nr_cpu, &total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
			if (_My_PID == 0) {
			stream << "\n\n\n\n%%MatrixMarket matrix coordinate real general\n";
			stream << _Nr_CPU << " " << _Nr_CPU << " " << total << std::endl;
			}
				for(int j =0; j < _Nr_CPU; j++) {
					if (_My_PID == j) {
						for(unsigned r =0; r < _sending_elements_pid.size(); r++) {
							stream << j+1 << "\t" << _sending_elements_pid[r]+1 << "\t" << _sending_elements_pos[r+1] -_sending_elements_pos[r] << std::endl;
						}
						for(unsigned r =0; r < _remote_elements_pid.size(); r++) {
							stream << j+1 << "\t" << _remote_elements_pid[r]+1 << "\t" << _remote_elements_pos[r+1] -_remote_elements_pos[r] << std::endl;
					}

					}
					MPI_Barrier(MPI_COMM_WORLD);
				}

		}

		/** 
		 * @brief Waits of all cpu in the import process
		 */
		void Wait_import_Ghost() {
			MPI_Waitall(_remote_elements_pid.size(), _remote_request, _remote_status);
			MPI_Waitall(_sending_elements_pid.size(), _sending_request, _sending_status);
		}

		/** 
		 * @brief Sends the data of the import process
		 */
		void Send_import_Ghost(VectorXd &disp)
		{
			t_index ind;
			for(unsigned i=0;i < _remote_elements_index.size(); i++)
			{
				ind = _remote_elements_index[i];
				_remote_elements_buf[3*i+0] = disp[3*ind+0];
				_remote_elements_buf[3*i+1] = disp[3*ind+1];
				_remote_elements_buf[3*i+2] = disp[3*ind+2];
			}

			for(unsigned i=0; i < _remote_elements_pid.size(); i++) {
				double* buff = &_remote_elements_buf[3*_remote_elements_pos[i]];
				int size = _remote_elements_pos[i+1] - _remote_elements_pos[i];
				MPI_Isend(buff, 3*size, MPI_DOUBLE, _remote_elements_pid[i], 99, MPI_COMM_WORLD, &_remote_request[i]);
			}
			return;
		}

		/** 
		 * @brief Post the receives of the import process
		 */
		void Recv_import_Ghost(VectorXd &disp)
		{
			int offset = _GridIteratorEnd - _OctreeGrid.begin();
			for(unsigned i =0;i < _sending_elements_pid.size(); i++) {
				double *recvbuf = &disp[3*(_sending_elements_pos[i]+offset)];
				int size = _sending_elements_pos[i+1] - _sending_elements_pos[i];
				MPI_Irecv(recvbuf, 3*size, MPI_DOUBLE, _sending_elements_pid[i], 99, MPI_COMM_WORLD, &_sending_request[i]);
			}
			return;
		}

		/** 
		 * @brief Sends the data of the export process
		 */
		void Send_export_Ghost(VectorXd &disp)
		{
			int offset = _GridIteratorEnd-_OctreeGrid.begin();
			for(unsigned i =0;i < _sending_elements_pid.size(); i++) {
				double *buf = &disp[3*(_sending_elements_pos[i]+offset)];
				int size = _sending_elements_pos[i+1] - _sending_elements_pos[i];
				MPI_Isend(buf, 3*size, MPI_DOUBLE, _sending_elements_pid[i], 10, MPI_COMM_WORLD, &_sending_request[i]);
			}
			return;
		}

		/** 
		 * @brief Post the receives of the export process
		 */
		void Recv_export_Ghost()
		{
			for(unsigned i=0; i < _remote_elements_pid.size(); i++) {
				double* buff = &_remote_elements_buf[3*_remote_elements_pos[i]];
				int size = _remote_elements_pos[i+1] - _remote_elements_pos[i];
				MPI_Irecv(buff, 3*size, MPI_DOUBLE, _remote_elements_pid[i], 10, MPI_COMM_WORLD, &_remote_request[i]);
			}
			return;
		}

		/** 
		 * @brief Waits of all cpu in the export process and sum it into the local nodes
		 */
		void WaitAndCopy_export_Ghost(VectorXd &disp) {
			MPI_Waitall(_remote_elements_pid.size(), _remote_request, _remote_status);
			MPI_Waitall(_sending_elements_pid.size(), _sending_request, _sending_status);
			t_index ind;
			for(unsigned i=0;i < _remote_elements_index.size(); i++)
			{
				ind = _remote_elements_index[i];
				//cout << "copying " << ind << endl;
				disp[3*ind+0] += _remote_elements_buf[3*i+0]; 
				disp[3*ind+1] += _remote_elements_buf[3*i+1];
				disp[3*ind+2] += _remote_elements_buf[3*i+2];
			}
		}


		/** 
		 * @brief Computes the dot product of two vectors that correspond to the grid. 
		 */
		double dot(VectorXd &vx, VectorXd &vy)
		{
			double localdot =0;
			double globaldot =0;
			t_index length = 3*(_GridIteratorEnd - _OctreeGrid.begin());

#define BLOCK 8192
            //compute the dot product block wise
			const unsigned int blocksize = BLOCK; //1 << 14; //2 ^ 15
			int start = 0;
			while (length > blocksize) {
				localdot += vx.segment<BLOCK>(start).dot(vy.segment<BLOCK>(start));
				start += blocksize;
				length -= blocksize;
			}
			if (length !=0)
				localdot += vx.segment(start,length).dot(vy.segment(start,length));

			MPI_Allreduce(&localdot, &globaldot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			return globaldot;
		}
        //!@}


//	protected:
	public:
		T KeyGen;
		void init_Members();
		t_index _nr_elem;
		t_octree_key _nr_elem_global;
		t_octree_key _my_elem_offset;
		t_octree_key _nr_nodes_global;
		t_octree_key _my_nodes_offset;
		t_index _nr_nodes;
		
		t_octree_key _neighbours[8];

		std::vector<OctreeNode> _OctreeGrid;
		t_octree_iterator _GridIterator;
		t_octree_iterator _GridIteratorEnd;
		


		std::vector< t_boundary_node> _BC_list;

		//! Convert an image to a grid
		void GenerateGrid(double elas[]);
		//! Convert an gird to a Octree
		void GenerateOctree();
		//some print function
		std::ostream& print(std::ostream& stream) const;


		void Compute_min_max(t_octree_key min, t_octree_key max) {
			t_octree_key min_max_in[2];
			t_octree_key min_max_out[2];
			min_max_in[0] = -min;
			min_max_in[1] = max;

			MPI_Allreduce(min_max_in, min_max_out, 2, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
			min_max_out[0] = -min_max_out[0];
			_min_max[0] = min_max_out[0];
			_min_max[1] = min_max_out[1]+1;
			return;
	}
		/** @name Parallel Setup */
		//@{
		/** 
		 * @brief Refine the buckets with more than the target size by an
         * factor. It goes only one time through all buckets.
		 * 
		 * @param[in, out] bucketsboarder Holds the lower and the upper key
         * of the bucket. First bucket bucketboarder[0] - bucketsboarder[1],
         * second bucket bucketboarder[1] - bucketboarder[2], (it on all proc the same)
		 * @param[in, out] buckets start and end iterator of the buckets.
         * Stored in the same schema as the bucket boarders.
		 * @param[out] bucket_sizes The sizes of the target buckets.
         * @param[in] max_size maximum size of the buckets.
         * @return true if it was refined by the factor
		 */
bool Refine_Buckets(std::set<t_octree_key> &bucketsboarder, std::vector<t_octree_iterator> &buckets, std::vector<t_index> &bucket_sizes, t_index max_size)
		{
			const int refinement_factor = 2;
			bool res =false;
			Compute_Nr_Element_in_buckets(buckets, bucket_sizes);

			std::set<t_octree_key>::iterator elem, old_elem;

			double new_bucket_size;
			t_octree_key new_boundary;

			elem = bucketsboarder.begin();
			old_elem = elem;
			++elem;
			int i=0;
		
			while(elem != bucketsboarder.end()) 
			{
				if (max_size < bucket_sizes[i]) {
					res = true;
					new_bucket_size = ((double)(*elem - *old_elem))/refinement_factor;
					for(int i = 1; i < refinement_factor; i++) {
						new_boundary = *old_elem + new_bucket_size * i;
						bucketsboarder.insert(new_boundary);
					}
				}
				old_elem = elem;
				++elem;
				++i;
			}
			if (res)
				Compute_Start_and_End_of_the_Buckets(bucketsboarder, buckets);

			return res;
		}

		/** 
		 * @brief Computes an even distribution over localsizes.size() items.
         *
		 * @param[in] globalsize Global number that is divided into chunks
		 * @param[in, out] localsizes vector of t_index that will be filled
		 */
		void Compute_equally_distribution(t_octree_key globalsize, std::vector<t_index> &localsizes) {
			t_index nr_buckets = localsizes.size();
			double average_nr_elem = ((double) globalsize)/nr_buckets;
			for(t_index i = 0; i < nr_buckets; i++) {
				if (i == nr_buckets) {
					localsizes[i] = globalsize - (t_octree_key)(average_nr_elem*i);
				} else {
					localsizes[i] = (t_octree_key)(average_nr_elem * (i+1)) - (t_octree_key)(average_nr_elem * i);
				}
			}
			return;
		}


		/** 
		 * @brief Fill in the target Buckets. It tries to achieve equal size per Bucket.
		 * 
		 * @param[in] histogram The buckets are stored with iterators for the left and right boundary
		 * @param[in] histogram_sizes The sizes of the buckets in the histogram.
		 * @param[out] buckets The target buckets. They are stored in the same format as the histogram.
		 * @param[in,out] bucket_sizes The sizes of the target buckets. As input it holds the sizes to achieve. Finally the actual size of the buckets.
		 */
		void Fill_in_the_Buckets(std::vector<t_octree_iterator> &histogram, std::vector<t_index> &histogram_sizes, std::vector<t_octree_iterator> &buckets, std::vector<t_index> &bucket_sizes)
		{
			int cur_bucket_size;
			int h=0;
			int nr_buckets_in_histogram = histogram_sizes.size();
			int diff =0;
			int nr_target_buckets = bucket_sizes.size();

			buckets[0] = histogram[0];

			for(int i =0; i < nr_target_buckets; i++) {
				cur_bucket_size = 0;
				while (
						(h != nr_buckets_in_histogram) && 
						(cur_bucket_size + histogram_sizes[h]  <= bucket_sizes[i] - diff)) 
				{
					cur_bucket_size += histogram_sizes[h];
					h++;
				}
				//upper boarder
				buckets[i+1] = histogram[h];
				diff += cur_bucket_size - bucket_sizes[i];
				bucket_sizes[i] = cur_bucket_size; //fill up the sizes
			}

			while(h < nr_buckets_in_histogram)
			{
				bucket_sizes[nr_target_buckets-1] += bucket_sizes[h];
				h++;
			}
			//the last cpubuckets end on the last bucket;
			buckets[nr_target_buckets] = histogram.back();
			return;
		}

		/** 
		 * @brief Compute the global number in buckets. If nodes a double on
         * the processor then they are counted twice.
         * 
         * @param[in] buckets Holds the start end end iterator of the buckets
         * @param[out] bucketsize Give the global size back
		 */
		void Compute_Nr_Element_in_buckets(std::vector<t_octree_iterator> &buckets, std::vector<t_index> &bucketsize)
		{
			t_index number_bucket = buckets.size()-1;

			std::vector< t_index > bucketsize_in(number_bucket);
			bucketsize.resize(number_bucket);

			for(unsigned i=0; i < number_bucket; i++) {
				bucketsize_in[i] = buckets[i+1] - buckets[i]; 
			}
			MPI_Allreduce(&bucketsize_in[0], &bucketsize[0], number_bucket, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
			return;
		}

		/** 
		 * @brief Computes the iterator to the gird array from the buckets
         * boarder that are given by the global keys. 
         * 
         * @param[in] bucketsboarder Holds the global keys of the buckets
         * boarder
         * @param[out] buckets Stores the iterator of the boarders. 
		 */
		void Compute_Start_and_End_of_the_Buckets(std::set<t_octree_key> &bucketsboarder, std::vector<t_octree_iterator> &buckets)
		{
			int number_bucket = bucketsboarder.size() -1;
			std::set<t_octree_key>::iterator elem;
			elem = bucketsboarder.begin();
			++elem;


			buckets.resize(number_bucket+1);
			buckets[0] = _OctreeGrid.begin();
			buckets[number_bucket] = _OctreeGrid.end();
			for(int i = 1; i < number_bucket; i++) {
				buckets[i] = lower_bound(buckets[i-1], buckets[number_bucket], *elem);
				++elem;
			}

			return;
		}

		/** 
		 * @brief Sorting the nodes according the global key. CPU 0 holds the
         * smallest CPU n the biggest
         */
		void Sort() {
			int Nr_CPU;
			MPI_Comm_size (MPI_COMM_WORLD, &Nr_CPU);
			int MyPID;
		    MPI_Comm_rank (MPI_COMM_WORLD, &MyPID);
				
			//Compute globalsize, Maximum and Minimum
			t_octree_key localsize = _OctreeGrid.size();
			t_octree_key globalsize;
			MPI_Allreduce(&localsize, &globalsize, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

			t_octree_key min = (((t_octree_key) 1)<< 50);
			t_octree_key max = -1; 
			if (_OctreeGrid.size() > 0) {
				min = _OctreeGrid.front().key;
				max = _OctreeGrid.back().key;
			}
			Compute_min_max(min, max);

			//Create equally sized buckets. Boarders are stored with values in the keyspace
			//The buckets are stored with an iterator on the octree
			//bucket 0 goes from buckets[0] up to (without) buckets[1]
			// Initial bucket is all
			std::set<t_octree_key> bucketsboarder;
			bucketsboarder.insert(_min_max[0]); //insert the minimum
			bucketsboarder.insert(_min_max[1]); //insert the maximum (boarder is allways without number
			
			std::vector<t_octree_iterator> buckets(2);
			buckets[0]=_OctreeGrid.begin();
			buckets[1]=_OctreeGrid.end();

			std::vector<t_index> bucket_sizes;

			double average_nr_elem_per_cpu =  ((double) globalsize)/Nr_CPU;
			double max_nr_elem_per_bucket =  (((double) globalsize)/Nr_CPU/32);
			if (max_nr_elem_per_bucket < 4)
				max_nr_elem_per_bucket = 4.1;

			int count_refinement_steps =0;

			while (Refine_Buckets(bucketsboarder, buckets, bucket_sizes, max_nr_elem_per_bucket)) {
				count_refinement_steps++;
			}

			if (MyPID == 0) {
				COUT << "Refinement steps: " << count_refinement_steps << std::endl;
				COUT << "Number of Buckets: " << bucket_sizes.size() << std::endl;
			}

			std::vector<t_index> CPU_buckets(Nr_CPU);
			Compute_equally_distribution(globalsize, CPU_buckets);

			if (MyPID ==0) {
				COUT << "Average eleme per cpu: " << average_nr_elem_per_cpu << std::endl;
			}

			std::vector<t_octree_iterator> cpu_buckets(Nr_CPU+1);
			Fill_in_the_Buckets(buckets, bucket_sizes, cpu_buckets, CPU_buckets); 

			//Send buckets to the destination processor
			OctreeNode *sending_addr;
			t_index sending_count;
			std::vector<OctreeNode> newgrid(CPU_buckets[MyPID]);

			std::vector<t_index> recv_count(Nr_CPU);
			std::vector<t_index> recv_pos(Nr_CPU);

			for(int recv_cpu =0; recv_cpu < Nr_CPU; recv_cpu++) {
				if (cpu_buckets[recv_cpu] == _OctreeGrid.end())
					sending_addr = NULL;
				else
					sending_addr = &*cpu_buckets[recv_cpu];

				sending_count = cpu_buckets[recv_cpu+1]-cpu_buckets[recv_cpu];
				sending_count *= sizeof(OctreeNode);
				MPI_Gather(&sending_count,1, MPI_INT, &recv_count[0], 1, MPI_INT, recv_cpu, MPI_COMM_WORLD);

				if (MyPID == recv_cpu) {
					recv_pos[0] = 0;
					for (int i =1; i < Nr_CPU; i++)
						recv_pos[i] = recv_pos[i-1] + recv_count[i-1];
				}
				

				MPI_Gatherv(sending_addr, sending_count, MPI_BYTE,
							&newgrid[0], (int*) &recv_count[0], (int*) &recv_pos[0], MPI_BYTE,
							recv_cpu, MPI_COMM_WORLD);
			}

			//sort
			Delete_a_vector_Octreenodes(_OctreeGrid);
			sort(newgrid.begin(), newgrid.end());
			

			//Delete Duplicates
			int i_dubfree =0;
			unsigned i=0;

			while (i < newgrid.size() && newgrid[i].w <=0)
				i++;
			if (i >0) {
				newgrid[0] = newgrid[i];
			}

			
			for(; i < newgrid.size(); i++) {
				if (newgrid[i_dubfree] < newgrid[i]) {
					i_dubfree++;  //go to the next elem
					newgrid[i_dubfree] = newgrid[i];
				} else {
					if (newgrid[i].w > 0) {
						newgrid[i_dubfree] = newgrid[i];
					}
				}
			}
			i_dubfree++;
			newgrid.resize(i_dubfree);


			_OctreeGrid.resize(i_dubfree);
			std::copy(newgrid.begin(), newgrid.end(), _OctreeGrid.begin()); 
			Delete_a_vector_Octreenodes(newgrid);

	
			return;
		}

		/** 
		 * @brief Add missing local nodes that belongs to other processors.
         */
		void Add_Missing_Nodes() {
			t_octree_iterator iter;
			std::set<OctreeNode>::iterator it;
			std::set<OctreeNode>::iterator it_inserted;
			std::set<OctreeNode> boundarynodes;
			t_octree_key nodekeys[8];
			OctreeNode tmp_node;

			for( iter = _OctreeGrid.begin(); iter != _OctreeGrid.end(); iter++) {
				if (iter->w >0) {
					ComputeKeysofNodes(nodekeys, iter->key);
					for(int i=0; i<8; i++) {
						if (GetIndexOfKey(nodekeys[i]) < 0) {
							tmp_node = OctreeNode(nodekeys[i], -2.0);
							boundarynodes.insert(tmp_node);
						}
					}
				}
			}

			//finally append it
			_nr_nodes = _OctreeGrid.size();
			_OctreeGrid.resize(_nr_nodes+boundarynodes.size());
			std::copy(boundarynodes.begin(), boundarynodes.end(), _OctreeGrid.begin()+_nr_nodes);

			_nr_nodes = _OctreeGrid.size();
			std::vector<OctreeNode> tmp(_nr_nodes);
			std::copy(_OctreeGrid.begin(), _OctreeGrid.end(), tmp.begin());
			_OctreeGrid.swap(tmp);

			return;
		}

		void Update_OctreesHelper() {
			_GridIteratorEnd = lower_bound(_OctreeGrid.begin(), _OctreeGrid.end(), _smallest_elements[_My_PID+1]);
			_nr_nodes = _OctreeGrid.size();
			_nr_elem = 0;
			t_octree_iterator iter;
			for(iter = _OctreeGrid.begin(); iter < _GridIteratorEnd; ++iter) {
				if (iter->w > 0)
					_nr_elem++;
			}
			t_octree_key local_elems = _nr_elem;
			MPI_Allreduce(&local_elems, &_nr_elem_global, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
			MPI_Scan(&local_elems, &_my_elem_offset,1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
			_my_elem_offset -= local_elems; //include it selfs
			local_elems = _GridIteratorEnd - _OctreeGrid.begin();
			MPI_Allreduce(&local_elems, &_nr_nodes_global, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
			MPI_Scan(&local_elems, &_my_nodes_offset,1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
			_my_nodes_offset -= local_elems; //include it selfs
		}

		/** 
		 * @brief Compute the n smallest element number on each cpu and stores
         * locally in the smallest_elem vector. 
         *  _min_max have to be correctly set
         */
		void Compute_Smallest_Elementnumber() {
			//find the smallest number of each cpu
			t_octree_key my_smallest_elem = -1;
			if (_OctreeGrid.size() > 0)
				my_smallest_elem = _OctreeGrid.begin()->key;
			_smallest_elements.resize(_Nr_CPU+1);
			MPI_Allgather(&my_smallest_elem, 1, MPI_LONG,
					&_smallest_elements[0], 1, MPI_LONG, MPI_COMM_WORLD);
			_smallest_elements[_Nr_CPU] = _min_max[1];
			for(unsigned i =0; i < _smallest_elements.size() -1; i++)
				if (_smallest_elements[i] == -1)
					_smallest_elements[i] = _smallest_elements[i+1];

		}
	



		void Prepare_Communication() {

			Compute_Smallest_Elementnumber();

			Update_OctreesHelper();

			t_octree_key localsize = _GridIteratorEnd - _OctreeGrid.begin();
			t_octree_key globalsize;
			MPI_Allreduce(&localsize, &globalsize, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

//			if (_My_PID == 0) {
//				COUT << "Grid has a globalsize of: " <<globalsize  << std::endl;
//				COUT << "Min-max: " << _min_max[0] << " - " << _min_max[1] << std::endl << std::flush;
//			}
		
			//Compute boundaries to send
			std::vector< t_octree_iterator> data_send(_Nr_CPU);
			for(int i =0; i < _Nr_CPU; i++) {
				data_send[i] =  lower_bound(_OctreeGrid.begin(), _OctreeGrid.end(), _smallest_elements[i]);
			}
			std::vector<t_index> nr_remote_elem_per_cpu(_Nr_CPU);

			//Copy OctreeKeys to send 
			std::vector<t_octree_key> keys_of_sending_elem(_OctreeGrid.end() - _GridIteratorEnd);
			t_octree_key * sending_addr;

			t_index c =0;
			for(t_octree_iterator iter = _GridIteratorEnd; iter != _OctreeGrid.end(); ++iter) {
				keys_of_sending_elem[c++] = iter->key;
			}
			std::vector<t_index> nr_sending_elements(_Nr_CPU);

			//receiving buffer
			std::vector<t_octree_key> keys_of_remote_elem;
			std::vector<int> keys_of_remote_elem_pos(_Nr_CPU);
			t_octree_key * recv_addr;
			t_index nr_elem;
			t_index nr_elem_recv;

			for(int recv_cpu =0; recv_cpu < _Nr_CPU; recv_cpu++) {
				if (recv_cpu <= _My_PID)
					nr_elem =0;
				else if (recv_cpu != _Nr_CPU-1)
					nr_elem = data_send[recv_cpu +1] -data_send[recv_cpu];
				else
					nr_elem = _OctreeGrid.end() - data_send[recv_cpu];
				
				nr_sending_elements[recv_cpu] = nr_elem;

				MPI_Gather(&nr_elem,1, MPI_INT, &nr_remote_elem_per_cpu[0], 1, MPI_INT, recv_cpu, MPI_COMM_WORLD);


				recv_addr = NULL;
				if (_My_PID == recv_cpu) {
					keys_of_remote_elem_pos[0] = 0;
					for (int i =1; i < _Nr_CPU; i++)
						keys_of_remote_elem_pos[i] = keys_of_remote_elem_pos[i-1] + nr_remote_elem_per_cpu[i-1];

					nr_elem_recv = keys_of_remote_elem_pos[_Nr_CPU-1]+nr_remote_elem_per_cpu[_Nr_CPU-1];
					if (nr_elem_recv > 0 ) {
						keys_of_remote_elem.resize(nr_elem_recv);
						recv_addr = &keys_of_remote_elem[0];
					} else {
						recv_addr = NULL;
					}
				}


				if (0 == nr_elem) {
					sending_addr = NULL;
				} else {
					sending_addr = &keys_of_sending_elem[data_send[recv_cpu] - _GridIteratorEnd];
				}

				MPI_Gatherv(sending_addr, nr_elem, MPI_LONG,
							recv_addr, (int*) &nr_remote_elem_per_cpu[0], (int*) &keys_of_remote_elem_pos[0], MPI_LONG,
							recv_cpu, MPI_COMM_WORLD);

			}
			//add key which should be on this cpu
			std::vector<t_octree_key>::iterator key_iter; 
			std::set<t_octree_key> missingkeys;
			int index;

			for (key_iter = keys_of_remote_elem.begin(); key_iter != keys_of_remote_elem.end(); ++key_iter) {
				index = GetIndexOfKey(*key_iter);
				if (index < 0)
					missingkeys.insert(*key_iter);
			}
			t_index o_size = _OctreeGrid.size();
			_OctreeGrid.reserve(o_size+missingkeys.size());
			std::set< t_octree_key >::iterator iter;
			OctreeNode tmpnode;
			for (iter = missingkeys.begin(); iter != missingkeys.end(); ++iter) {
				tmpnode = OctreeNode(*iter, -3);
				_OctreeGrid.push_back(tmpnode);
			}
			inplace_merge(_OctreeGrid.begin(), _OctreeGrid.begin()+o_size, _OctreeGrid.end());

			Update_OctreesHelper();


			//count how many cpu are sending
			int nr_procs = Count_Nr_CPU(nr_remote_elem_per_cpu);

			_remote_elements_pid.resize(nr_procs);
			_remote_elements_pos.resize(nr_procs+1);
			_remote_elements_pos[0]=0;

			Fill_in_pid_and_positions(nr_remote_elem_per_cpu, _remote_elements_pid, _remote_elements_pos);

			nr_procs = Count_Nr_CPU(nr_sending_elements);
			_sending_elements_pid.resize(nr_procs);
			_sending_elements_pos.resize(nr_procs+1);
			_sending_elements_pos[0]=0;


			Fill_in_pid_and_positions(nr_sending_elements, _sending_elements_pid, _sending_elements_pos);
					
			_remote_elements_index.resize(keys_of_remote_elem.size());
			std::vector<t_index>::iterator remote_iter = _remote_elements_index.begin();

			for (key_iter = keys_of_remote_elem.begin(); key_iter != keys_of_remote_elem.end(); ++key_iter) {
				index = GetIndexOfKey(*key_iter);
				*remote_iter = index;
				remote_iter++;

			}

			_remote_elements_buf.resize(3*_remote_elements_index.size());
			_remote_status = new MPI_Status[_remote_elements_pid.size()];
			_remote_request = new MPI_Request[_remote_elements_pid.size()];
			_sending_status = new MPI_Status[_sending_elements_pid.size()];
			_sending_request = new MPI_Request[_sending_elements_pid.size()];

			return;
		}


		//@}
		int _Nr_CPU;
		int _My_PID;
		t_octree_key _min_max[2];
		std::vector<t_octree_key> _smallest_elements;
		std::vector<t_index> _remote_elements_index;
		std::vector<t_index> _remote_elements_pos;
		std::vector<int> _remote_elements_pid;
		std::vector<double> _remote_elements_buf;

		std::vector<t_index> _sending_elements_pos;
		std::vector<int> _sending_elements_pid;

		MPI_Request *_remote_request;
		MPI_Request *_sending_request;
		MPI_Status *_remote_status;
		MPI_Status *_sending_status;

	private:
		/** 
		 * @brief Copy 3 values from a dynamic Vector to a vector with size 24
		 * 
		 * @param x dynamic vector
		 * @param k index the the block with size 3 in x
		 * @param fix fixed size vector
		 * @param l index of the block with size 3 in fix
		 */
		void Copy3(VectorXd &x, t_index k, Matrix<double, 24, 1> &fix, t_index l)
		{
			for(t_index i=0; i < 3; ++i) 
				fix[3*l+i] = x[3*k+i];
		}

		/** 
		 * @brief Counts the number of CPU that have more then 0 element in
         * the array.
         *
         * @param nr_elem_per_cpu Holds the nr element per CPU
		 */ 
		template <class G>
		int Count_Nr_CPU(const std::vector<G> &nr_elem_per_cpu)
		{
			int nr_procs =0;
			for (unsigned i =0; i < nr_elem_per_cpu.size(); i++) {
				if (nr_elem_per_cpu[i] > 0)
					nr_procs++;
			}
			return nr_procs;
		}

		/** 
		 * @brief Fill in the pid oft the proc and the indexes to send. It is
         * used for the communication.
         *
         * @param[in] nr_elem_per_cpu  number element per CPU
         * @param[out] elements_pid  pid of the CPU
         * @param[out] elements_pos offset of the recv/sending CPU
		 */ 
		void Fill_in_pid_and_positions(const std::vector<t_index> &nr_elem_per_cpu, std::vector<int> &elements_pid, std::vector<t_index> &elements_pos) {
			int k =0;
			for (unsigned i =0; i < nr_elem_per_cpu.size(); i++) {
				if (nr_elem_per_cpu[i] > 0) {
					elements_pid[k] = i;
					elements_pos[k+1] = elements_pos[k]+nr_elem_per_cpu[i];
					k++;
				}
			}
		}

		void Delete_a_vector_Octreenodes(std::vector<OctreeNode> &x) {
			std::vector<OctreeNode> tmp;
			x.swap(tmp);
		}
};

	template <class T>
OctreeGrid<T>::OctreeGrid(ImageReader& ir)
{
	//just some testing for one cpu
	
	int MyPID; 
	MPI_Comm_rank (MPI_COMM_WORLD, &MyPID);
	_My_PID = MyPID;
	MPI_Comm_size(MPI_COMM_WORLD, &_Nr_CPU);

	//Read the data in
	ir.Scan((BaseGrid*) this);
	
	
	GenerateOctree();
	Sort();
	Add_Missing_Nodes();
	Prepare_Communication();


	init_Members();
	if (_grid !=0) {
		delete[] _grid;
		_grid = 0;
	}

}

	template <class T>
void OctreeGrid<T>::init_Members()
{
	for(int i = 0; i< 8; i++)
		_neighbours[i] =0;

	_GridIterator = _OctreeGrid.begin();
	return;
}

// Print out the image in ASCI ART
// only debug reason
	template <class T>
void OctreeGrid<T>::Print()
{
	COUT << "ldim: " << ldim_x << " "<< ldim_y << " "<< ldim_z << "\n";
	COUT << "gdim: " << gdim_x << " "<< gdim_y << " "<< gdim_z << "\n";
	COUT << "Nr. Nodes: " << _nr_nodes << " Nr. Elements"<< _nr_elem << "\n";
	COUT << "Octree length: " << _OctreeGrid.size() << std::endl;
	COUT << "Octree capacity: " << _OctreeGrid.capacity() << "\n";
	COUT << "elemental density: " << ((float) _nr_elem) / ldim_z / ldim_x / ldim_y << "\n";
	COUT << "nodal density: " <<((float) _nr_nodes) / (ldim_z+1) / (ldim_x+1) / (ldim_y+1) << "\n";


	COUT << std::endl << "----------------------------------------------" << std::endl << std::endl;

	for(int z=0; z < ldim_z; z++) {
		for(int y=0; y < ldim_y; y++) {
			for(int x=0; x < ldim_x; x++) {
				COUT <<  _grid[z*ldim_y*ldim_x+y*ldim_x+x] << " ";
			}
			COUT << std::endl;
		}
		COUT << std::endl << "----------------------------------------------" << std::endl << std::endl;
	}
}

template <class T>
void OctreeGrid<T>::GenerateOctree() {
	// some members
	_nr_elem =0;
	_nr_nodes = 0;
	_BC_list.resize(0);
	if (_grid == 0)
	{
		return;
	}


	long ldimxy = ldim_x*ldim_y;
	// copy the emodul in to the picture
	std::pair<std::set<OctreeNode>::iterator, bool> ret;
	std::set<OctreeNode>::iterator it;
	std::set<OctreeNode> octset;
	long tmpx, tmpy, tmpz;
	double w;
	OctreeNode tmp_node(0,0);
	for(long z=0; z < ldim_z; z++) 
		for(long y=0; y < ldim_y; y++) 
			for(long x=0; x < ldim_x; x++) {
				w = _grid[z*ldimxy+y*ldim_x+x]; 
				if  (w > 0) {
					_nr_elem++;
					tmp_node = OctreeNode(KeyGen(x+corner_x,y+corner_y,z+corner_z),w);
					ret = octset.insert(tmp_node);
					for(long tz=0; tz < 2; tz++) 
						for(long ty=0; ty < 2; ty++) 
							for(long tx=0; tx < 2; tx++) {
								if (tz+ty+tx == 0) {
									continue;
								} else {
									// compute temp Coordinate
									tmpx = x +tx;
									tmpy = y +ty;
									tmpz = z +tz;
									it = ret.first;
									tmp_node = OctreeNode(KeyGen(tmpx+corner_x, tmpy+corner_y, tmpz+corner_z), -2.0);
									if  (tmpz == ldim_z ||
										 tmpy == ldim_y ||
										 tmpx == ldim_x ||
										 _grid[tmpz*ldimxy + tmpy*ldim_x+ tmpx] == 0)
									{
										++it;
										tmp_node = OctreeNode(KeyGen(tmpx+corner_x, tmpy+corner_y, tmpz+corner_z), -2.0);
										it = octset.insert(it,tmp_node);
									}
								}
							}
				}
			}
	_OctreeGrid.reserve(octset.size());
	_OctreeGrid.insert(_OctreeGrid.begin(), octset.begin(), octset.end());

	std::set<t_boundary_node> bcset;
	std::vector<unsigned short>::iterator it_coord; 
	std::vector<float>::iterator it_val = fixed_nodes_values.begin();
	for(it_coord = fixed_nodes_coordinates.begin(); it_coord != fixed_nodes_coordinates.end(); ) {
		unsigned short x,y,z,d;
		float disp;
		x = *it_coord; ++it_coord;
		y = *it_coord; ++it_coord;
		z = *it_coord; ++it_coord;
		d = *it_coord; ++it_coord;
		tmp_node = OctreeNode(KeyGen(x,y,z),-1);
		disp = *it_val; ++it_val;
		bcset.insert(t_boundary_node(tmp_node, boundary_disp(d,disp)));
	}
	_BC_list.reserve(bcset.size());
	_BC_list.insert(_BC_list.begin(), bcset.begin(), bcset.end());
	_nr_nodes = _OctreeGrid.size();
	return;
}

template <class T>
std::ostream& OctreeGrid<T>::print(std::ostream& stream) const
{
	stream << "OctreeGrid:\n";
	stream << "   local Dimension: " << ldim[0] << ", " << ldim[1] << ", " << ldim[2]  << "\n";
	stream << "   global Dimension: " <<gdim[0] << ", " << gdim[1] << ", " << gdim[2]  << "\n";
	stream << "   Resolution: " << res[0] << ", " << res[1] << ", " << res[2]  << "\n";
	stream << "   Nr. Nodes: " << _nr_nodes << " Nr. Elements: "<< _nr_elem << "\n";
	stream << "   elemental density: " << ((float) _nr_elem) / ldim_z / ldim_x / ldim_y << "\n";
	stream << "   nodal density: " <<((float) _nr_nodes) / (ldim_z+1) / (ldim_x+1) / (ldim_y+1) << "\n";
	//stream << "   Number BCs: " << _BC_list.size() << " capacity  " << _BC_list.capacity() << std::endl;
	stream << bc;
	return stream;
}


#endif /* OCTREEGRID_H */

