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

#ifndef MLOCTREEGRID_H
#define MLOCTREEGRID_H

#include <ostream>

#include <vector>

#include <utility>
#include <algorithm>

#include "OctreeGrid.h"


/*! The MlOctreeGrid is base on the standard grid. It adds the restriction and the prolongation operations.
 * The corresponding key of the coarse grid is the fine grid key divided by 8.
 */
template <typename T>
class MlOctreeGrid : public OctreeGrid<T>
{
	public:
		MlOctreeGrid(OctreeGrid<T> *finegrid):_finegrid(*finegrid) 
		{
			OctreeGrid<T>::_My_PID = finegrid->_My_PID;
			OctreeGrid<T>::_Nr_CPU = finegrid->_Nr_CPU;
			Coarse();
		}

		virtual ~MlOctreeGrid()
		{
		}

		void Restrict(VectorXd &fine, VectorXd &coarse);
		void Prolongate(VectorXd &coarse, VectorXd &fine);

	private:
		void Coarse();

		OctreeGrid<T> &_finegrid;

		void coarse_member() {
			t_coord dims[3];
			_finegrid.GetGlobalDim(dims);
			for(int i =0; i < 3; i++)
				OctreeGrid<T>::gdim[i] = (dims[i]+1)/2; 

			_finegrid.GetLocalDim(dims);
			for(int i =0; i < 3; i++)
				OctreeGrid<T>::ldim[i] = (dims[i]+1)/2; 

			double res[3];
			_finegrid.GetRes(res);
			for(int i =0; i < 3; i++)
				OctreeGrid<T>::res[i] = res[i]*2;

			OctreeGrid<T>::_smallest_elements.resize(OctreeGrid<T>::_Nr_CPU+1);
			for(int i =0; i < OctreeGrid<T>::_Nr_CPU+1; i++)
				OctreeGrid<T>::_smallest_elements[i] = _finegrid._smallest_elements[i]/8;
			OctreeGrid<T>::_smallest_elements[OctreeGrid<T>::_Nr_CPU] +=8; //round up

		}

		int SetBoundaryNodes(BoundaryCondition &bc , VectorXd &x)
		{
			t_index num_ind = bc.FixedNodes_Ind.size();
			// write the disp in the vector 
			for (t_index i=0; i<num_ind; ++i) {
				x[bc.FixedNodes_Ind[i]] = 0; 
			}
			return 0;
		}
	
	void Distribute_Elements(std::set<OctreeNode> &octset) {
		std::set<OctreeNode>::iterator it;
		std::pair<std::set<OctreeNode>::iterator, bool> ret;
		OctreeNode tmp_node;

		//copy things in to an vector
		it = octset.lower_bound(OctreeNode(OctreeGrid<T>::_smallest_elements[OctreeGrid<T>::_My_PID+1],0));
		OctreeGrid<T>::_OctreeGrid.resize(0);
		OctreeGrid<T>::_OctreeGrid.reserve(octset.size());
		OctreeGrid<T>::_OctreeGrid.insert(OctreeGrid<T>::_OctreeGrid.begin(), it, octset.end());

		std::vector<OctreeNode> offcpu_elem;

		std::vector<OctreeNode>::iterator Oct_iter = OctreeGrid<T>::_OctreeGrid.begin();
		std::vector< std::vector<OctreeNode>::iterator> send_pos_iters(OctreeGrid<T>::_Nr_CPU+1);
		std::vector<int> send_count(OctreeGrid<T>::_Nr_CPU);
		std::vector<int> send_pos(OctreeGrid<T>::_Nr_CPU);
		std::vector<int> recv_count(OctreeGrid<T>::_Nr_CPU);
		std::vector<int> recv_pos(OctreeGrid<T>::_Nr_CPU);

		send_pos_iters[0] = OctreeGrid<T>::_OctreeGrid.begin();
		send_pos_iters[OctreeGrid<T>::_Nr_CPU] = OctreeGrid<T>::_OctreeGrid.end();

		//Prepare Position which belongs to other CPU
		for(int i = 0; i < OctreeGrid<T>::_Nr_CPU; i++)
		{
			while( (Oct_iter != OctreeGrid<T>::_OctreeGrid.end()) && (Oct_iter->key < OctreeGrid<T>::_smallest_elements[i+1]))
				Oct_iter++;
			send_pos_iters[i+1] = Oct_iter;
		}

		//Compute how much date to send
		for(int i = 0; i < OctreeGrid<T>::_Nr_CPU; i++) {
			send_count[i] = send_pos_iters[i+1] - send_pos_iters[i];
		}
		for(int i = 0; i < OctreeGrid<T>::_Nr_CPU; i++) {
			send_count[i] *= sizeof(OctreeNode);
		}

		MPI_Alltoall (&send_count[0], 1, MPI_INT, &recv_count[0], 1, MPI_INT, MPI_COMM_WORLD);

		//Compute the Positions
		send_pos[0] = 0;
		recv_pos[0] = 0;
		for(int i=1; i < OctreeGrid<T>::_Nr_CPU; i++) {
			send_pos[i] = (send_pos[i-1] + send_count[i-1]);
			recv_pos[i] = recv_pos[i-1] + recv_count[i-1];
		}

		int newsize = (recv_pos[OctreeGrid<T>::_Nr_CPU-1] + recv_count[OctreeGrid<T>::_Nr_CPU-1])/sizeof(OctreeNode);
		offcpu_elem.resize(newsize);

		void * send_addr;
		void * recv_addr;
		if (OctreeGrid<T>::_OctreeGrid.size() > 0)
			send_addr = &OctreeGrid<T>::_OctreeGrid[0];
		else
			send_addr = NULL;

		if (offcpu_elem.size() > 0)
			recv_addr = &offcpu_elem[0];
		else
			recv_addr = NULL;

		MPI_Alltoallv(send_addr, &send_count[0], &send_pos[0], MPI_CHAR,
				recv_addr, &recv_count[0], &recv_pos[0], MPI_CHAR,
				MPI_COMM_WORLD);

		//include int to the set

		for(t_index i =0; i< offcpu_elem.size(); i++) {
			tmp_node = offcpu_elem[i];
			ret = octset.insert(tmp_node);
			if (ret.second) {   //inserted
				if(tmp_node.w > 0) {
					OctreeGrid<T>::_nr_elem++;
				}
			} else {     // already in the set update the weight when it is an element
				if (tmp_node.w >0) {
					if (ret.first->w >0) {   // add up only if an element is in the set
						tmp_node.w += ret.first->w;
					} 
					octset.erase(ret.first);
					ret = octset.insert(tmp_node);
				}
			}
		}
		OctreeGrid<T>::_OctreeGrid.resize(0);
		return;
	}


};
template <class T>
void MlOctreeGrid<T>::Coarse()
{
	coarse_member();
	t_octree_iterator iter;

	std::pair<std::set<OctreeNode>::iterator, bool> ret;
	std::set<OctreeNode>::iterator it;
	std::set<OctreeNode>::iterator it_inserted;
	std::set<OctreeNode> octset;
	OctreeNode tmp_node(0,0);
	t_octree_key coarsekey;
	t_octree_key nodekeys[8];
	double weight;

	OctreeGrid<T>::_nr_elem =0;
	//coarse the grid
	for(iter = _finegrid._OctreeGrid.begin(); iter != _finegrid._GridIteratorEnd; ++iter) {
		if (iter->w > 0) {
			weight = iter->w/8;
			// cut the old coordinate away
			coarsekey = iter->key >> 3; 
			tmp_node = OctreeNode(coarsekey, weight);
			ret = octset.insert(tmp_node);
			if (ret.second) {   //inserted
				OctreeGrid<T>::_nr_elem++;
			} else {     // already in the set update the weight
				tmp_node.w += ret.first->w;
				octset.erase(ret.first);
				ret = octset.insert(tmp_node);
			}
		}
	}

	// Add missing nodes
	for(it = octset.begin(); it != octset.end(); ++it) {
			if (it->w >0) {
				OctreeGrid<T>::ComputeKeysofNodes(nodekeys, it->key);
				it_inserted = it;
				for(int i = 1; i < 8; i++) {
					tmp_node =  OctreeNode(nodekeys[i], -2.0);
					octset.insert(tmp_node);
				}
			}

	}

	t_octree_key min = (((t_octree_key) 1)<< 50);
	t_octree_key max = -1;
	if (octset.size() > 0) {
		min = octset.begin()->key;
		max = octset.rbegin()->key;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	OctreeGrid<T>::Compute_min_max(min, max);

	Distribute_Elements(octset);

// adding missing nodes from prolongation

	for(iter = _finegrid._OctreeGrid.begin(); iter < _finegrid._GridIteratorEnd; ++iter) {
		unsigned int LocalFineCoordinate = (iter->key) & 7;
		t_octree_key newkey;

		short max_x = (LocalFineCoordinate & 1);
		short max_y = (LocalFineCoordinate & 2);
		short max_z = (LocalFineCoordinate & 4);
		for(short z =0; z <= max_z; z+=4) {
			for(short y =0; y <= max_y; y+=2) {
				for(short x =0; x <= max_x; x++) {
					newkey = (iter->key)/8 ;
					if (z!=0)
						newkey = OctreeGrid<T>::KeyGen.incZ(newkey);
					if (y!=0)
						newkey = OctreeGrid<T>::KeyGen.incY(newkey);
					if (x!=0)
						newkey = OctreeGrid<T>::KeyGen.incX(newkey);

					tmp_node = OctreeNode(newkey, -1.1111);

					octset.insert(tmp_node);
				}
			}
		}
	}

	OctreeGrid<T>::_OctreeGrid.reserve(octset.size());
	OctreeGrid<T>::_OctreeGrid.resize(0);
	OctreeGrid<T>::_OctreeGrid.insert(OctreeGrid<T>::_OctreeGrid.begin(), octset.begin(), octset.end());
    OctreeGrid<T>::_nr_nodes = OctreeGrid<T>::_OctreeGrid.size();

	OctreeGrid<T>::Prepare_Communication();

	std::vector<t_index>::iterator bc_it;
	t_octree_key finekey;
	typedef t_boundary_node bcpair;
	std::set<bcpair> bcset;
	tmp_node.w = -1; // Does it work

	short LocalFineCoordinate;
	short max_x, max_y, max_z;
	t_octree_key new_bckey;
	boundary_disp bc_value;
	short direction;

	for(bc_it = _finegrid.bc.FixedNodes_Ind.begin(); bc_it != _finegrid.bc.FixedNodes_Ind.end(); ++bc_it){ // , ++i) {
		finekey  = _finegrid._OctreeGrid[*bc_it/3];
		direction = *bc_it%3;
		
		tmp_node.key = finekey >> 3;
		LocalFineCoordinate = finekey & 7;

		max_x = (LocalFineCoordinate & 1);
		max_y = (LocalFineCoordinate & 2);
		max_z = (LocalFineCoordinate & 4);

		for(short z =0; z <= max_z; z+=4) {
			for(short y =0; y <= max_y; y+=2) {
				for(short x =0; x <= max_x; x++) {
					new_bckey = finekey >> 3;
					if (z!=0)
						new_bckey = OctreeGrid<T>::KeyGen.incZ(new_bckey);
					if (y!=0)
						new_bckey = OctreeGrid<T>::KeyGen.incY(new_bckey);
					if (x!=0)
						new_bckey = OctreeGrid<T>::KeyGen.incX(new_bckey);

					//TODO only for z
					tmp_node.key = new_bckey;
				   	
					bcset.insert(bcpair(tmp_node, boundary_disp(direction, 0)));
				}
			}
		}
	}

	OctreeGrid<T>::_BC_list.reserve(bcset.size());
	OctreeGrid<T>::_BC_list.insert(OctreeGrid<T>::_BC_list.begin(), bcset.begin(), bcset.end());

	OctreeGrid<T>::init_Members();
}

	template <class T>
void MlOctreeGrid<T>::Restrict(VectorXd &fine, VectorXd &coarse)
{
	_finegrid.Recv_import_Ghost(fine);
	_finegrid.Send_import_Ghost(fine);
	t_octree_vec& finegrid = _finegrid.GetOctGrid();
	t_octree_iterator fine_iter;
	t_index fineindex =0;

	t_octree_key coarsekey = -1; // biggest point. the first point is always smaller
	t_octree_key new_coarsekey;
	t_octree_iterator coarse_iter = OctreeGrid<T>::_OctreeGrid.begin();
	t_octree_iterator new_coarse_iter = OctreeGrid<T>::_OctreeGrid.begin();
	int new_coarseindex;
	

	unsigned short LocalFineCoordinate;
	int max_x, max_y, max_z;
	double factor;
	_finegrid.Wait_import_Ghost();
	OctreeGrid<T>::Recv_export_Ghost();

	for(fine_iter = finegrid.begin(); fine_iter != _finegrid._GridIteratorEnd; ++fine_iter, fineindex +=3) {
		coarsekey = new_coarsekey = (fine_iter->key >> 3);
		LocalFineCoordinate = fine_iter->key & 7;

		t_octree_key diffx = OctreeGrid<T>::KeyGen.incX(new_coarsekey) ^ new_coarsekey;
		t_octree_key diffy = OctreeGrid<T>::KeyGen.incY(new_coarsekey) ^ new_coarsekey;
		t_octree_key diffz = OctreeGrid<T>::KeyGen.incZ(new_coarsekey) ^ new_coarsekey;

		max_x = (LocalFineCoordinate & 1);
		max_y = (LocalFineCoordinate & 2);
		max_z = (LocalFineCoordinate & 4);

		factor = ((max_z >> 2) + 1) * ((max_y >> 1) + 1) * (max_x  + 1);
		factor = 1.0/ (factor);//*8); 

		for(int z =0; z <= max_z; z+=4) {
			for(int y =0; y <= max_y; y+=2) {
				for(int x =0; x <= max_x; x++) {
					new_coarsekey = coarsekey;
					if (z!=0)
						new_coarsekey = new_coarsekey ^ diffz;
					if (y!=0)
						new_coarsekey = new_coarsekey ^ diffy;
					if (x!=0)
						new_coarsekey = new_coarsekey ^ diffx;
					new_coarse_iter = coarse_iter;
					new_coarseindex = OctreeGrid<T>::GetIndexOfKeyA(new_coarse_iter, new_coarsekey);
					for(int i = 0; i < 3; i++) {
						coarse[new_coarseindex*3+i] += fine[fineindex +i]*factor;
					}
				}
			}
		}
	}
	OctreeGrid<T>::Send_export_Ghost(coarse);
	OctreeGrid<T>::WaitAndCopy_export_Ghost(coarse);

	SetBoundaryNodes(OctreeGrid<T>::bc, coarse);
	return;
}

template <class T>
void MlOctreeGrid<T>::Prolongate(VectorXd &coarse, VectorXd &fine)
{
	OctreeGrid<T>::Recv_import_Ghost(coarse);
	SetBoundaryNodes(OctreeGrid<T>::bc, coarse);
	t_octree_vec& finegrid = _finegrid.GetOctGrid();
	t_octree_iterator fine_iter;
	t_index fineindex =0;

	t_octree_key coarsekey = -1; // biggest point. the first point is always smaller
	t_octree_key new_coarsekey;
	t_octree_iterator coarse_iter = OctreeGrid<T>::_OctreeGrid.begin();
	t_octree_iterator new_coarse_iter = OctreeGrid<T>::_OctreeGrid.begin();
	t_index new_coarseindex;

	unsigned short LocalFineCoordinate;
	int max_x, max_y, max_z;
	double factor;


	OctreeGrid<T>::Send_import_Ghost(coarse);
	OctreeGrid<T>::Wait_import_Ghost();

	for(fine_iter = finegrid.begin(); fine_iter !=_finegrid._GridIteratorEnd; ++fine_iter, fineindex +=3) {
		coarsekey = new_coarsekey = (fine_iter->key >> 3);
		LocalFineCoordinate = fine_iter->key & 7;

		t_octree_key diffx = OctreeGrid<T>::KeyGen.incX(new_coarsekey) ^ new_coarsekey;
		t_octree_key diffy = OctreeGrid<T>::KeyGen.incY(new_coarsekey) ^ new_coarsekey;
		t_octree_key diffz = OctreeGrid<T>::KeyGen.incZ(new_coarsekey) ^ new_coarsekey;

		max_x = (LocalFineCoordinate & 1);
		max_y = (LocalFineCoordinate & 2);
		max_z = (LocalFineCoordinate & 4);

		factor = ((max_z >> 2) + 1) * ((max_y >> 1) + 1) * (max_x  + 1);
		factor = 1.0/ factor; 


		for(int z =0; z <= max_z; z+=4) {
			for(int y =0; y <= max_y; y+=2) {
				for(int x =0; x <= max_x; x++) {
					new_coarsekey = coarsekey;
					if (z!=0)
						new_coarsekey = new_coarsekey ^ diffz;
					if (y!=0)
						new_coarsekey = new_coarsekey ^ diffy;
					if (x!=0)
						new_coarsekey = new_coarsekey ^ diffx;
					new_coarse_iter = coarse_iter;
					new_coarseindex = OctreeGrid<T>::GetIndexOfKeyA(new_coarse_iter, new_coarsekey);
					for(int i = 0; i < 3; i++) {
						fine[fineindex +i] += coarse[new_coarseindex*3+i]*factor;
					}
				}
			}
		}
	}
//	_finegrid.Send_export_Ghost(fine);
//	_finegrid.WaitAndCopy_export_Ghost(fine);

	SetBoundaryNodes(_finegrid.bc, fine);
	return;
}
#endif 

