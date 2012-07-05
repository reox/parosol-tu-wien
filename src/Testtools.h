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


#ifndef TESTTOOLS_H
#define TESTTOOLS_H

#include "Config.h"
#include "OctreeGrid.h"


/*! PrintVector is a helper class to display vectors with different orderings of the items.
 *
 * It is used for debung reason. It is not testet in parallel. 
 */
template <class T>
class PrintVector {
	public:
		/** 
		 * @brief Constructor
		 * 
		 * @param grid the finite element grid
		 */
		PrintVector(OctreeGrid<T> &grid): _grid(grid)
		{
			t_coord dims[3];
			grid.GetGlobalDim(dims);
			x = dims[0] + 1;
			y = dims[1] + 1;
			z = dims[2] + 1;
		}
		
	  	/** 
		 * @brief Converts coordinates into a morton key.
		 * 
		 * @param x x coordinate
		 * @param y y coordinate
		 * @param z z coordinate
		 * @return morton key
		 */
		long nTOo(int x, int y, int z)
		{
			int mask = 1;
			long res = 0;
			for(int i = 0; i < 16; i++) {
				res += ((long)x & mask) << (2*i);  // 3*i - i;
				res += ((long)y & mask) << (2*i+1);
				res += ((long)z & mask) << (2*i+2);
				mask = mask << 1;
			}
			return res;
		}

	  	/** 
		 * @brief Converts a morton key into coordinates.
		 * 
		 * @param [out] x x coordinate
		 * @param [out] y y coordinate
		 * @param [out] z z coordinate
		 * @param [in] key morton key
		 */
		void oTOn(long key, int &x, int &y, int &z)
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

	  	/** 
		 * @brief Changes a vector from lexicographical ordering to morton ordering.
		 * 
		 * @param [out] oct Vector in morton ordering
		 * @param [in] block number of points per node
		 * @param [in] normal Vector in lexicographical order
		 */
		void NormalToOctgrid(const Eigen::VectorXd &normal, Eigen::VectorXd &oct, t_index block)
		{
			t_index ind;
			for(t_coord zi =0; zi < z; zi++) {
				for(t_coord yi =0; yi < y; yi++) {
					for(t_coord xi =0; xi < x; xi++) {
						ind = _grid.GetIndexOfKey(nTOo(xi, yi,zi)); 
						if (ind >=0) 
							for(t_index i =0; i < block; i++) {
								oct[ind*block+i] = normal[((zi*y+yi)*x+xi)*block+i];
						}
					}
				}
			}
		}

	  	/** 
		 * @brief Changes a vector from morton ordering to lexicographical ordering.
		 * 
		 * @param [in] oct Vector in morton ordering
		 * @param [in] block number of points per node
		 * @param [out] normal Vector in lexicographical order
		 */
		void OctgridToNormal(const Eigen::VectorXd &oct, Eigen::VectorXd &normal, int block)
		{
			normal.setZero(x*y*z*block);
			t_octree_key k;
			int xi, yi, zi;
			int ind =0;
			int normal_index;
			t_octree_iterator iter;
			/*cout << "block " << block << endl;
			cout << "x: " << x << " ";
			cout << "y: " << y << " ";
			cout << "z: " << z << " ";
			cout << endl; */
			std::vector<OctreeNode> &grid = _grid.GetOctGrid();
			for(iter = grid.begin(); iter != grid.end(); ++iter ) {
				k = iter->key;
				oTOn(k, xi, yi, zi);
			/*	cout << "xi: " << xi << " ";
				cout << "yi: " << yi << " ";
				cout << "zi: " << zi << " ";
				cout << "key: " << k << " ";
				cout << endl;*/
				normal_index = (zi*y+yi)*x+xi;
			//	cout << "normal_index: " << normal_index << " ind: " << ind;
			//	cout << " block " << block << endl;

				for(int i =0; i < block; i++) {
					normal[normal_index*block+i] = oct[ind*block +i];
				}
				ind++;
			}
		}

	  	/** 
		 * @brief Displays a morton ordered vector in a lexicographical ordering.
		 * 
		 * @param [in] v Vector in morton ordering
		 * @param [in] block number of points per node
		 */
		void PrintVectorOct(const Eigen::VectorXd v, t_coord block) {
            Eigen::VectorXd norm(x*y*z*block);
			OctgridToNormal(v, norm, block);
			PrintVectorNormal(norm, block);
			return;
		}

	  	/** 
		 * @brief Displays a vector without changing the ordering.
		 * 
		 * @param [in] v Vector in morton ordering
		 * @param [in] block number of points per node
		 */
		void PrintVectorNormal(const Eigen::VectorXd v, t_coord block)
		{
			int i;
			for(t_coord zi =0; zi < z; zi++) {
				for(t_coord yi =0; yi < y; yi++) {
					for(t_coord xi =0; xi < x; xi++) {
						i = (zi*y+yi)*x+xi;
						for(t_coord j = 0; j < block; j++) {
							std::cout << v[3*i+j] << " ";
						}
						std::cout << std::endl;
					}
				}
				std::cout << std::endl;
			}
			return;
		}


		OctreeGrid<T> &_grid;

		//! Nr nodes in each dimension
		t_coord x,y,z;
};
#endif /* TESTTOOLS_H */
