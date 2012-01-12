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



#ifndef KEYGENERATOR_H
#define KEYGENERATOR_H

#include "Config.h"
#include <bitset>
#include <iostream>


/*! The OctreeKey_Loop implements the mapping of a node from the coordinates to
 *  its number by a functor with a loop
 */
class OctreeKey_Loop {
	public:
		t_octree_key operator()(t_coord x, t_coord y, t_coord z)
		{
			t_coord mask = 1;
			t_octree_key res = 0;
			for(int i = 0; i < 16; i++) {
				res += ((t_octree_key)x & mask) << (2*i);  // 3*i - i;
				res += ((t_octree_key)y & mask) << (2*i+1);
				res += ((t_octree_key)z & mask) << (2*i+2);
				mask = mask << 1;
			}
			return res;
		}
};

/*! This class implements the mapping from the coordinates to node value with
 *  fast loopkup. It also implements to compute the new node value from the old
 *  if a coordinate in one direction is increased by one.
 */
class OctreeKey_Lookup {
	public:
		OctreeKey_Lookup()
		{
			unsigned int tmp =0, mask =1;
			for(unsigned short x =0; x < 256; x++) {
				tmp =0;
				mask = 1;
				for(int i = 0; i < 8; i++) {
					tmp += ((unsigned int)x & mask) << 2*i;
					mask = mask << 1;
				}
				table[x] = tmp;
			}
			t_coord y = -1;
			t_coord z = -1;
			maskx = operator()(0,y,z);
			masky = 1+ (maskx<< 1); 
			maskz = 3+ (maskx<< 2); 
		}

		t_octree_key operator()(t_coord x, t_coord y, t_coord z)
		{
			t_coord mask = 255;
			t_octree_key res = 0;
			for(int i = 0; i < 2; i++) {
				res += ((t_octree_key) table[x & mask]) << (8*3*i);
				res += ((t_octree_key) table[y & mask]) << (8*3*i+1);
				res += ((t_octree_key) table[z & mask]) << (8*3*i+2);
				x = x >> 8;
				y = y >> 8;
				z = z >> 8;
			}
			return res;
		}

		t_octree_key incX(t_octree_key key)
		{
			t_octree_key tmp = (key | maskx)+1;
			return (tmp&~maskx) | (key&maskx);
		}

		t_octree_key incY(t_octree_key key)
		{
			t_octree_key tmp = (key | masky)+1;
			return (tmp&~masky) | (key&masky);
		}

		t_octree_key incZ(t_octree_key key)
		{
			t_octree_key tmp = (key | maskz)+1;
			return (tmp&~maskz) | (key&maskz);
		}
		unsigned int table[256]; 
		t_octree_key maskx;
		t_octree_key masky;
		t_octree_key maskz;
};

#endif /*KEYGENERATOR_H*/
