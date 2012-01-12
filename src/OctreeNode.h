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

#ifndef OCTREENODE_H
#define OCTREENODE_H


#include "Config.h"
/*! This class implements the item of a grid.
 *  Such a item has to implement following operators: (), <, == and
 *  must have a key and a w member.
 */
class OctreeNode
{
public:
	OctreeNode(): key(0), w(0)
			{
			}
	OctreeNode(t_octree_key key, double w): key(key), w(w)
	{
	}

	operator t_octree_key() const
	{
		return key;
	}

	bool operator <(const OctreeNode &b) const
	{
		return key < b.key;
	}	
	bool operator ==(const OctreeNode &b) const
	{
		return key == b.key;
	}	
		
	t_octree_key key;
	double w;
	
private:
};

#endif
