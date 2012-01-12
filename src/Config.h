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


#ifndef CONFIG_H
#define CONFIG_H

#include <utility>

#define COUT std::cout

/* 32bit is enough for the index */
/// type of the local array index
typedef unsigned int t_index;

/// type of the coordinate
typedef unsigned int t_coord;

/// type of the octree key
typedef long t_octree_key;


/// holds direction and displacement of the BC
struct boundary_disp
{
	boundary_disp(short dir=0, float disp=0):dir(dir), disp(disp)
	{
	}

	short dir;
	float disp;
	bool operator < (const boundary_disp &rhs) const {
		return dir < rhs.dir;
	}
};

/// A boundary node is consist of an node and bc properties
typedef std::pair<t_octree_key, boundary_disp> t_boundary_node;



//some Marcros

#define COUT_RANK0 std::cout
//cout only on cpunode 0
#define PCOUT(PID, msg) \
	    if (PID == 0) COUT_RANK0 << msg;

//cout on all cpus with barrier 
#define ALL_OUT(PID, SIZE, msg) \
{ \
    int i=0; \
    MPI_Status status; \
    MPI_Barrier(MPI_COMM_WORLD); \
    if ( 0 != PID) { \
        MPI_Recv(&i, 1, MPI_INT, PID -1,0, MPI_COMM_WORLD, &status); \
    } \
 \
    COUT << msg; \
    COUT.flush(); \
    if (PID != SIZE -1) { \
        MPI_Send(&i, 1, MPI_INT, PID+1, 0, MPI_COMM_WORLD);  \
    } \
    MPI_Barrier(MPI_COMM_WORLD); \
}

#endif /*CONFIG_H*/
