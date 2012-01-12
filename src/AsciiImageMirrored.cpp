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



#include "AsciiImageMirrored.h"
#include <mpi.h>
#include <cstring>
#include <iostream>
#include <fstream>

AsciiImageMirrored::AsciiImageMirrored(const char *fi, CPULayout &layout): _file(fi),_layout(layout), MyPID(0)
{
}

AsciiImageMirrored::~AsciiImageMirrored()
{
}

int AsciiImageMirrored::Scan(BaseGrid* grid)
{
	int mpi_rank, mpi_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MyPID = mpi_rank;
	PCOUT(MyPID, " Generatin mesh.... \n")


	PCOUT(MyPID, " compute dimension.... \n")

	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	_procx = _layout.CPUGrid()[0];
	_procy = _layout.CPUGrid()[1];
	_procz = _layout.CPUGrid()[2];

	int pdimx = _layout.CPUCoord()[0];
	int pdimy = _layout.CPUCoord()[1];
	int pdimz = _layout.CPUCoord()[2];

	grid->corner_x = 0; 
	grid->corner_y = 0;
	grid->corner_z = 0;

	// only cpu one
	if (mpi_rank == 0) {
		std::fstream fb(_file, std::fstream::in);
		fb >> grid->ldim_x;
		fb >> grid->ldim_y;
		fb >> grid->ldim_z;

		for (int i = 0; i<3; i++)
			grid->gdim[i] = grid->ldim[i]*_proc[i];

		//for (int i = 0; i<3; i++)
		//cout << grid->ldim[i] << " ";

		// copy gdim to ldim
		// read in mat
		int m;
		fb >> m;
		//cout << "m: " << m << endl;
		double *elas = new double[m];
		for( int i=0; i < m; i++) {
			fb >> elas[i];
			//if (elas[i] == 0)
			//  elas[i] = 0.001;
		}

		//cout << "\nelas: ";
		//for( int i=0; i < m; i++) 
		//  cout <<  elas[i] << " ";
		//cout << endl;


		// Setting sizes 
		// 22 on 5 proc will be 4 4 5 4 5
		// for (int i = 0; i < 3; i++)
		//   grid->ldim[i] = grid->gdim[i];
		//    grid->ldim[i] = long((float)grid->gdim[i]/(float)grid->proc[i] * (float)(pdimx+1))- long((float)grid->gdim[i]/(float)grid->proc[i] *(float)(pdimx));

		// Set cornercord
		//  for(int i = 0; i < 3; i++)
		//   grid->cornercord[i] = long((float)_dim/(float)grid->proc[i] *(float)(grid->cpucord[i])); 

		// if the there is en upper neighbour, import ghostlayer
		//  if (grid->cpu_neighbour[1] > -1) grid->ldim_x++; //ghost layer
		//  if (grid->cpu_neighbour[3] > -1) grid->ldim_y++; //ghost layer
		//  if (grid->cpu_neighbour[5] > -1) grid->ldim_z++; //ghost layer

		//Cube has size 1
		for(int i = 0; i < 3; i++)
			grid->res[i] = 1.0/grid->gdim[i];

		long slice = grid->ldim_x*grid->ldim_y;
		long size  = slice*grid->ldim_z;
		unsigned short r;

		PCOUT(MyPID, " allocate image (double).... \n")
			long ldim_x=grid->ldim_x;
		long ldim_y=grid->ldim_y;
		long ldim_z=grid->ldim_z;
		grid->_grid = new double[size];
		for(long z=0; z < ldim_z; z++)
			for(long y=0; y < ldim_y; y++) {
				for(long x=0; x < ldim_x; x++) {
					fb >> r;
					grid->_grid[z*slice+y*ldim_x+x] = elas[r];
				}
			}
		fb.close();
		if (mpi_size > 1) {
			MPI_Bcast(grid->ldim, 3,MPI_LONG, 0, MPI_COMM_WORLD);
			MPI_Bcast(grid->gdim, 3,MPI_LONG, 0, MPI_COMM_WORLD);
			MPI_Bcast(grid->res, 3,MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(grid->_grid, size,MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}
		delete[] elas;
	}
	else
	{
		MPI_Bcast(grid->ldim, 3,MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(grid->gdim, 3,MPI_LONG, 0, MPI_COMM_WORLD);
		MPI_Bcast(grid->res, 3,MPI_DOUBLE, 0, MPI_COMM_WORLD);

		grid->corner[0] = grid->ldim[0]*_layout.CPUCoord()[0]; 
		grid->corner[1] = grid->ldim[1]*_layout.CPUCoord()[1]; 
		grid->corner[2] = grid->ldim[2]*_layout.CPUCoord()[2]; 

		long ldim_x=grid->ldim_x;
		long ldim_y=grid->ldim_y;
		long ldim_z=grid->ldim_z;
		long slice = grid->ldim_x*grid->ldim_y;
		long size  = slice*grid->ldim_z;
		grid->_grid = new double[size];
		double *tmp = new double[size];

		MPI_Bcast(tmp, size ,MPI_DOUBLE, 0, MPI_COMM_WORLD);
		long nz, ny, nx;
		for(long z=0; z < ldim_z; z++) {
			nz = z;
			if (pdimz % 2 == 1)
				nz = ldim_z-1 - z;
			for(long y=0; y < ldim_y; y++) {
				ny = y;
				if (pdimy % 2 == 1)
					ny = ldim_y-1 - y;
				for(long x=0; x < ldim_x; x++) {
					nx = x;
					if (pdimx % 2 == 1)
						nx = ldim_x-1 - x;
					grid->_grid[nz*slice+ny*ldim_x+nx] = tmp[z*slice+y*ldim_x+x];
				}
			}
		}
		delete[] tmp;
	}

	PCOUT(MyPID, "AsciiReader: \n")
		PCOUT(MyPID, "  local Dimension: " << grid->ldim[0] << " " << grid->ldim[1] << " " << grid->ldim[2] << " Resolution: " << grid->res[0] << "\n")

		return 0;
}

