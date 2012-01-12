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


#include "AsciiImage.h"
#include <mpi.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <map>

// A sample Boundary condition. Top and bottom node are fixed.
class BCTopDown {
public:
    BCTopDown(short iz, double disp): _iz(iz), _disp(disp)
    {
    }

    ~BCTopDown() {}

    bool operator()(short x, short y, short z, short d,  double &v) {
        v =0;
        if ( z == 0) {
			if (d == 2) {
				v = _disp;
            	return true;
			}
			return true;

        } else if ( z == _iz ) {
        // only fix in z direction
		//	if (d ==2)
              return true;
		}
        return false;
    }
    short _iz;
    double _disp;

};


AsciiImage::AsciiImage(const char *fi, CPULayout &layout): _file(fi),_layout(layout), MyPID(0)
{
}

AsciiImage::~AsciiImage()
{
}


int AsciiImage::Scan(BaseGrid* grid)
{
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MyPID = mpi_rank;
  PCOUT(MyPID, " Generatin mesh.... \n")

  PCOUT(MyPID, " compute dimension.... \n")

  _procx = _layout.CPUGrid()[0];
  _procy = _layout.CPUGrid()[1];
  _procz = _layout.CPUGrid()[2];

  //setting left corner
  grid->corner_x = 0;
  grid->corner_y = 0;
  grid->corner_z = 0;

  std::map<bcitem, double> bcmap;
  
  //read in only on CPU 0
  if (mpi_rank == 0) {
    std::fstream fb(_file, std::fstream::in);
    fb >> grid->ldim_x;
    fb >> grid->ldim_y;
    fb >> grid->ldim_z;

    for (int i = 0; i<3; i++)
      grid->gdim[i] = grid->ldim[i];

    // copy gdim to ldim
    // read in mat
    int m;
    fb >> m;
    //cout << "m: " << m << endl;
    double *elas = new double[m];
    for( int i=0; i < m; i++) {
      fb >> elas[i];
    }

    //Cube has size 1
    for(int i = 0; i < 3; i++)
      grid->res[i] = 1.0/grid->gdim[1];

    long slice = grid->ldim_x*grid->ldim_y;
    long size  = slice*grid->ldim_z;
    unsigned short r;

    PCOUT(MyPID, " allocate image (double).... \n")
    long ldim_x=grid->ldim_x;
    long ldim_y=grid->ldim_y;
    long ldim_z=grid->ldim_z;
    grid->_grid = new double[size];
    BCTopDown bc(grid->gdim[2],0.04);
    for(long z=0; z < ldim_z; z++)
      for(long y=0; y < ldim_y; y++) {
        for(long x=0; x < ldim_x; x++) {
          fb >> r;
          grid->_grid[z*slice+y*ldim_x+x] = elas[r];
		  if (z == 0 || z == (ldim_z-1)) {
			  if (grid->_grid[z*slice+y*ldim_x+x] != 0)
			  {
				  for (int zi =0; zi < 2; zi++) {
					  for (int yi =0; yi < 2; yi++) {
						  for (int xi =0; xi < 2; xi++) {
							  for (int d = 0; d <3; d++) {
								  double disp =0;
								  if (bc(x+xi,y+yi,z+zi+grid->corner_z, d, disp)) {
									  bcmap[bcitem(x+xi,y+yi,z+zi+grid->corner_z,d)] = disp;
								  }
							  }
						  }
					  }
				  }
			  }
		  }
        }
      }
    fb.close();
    if (mpi_size > 1) {
      MPI_Bcast(grid->ldim, 3,MPI_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(grid->gdim, 3,MPI_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(grid->res, 3,MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    delete[] elas;
  }
  else
  {
    MPI_Bcast(grid->ldim, 3,MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(grid->gdim, 3,MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(grid->res, 3,MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
    grid->_grid =0;
  }
  for(std::map<bcitem, double>::iterator it =  bcmap.begin(); it != bcmap.end(); ++it) {
	  grid->fixed_nodes_coordinates.push_back(it->first.x);
	  grid->fixed_nodes_coordinates.push_back(it->first.y);
	  grid->fixed_nodes_coordinates.push_back(it->first.z);
	  grid->fixed_nodes_coordinates.push_back(it->first.d);
	  grid->fixed_nodes_values.push_back(it->second);
	  
  }



    PCOUT(MyPID, "AsciiReader: \n")
    PCOUT(MyPID, "  local Dimension: " << grid->ldim[0] << " " << grid->ldim[1] << " " << grid->ldim[2] << " Resolution: " << grid->res[0] << "\n")
    PCOUT( MyPID, "  BC Size: " << grid->fixed_nodes_values.size() << "\n")

      return 0;
  }
