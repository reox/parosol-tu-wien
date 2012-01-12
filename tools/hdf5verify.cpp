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

//mpic++ -Wall -lhdf5 -o asfasf hdf5verify.cpp ../src/GReader.cpp

//check the BC of an HDF5 file
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <string>


#include "../src/GReader.hpp"


int main (int argc, char * argv[])
{
	if (argc < 2) {
		std::cout << "usage " << argv[0] << " hdf5file" << std::endl;
		return 1;
	}
	char * hdf5filename = argv[1];
	std::string file(hdf5filename);
	HDF5_GReader reader(file);
	hsize_t global_dims_of_hdf5[3];
	hsize_t my_offset[3] = {};
	reader.Select("Image_Data");
	reader.GetSizeOfDataset("Image",global_dims_of_hdf5, 3);
	long imagesize = global_dims_of_hdf5[0]*global_dims_of_hdf5[1]*global_dims_of_hdf5[2] ;
	float *image = new float[imagesize];
	reader.Read("Image", image, my_offset, global_dims_of_hdf5, 3);

	hsize_t bc_dims[2] = {};
	reader.GetSizeOfDataset("Fixed_Displacement_Coordinates", bc_dims, 2);
	std::vector<short> bccoords(bc_dims[0]*bc_dims[1]);
	reader.Read("Fixed_Displacement_Coordinates",&bccoords[0], my_offset, bc_dims, 2);
	
	short x,y,z;
	short ix, iy, iz;
	bool flag1,flag2;
	int count =0;
	long index;
	for(unsigned int i=0; i < bccoords.size();) {
		z = bccoords[i++];
		y = bccoords[i++];
		x = bccoords[i++];
//		std::cout << x << y << z << std::endl;
		i++;
		flag1 = false;
		for(short a =-1; a <1;a++) {
			for(short b =-1; b <1;b++) {
				for(short c =-1; c <1;c++) {
					ix = x+a;
					iy = y+b;
					iz = z+c;
					if (iz >= global_dims_of_hdf5[0])
						continue;
					if (iy >= global_dims_of_hdf5[1])
						continue;
					if (ix >= global_dims_of_hdf5[2])
						continue;
					if (iz < 0)
						continue;
					if (iy < 0 )
						continue;
					if (ix < 0)
						continue;
					index = (iz*global_dims_of_hdf5[1]+iy)*global_dims_of_hdf5[2]+ix;
				//	std::cout << "testing " << ix << " " <<  iy<< " " <<  iz << " testing at " << index << " image: " << image[index] << std::endl;
					if (image[index] != 0) {
//						std::cout << "flagtrue\n";
						flag1 = true;
						break;
					}
				}
				if(flag1)
					break;
			}
			if(flag1)
				break;
			
		}
		if (!flag1) {
			count++;
			std::cout << "no elemen for node: " << x << " " << y << " " << z << std::endl;
		}
	}
	std::cout << count << " bad BC found\n";
    return 0;
}
