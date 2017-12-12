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

//g++ -Wall -lhdf5 -o AsciiToHDF5 AsciiToHDF5.cpp

//Converts a .txt image into a hdf5 image and set the parameter, that are
//needed for the simulation. The height of the image is set to 1mm.
//
//Two example boundary condition are implemented


#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include "hdf5.h"

class BC
{
    public:
        virtual ~BC(){}
        virtual bool operator()(short, short, short, short, double &) =0;
};


class BCTopDown: public BC {
    public:
        BCTopDown(short iz, double disp): _iz(iz), _disp(disp)
    {
    }

        ~BCTopDown() {}

        bool operator()(short x, short y, short z, short d,  double &v) {
            v =0;
            //on bottom plate x,y=0 and  z = disp
            if ( z == 0) {
                if ( d == 2) {
                    v = _disp;
                }
                return true;
            }
            //on top plate all BCs are fixed 
            if ( z == _iz)
                //if (d == 2) // only z is fixed
                return true;
            return false;
        }
        short _iz;
        double _disp;

};

class BCTopDown_dispfromtop: public BC {
    public:
        BCTopDown_dispfromtop(short iz, double disp): _iz(iz), _disp(disp)
    {
    }

        ~BCTopDown_dispfromtop() {}

        bool operator()(short x, short y, short z, short d,  double &v) {
            v = 0;
            /*
               if ( z == 0 && x == 0 && y == 0 )
               return true;
               if ( z == 0 && x == _iz/2  && y == 0 && d == 1 )
               return true;
               */
            //on bottom plate x,y=0 and  z = 0
            //if ( z == 0 && d == 2 ) {
            if ( z == 0) {
                return true;
            }
            //on top plate all BCs are fixed 
            if ( z == _iz ) {
                if ( d == 2 ) {
                    v = -_disp;
                    return true;
                }
                return true;
            }
            return false;
        }
        short _iz;
        double _disp;
        };

        struct bcitem {
            bcitem(short lx=0, short ly=0, short lz=0, short ld=0):x(lx), y(ly), z(lz), d(ld)
            {
            }
            short x,y,z,d;

            bool operator<(const bcitem &rhs) const {
                if (z < rhs.z)
                    return true;
                if (z == rhs.z) {
                    if (y < rhs.y)
                        return true;
                    if ( y == rhs.y) {
                        if (x < rhs.x)
                            return true;
                        if (x == rhs.x && d < rhs.d)
                            return true;
                    }
                }
                return false;
            }
        };

        int WriteScalar(const std::string name, hid_t type, const void* data, hid_t group)
        {
            H5E_auto2_t old_func; 
            void *old_client_data;
            H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
            H5Eset_auto( H5E_DEFAULT, NULL, NULL );
            hid_t dataspace = H5Screate(H5S_SCALAR);
            hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);

            if (dataset < 0)
                dataset = H5Dcreate(group, name.c_str(), type, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            H5Dwrite( dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

            H5Dclose(dataset);
            H5Sclose(dataspace);
            H5Eset_auto( H5E_DEFAULT, old_func, old_client_data );

            return 1;
        }

        int main (int argc, char * argv[])
        {
            if (argc < 3) {
                std::cout << "usgage " << argv[0] << " asciifile hdf5file" << std::endl;
                return 1;
            }
            char * asciifilename = argv[1];
            char * hdf5filename = argv[2];

            std::fstream ifst(asciifilename, std::fstream::in);

            // Read in sizes
            hsize_t imagesize[3];
            ifst >> imagesize[0];
            ifst >> imagesize[1];
            ifst >> imagesize[2];

            double vsize = ((double)1)/imagesize[2];
            double disp = vsize*imagesize[2]*0.1; // compression 10% of the whole size
            BCTopDown_dispfromtop bc(imagesize[2],disp);
            std::map<bcitem, double> bcmap;

            //map for boundary condition

            //Read in materials
            int num_mat;
            ifst >> num_mat;
            float *elas = new float[num_mat];
            for( int i=0; i < num_mat; i++) {
                ifst >> elas[i];
            }

            std::cout << "Image has size of";
            for (int i=0; i < 3; i++)
                std::cout << " " << imagesize[i];
            std::cout << std::endl; 

            long int slice = ((long)imagesize[0])*imagesize[1];
            std::cout << "a slice has " << slice << " bytes\n";
            //int chunksize = imagesize[2]/10;
            unsigned long chunks =  ((((long)imagesize[0])*imagesize[1]*imagesize[2]-1)/ (1L<< 26))+1;
            if (chunks > imagesize[2])
                chunks = imagesize[2];
            unsigned long chunksize = (imagesize[2]-1)/chunks+1;

            std::cout << "using: " << (imagesize[2]-1)/chunksize +1  << " chunks with with size " << chunksize << " = "<<  slice*chunksize << " Bytes\n";

            float *chunkimage = new float[slice*chunksize];
            short int voxel;

            std::string datasetname("Image");


            //Suppresses error messages. If we can't create we open
            H5E_auto2_t old_func;
            void *old_client_data;
            H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
            H5Eset_auto( H5E_DEFAULT, NULL, NULL );
            // open file
            hid_t file;
            hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
            file = H5Fcreate(hdf5filename, H5F_ACC_EXCL, H5P_DEFAULT, plist);
            if (file < 0)
                file = H5Fopen(hdf5filename, H5F_ACC_RDWR, plist);
            H5Pclose(plist);

            hid_t group;
            group =  H5Gopen(file, "Image_Data", H5P_DEFAULT);
            if (group < 0)
                group = H5Gcreate(file, "Image_Data",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            // create 3d data size with current and maximal size imagesize
            hsize_t hdf_imagesize[3];
            hdf_imagesize[0] = imagesize[2]; 
            hdf_imagesize[1] = imagesize[1]; 
            hdf_imagesize[2] = imagesize[0]; 
            hid_t dataspace =  H5Screate_simple (3, hdf_imagesize, NULL);

            hid_t dataset = H5Dopen(group, datasetname.c_str(), H5P_DEFAULT);
            if (dataset < 0) 
                dataset = H5Dcreate(group, datasetname.c_str(), H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            hsize_t my_dims[3];
            hsize_t offset[3];
            hsize_t hdfdims_my_dims[3];
            hsize_t hdfdims_offset[3];
            hsize_t *stride = NULL; // No Stride -> continous junk
            hid_t memspace;
            herr_t status;

            int j;

            //turns on
            H5Eset_auto( H5E_DEFAULT, old_func, old_client_data );

            for (unsigned int i =0; i < imagesize[2];) {

                offset[0] = 0;
                offset[1] = 0;
                offset[2] = i;

                my_dims[0] = imagesize[0];
                my_dims[1] = imagesize[1];
                my_dims[2] = chunksize; 

                if (my_dims[2] + offset[2] > imagesize[2])
                    my_dims[2] = imagesize[2] -offset[2];

                std::cout << "writing " << my_dims[2] << " slices" << std::endl;

                for(unsigned long int z = 0; z < my_dims[2]; z++) {
                    for(unsigned long int y = 0; y < imagesize[1]; y++) {
                        for(unsigned long int x = 0; x < imagesize[0]; x++) { 
                            ifst >> voxel;
                            j = (z*imagesize[1]+y)*imagesize[0]+x; 
                            chunkimage[j] = elas[voxel];
                            if (chunkimage[j] != 0)
                            {
                                for (long int zi =0; zi < 2; zi++) {
                                    for (long int yi =0; yi < 2; yi++) {
                                        for (long int xi =0; xi < 2; xi++) {
                                            for (long int d = 0; d <3; d++) {
                                                double disp =0;
                                                if (bc(x+xi,y+yi,z+zi+i, d, disp))
                                                    bcmap[bcitem(x+xi,y+yi,z+i+zi,d)] = disp;
                                            }
                                        }
                                    }
                                }
                            }

                        }
                    }
                }

                //change row major
                hdfdims_offset[2] = offset[0];
                hdfdims_offset[1] = offset[1];
                hdfdims_offset[0] = offset[2];

                hdfdims_my_dims[2] = my_dims[0];
                hdfdims_my_dims[1] = my_dims[1];
                hdfdims_my_dims[0] = my_dims[2];

                status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hdfdims_offset, stride,  hdfdims_my_dims, NULL);
                if (status <0)
                    std::cout << "Error in Selection";
                memspace =  H5Screate_simple (3, hdfdims_my_dims, NULL);
                plist = H5P_DEFAULT;
                H5Dwrite( dataset, H5T_NATIVE_FLOAT,  memspace , dataspace, H5P_DEFAULT, chunkimage);
                H5Sclose(memspace);

                i += my_dims[2];
            }

            H5Sclose(dataspace);
            H5Dclose(dataset);
            std::vector<short> bccoords;
            std::vector<float> bcdisp;
            for(std::map<bcitem, double>::iterator it =  bcmap.begin(); it != bcmap.end(); ++it) {
                bccoords.push_back(it->first.z);
                bccoords.push_back(it->first.y);
                bccoords.push_back(it->first.x);
                bccoords.push_back(it->first.d);
                bcdisp.push_back(it->second);
            }
            hsize_t dims[2];
            dims[1] = 4;
            dims[0] = bcdisp.size();
            dataspace = H5Screate_simple(2,dims, NULL);
            dataset = H5Dcreate (group, "Fixed_Displacement_Coordinates", H5T_NATIVE_USHORT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dataset, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bccoords[0]);
            H5Sclose(dataspace);
            H5Dclose(dataset);

            dims[0] = bcdisp.size();
            dataspace = H5Screate_simple(1,dims, NULL);
            dataset = H5Dcreate (group, "Fixed_Displacement_Values", H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &bcdisp[0]);
            H5Sclose(dataspace);
            H5Dclose(dataset);

            double tmp = 0.3;
            WriteScalar("Poisons_ratio",  H5T_NATIVE_DOUBLE, &tmp, group);
            tmp = ((double)1)/imagesize[2];
            WriteScalar("Voxelsize",  H5T_NATIVE_DOUBLE, &tmp, group);


            H5Gclose(group);
            H5Fclose(file);
            delete[] chunkimage;
            delete[] elas;

            return 0;
        }
