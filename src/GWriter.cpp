/*
 * ParFE: a micro-FE solver for trabecular bone modeling
 * Copyright (C) 2006,2001, Uche Mennel and Marzio Sala, Cyril Flaig
 * 
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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  
 * 02110-1301, USA.
 */

#include "GWriter.hpp"
#include <unistd.h> //fsync
#include <cstring>

//#define H5_HAVE_PARALLEL
//#define __QK_USER__

//C functions
C_ASCII_GWriter::C_ASCII_GWriter(MPI_Comm a_comm)
  : file(NULL), comm(a_comm)
{
#ifdef HAVE_MPI
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
#endif
}

C_ASCII_GWriter::C_ASCII_GWriter(const std::string& a_filename, MPI_Comm a_comm)
  : filename(a_filename), comm(a_comm)   //save a copy of the filename
{
#ifdef HAVE_MPI
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
#endif
  //Truncate file
  Open();
  Close();
}

C_ASCII_GWriter::~C_ASCII_GWriter()
{
}

bool C_ASCII_GWriter::Open()
{
  file = fopen(filename.c_str(), "w");
  return (file > 0);
}

bool C_ASCII_GWriter::Open(const std::string& a_filename)
{
  filename = a_filename;
  return Open();
}

bool C_ASCII_GWriter::Append()
{
  file = fopen(filename.c_str(), "a");
  return file;
}

void C_ASCII_GWriter::Close()
{
  fclose(file);
}

bool C_ASCII_GWriter::operator!()
{
  return !file;
}

bool C_ASCII_GWriter::Select(const std::string& s)
{
  //nothing to do here
  return true;
}


bool C_ASCII_GWriter::Seek(const std::string& s)
{
  if (marks.count(s) == 0)
    fgetpos(file, &marks[s]);
  return (fsetpos(file, &marks[s]) == 0);
}

void C_ASCII_GWriter::Flush()
{
  fflush(file);
  fsync(fileno(file));
}

int C_ASCII_GWriter::Write(const long & i)
{
  return fprintf(file, "%ld ", i);
}

int C_ASCII_GWriter::Write(const double& d)
{
  return fprintf(file, "%e ", d);
}

int C_ASCII_GWriter::Write(const char* s)
{
  return fprintf(file, "%s ", s);
}

int C_ASCII_GWriter::Write(const std::string& s)
{
  return fprintf(file, "%s ", s.c_str());
}


//CPP functions
CPP_ASCII_GWriter::CPP_ASCII_GWriter(MPI_Comm a_comm)
  : comm(a_comm), mpi_size(0), mpi_rank(0)
{
#ifdef HAVE_MPI
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
#endif
}


CPP_ASCII_GWriter::CPP_ASCII_GWriter(const std::string& a_filename, MPI_Comm a_comm)
  : filename(a_filename), comm(a_comm), mpi_size(0), mpi_rank(0)
{
#ifdef HAVE_MPI
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
#endif
  //Truncate file
  Open();
  Close();
}


CPP_ASCII_GWriter::~CPP_ASCII_GWriter()
{}

bool CPP_ASCII_GWriter::Open()
{
  file.open(filename.c_str());
  return file.good();
}

bool CPP_ASCII_GWriter::Open(const std::string& a_filename)
{
  filename = a_filename;
  return Open();
}

bool CPP_ASCII_GWriter::Append()
{
  file.open(filename.c_str(), std::ios::app);
  return file.good();
}

void CPP_ASCII_GWriter::Close()
{
  file.close();
}


bool CPP_ASCII_GWriter::operator!()
{
  return (!file);
}


void CPP_ASCII_GWriter::Flush()
{
  file.flush();
}


bool CPP_ASCII_GWriter::Select(const std::string& s)
{
  //nothing to do here
  return true;
}


bool CPP_ASCII_GWriter::Seek(const std::string& s)
{
  if (marks.count(s) == 0)
    marks[s] = file.tellp();
  return file.seekp(marks[s]).good();
}


HDF5_GWriter::HDF5_GWriter(MPI_Comm a_comm)
    : groupname("/"), comm(a_comm), file_info(MPI_INFO_NULL), file(-1), plist(H5Pcreate(H5P_FILE_ACCESS)), mpi_size(0), mpi_rank(0), index(0)
{
#ifdef HAVE_MPI
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Info_create(&file_info);
#ifdef __QK_USER__
  //Optimization parameters needed if running on the CRAY XT3
  MPI_Info_set(file_info, "cb_config_list", "*:*");
  MPI_Info_set(file_info, "cb_buffer_size", "1048576");
  MPI_Info_set(file_info, "romio_cb_write", "enable");
#endif
#endif
#ifdef H5_HAVE_PARALLEL
  //Initialize driver specific properties
  H5Pset_fapl_mpio(plist, comm, file_info);
#endif
  offset[0] = 0;
  offset[1] = 0;
  stride[0] = 1;
  stride[1] = 1;
}


HDF5_GWriter::HDF5_GWriter(const std::string& a_filename, MPI_Comm a_comm)
    : filename(a_filename), groupname("/"), comm(a_comm),  file_info(MPI_INFO_NULL), file(-1),  plist(H5Pcreate(H5P_FILE_ACCESS)), mpi_size(0), mpi_rank(0), index(0)
{
#ifdef HAVE_MPI
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Info_create(&file_info);
#ifdef __QK_USER__
  //Optimization parameters needed if running on the CRAY XT3
  MPI_Info_set(file_info, "cb_config_list", "*:*");
  //MPI_Info_set(file_info, "cb_buffer_size", "67108864");
  MPI_Info_set(file_info, "romio_cb_write", "enable");
#endif
#endif
#ifdef H5_HAVE_PARALLEL
  //Initialize driver specific properties
  H5Pset_fapl_mpio(plist, comm, file_info);
#endif
  stride[0] = 1;
  stride[1] = 1;
  offset[0] = 0;
  offset[1] = 0;
#ifndef SEQUENTIAL_HDF5
  Open();
#endif
}


HDF5_GWriter::~HDF5_GWriter()
{
#ifndef SEQUENTIAL_HDF5
  Close();
#endif
#ifdef HAVE_MPI
  MPI_Info_free(&file_info);
#endif
}


bool HDF5_GWriter::Open(const std::string& a_filename)
{
  filename = a_filename;
  return Open();
}


bool HDF5_GWriter::Open()
{
  //Suppress error messges.
  H5Eset_auto( H5E_DEFAULT, NULL, NULL );
  file = H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, plist);
  if (file < 0)
    file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist);
  H5Pclose(plist);
  group = H5Gopen(file, groupname.c_str(), H5P_DEFAULT);
  if (group < 0)
    group = H5Gcreate(file, groupname.c_str(),H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  //independent I/O is the default
#ifdef SEQUENTIAL_HDF5
  plist = H5P_DEFAULT;
#else
//but collective I/O is more effective
  plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
#endif
  return (group >= 0);
}


void HDF5_GWriter::Close()
{
#ifndef SEQUENTIAL_HDF5
  H5Pclose(plist);
#endif
  H5Gclose(group);
  H5Fclose(file);
}


bool HDF5_GWriter::operator!()
{
  return (file < 0);
}


bool HDF5_GWriter::Select(const std::string& s)
{
  groupname = s;  
#ifndef SEQUENTIAL_HDF5
  hid_t temp;
  //Test group
  temp = group;
  group = H5Gopen(temp, groupname.c_str(), H5P_DEFAULT);
  if (group < 0) {
    group = H5Gcreate1(temp, groupname.c_str(),0);
  }
  H5Gclose(temp);
#endif
  return (group >= 0);
}


int HDF5_GWriter::Write(const std::string& name, hid_t type, const void* data, long num_global_elems, long num_my_elems, long elem_size, long my_offset)
{

  hsize_t len = num_my_elems*elem_size;
  hsize_t dims[2];
  dims[0]= num_global_elems;
  dims[1]= elem_size;
  hsize_t my_dims[2];
  my_dims[0] = num_my_elems;
  my_dims[1] = elem_size;
  
#ifdef SEQUENTIAL_HDF5
  long i=0;
  MPI_Status status;
  //synchroniztion
  MPI_Barrier(comm);
  if (mpi_rank != 0) 
    MPI_Recv(&i, 1, MPI_INT, mpi_rank-1, 0, comm, &status);
  Open();
#endif
  hid_t dataspace = H5Screate_simple(2, dims, NULL);

  //hid_t memspace = H5Screate_simple( 1, &len, NULL );
  hid_t memspace = H5Screate_simple( 2, my_dims, NULL );
  if (memspace <= 0) {
      hsize_t one=1;
      memspace = H5Screate_simple( 1, &one, NULL );
      H5Sselect_none(memspace);
  }

  hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
  if (dataset < 0)
    dataset = H5Dcreate(group, name.c_str(), type, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  offset[0] = my_offset;
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, stride,  my_dims, NULL);

  H5Dwrite( dataset, type, memspace, dataspace, plist, data );

  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

#ifdef SEQUENTIAL_HDF5
  Close();
  //send my_offset+num_my_elems to next CPU
  i = my_offset+num_my_elems;
  if (mpi_rank != mpi_size-1)
    MPI_Send(&i, 1, MPI_INT, mpi_rank+1, 0, comm);
#endif

  return len;
}

int HDF5_GWriter::Write(const std::string& name, hid_t type, const void* data)
{
    
#ifdef SEQUENTIAL_HDF5
  MPI_Barrier(comm);
  if (mpi_rank == 0) {
    Open();
#endif
    hid_t dataspace = H5Screate(H5S_SCALAR);
    hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
	
    if (dataset < 0)
      dataset = H5Dcreate(group, name.c_str(), type, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      
    H5Dwrite( dataset, type, H5S_ALL, H5S_ALL, plist, data);

    H5Dclose(dataset);
    H5Sclose(dataspace);
#ifdef SEQUENTIAL_HDF5
    Close();
  }
#endif

  return 1;
}


int HDF5_GWriter::Write(const std::string& name, const std::string& data)
{
  return Write(name, getNativeType(data), data.c_str());
}

hid_t HDF5_GWriter::getNativeType(long, hsize_t i) {
  if (i == 1)
    return H5T_NATIVE_LONG;
  return H5Tarray_create(H5T_NATIVE_LONG, 1, &i);
}

hid_t HDF5_GWriter::getNativeType(short, hsize_t i) {
  if (i == 1)
    return H5T_NATIVE_SHORT;
  return H5Tarray_create(H5T_NATIVE_SHORT, 1, &i);
}

hid_t HDF5_GWriter::getNativeType(int, hsize_t i) {
  if (i == 1)
    return H5T_NATIVE_INT;
  return H5Tarray_create(H5T_NATIVE_INT, 1, &i);
}

hid_t HDF5_GWriter::getNativeType(double, hsize_t i) {
  if (i == 1)
    return H5T_NATIVE_DOUBLE;
  return H5Tarray_create(H5T_NATIVE_DOUBLE, 1, &i);
}

hid_t HDF5_GWriter::getNativeType(float, hsize_t i) {
  if (i==1)
    return H5T_NATIVE_FLOAT; 
  return H5Tarray_create(H5T_NATIVE_FLOAT, 1, &i);
}

hid_t HDF5_GWriter::getNativeType(const std::string&, hsize_t i) {
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, 256);
  if (i==1)
    return type; 
  return H5Tarray_create(type, 1, &i);
}

void HDF5_GWriter::copy(void* dest, const void* src, long nblocks, long block_size, long stride)
{
  for(long i=0; i<nblocks; ++i) {
    memcpy(dest, src, block_size);
    dest = (void*)((unsigned long)(dest)+stride);
    src = (void*)((unsigned long)(src)+block_size);
  }
}
