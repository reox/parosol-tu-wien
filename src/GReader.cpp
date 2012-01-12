/*
 * ParFE: a micro-FE solver for trabecular bone modeling
 * Copyright (C) 2006,2011, Uche Mennel and Marzio Sala, Cyril Flaig
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

#include "GReader.hpp"
#include <cstring>

//#define H5_HAVE_PARALLEL
//#define __QK_USER__

//C functions
C_ASCII_GReader::C_ASCII_GReader(MPI_Comm a_comm)
    : file(NULL), comm(a_comm), mpi_rank(0), mpi_size(1)
{
#ifdef HAVE_MPI
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
#endif
}

C_ASCII_GReader::C_ASCII_GReader(const std::string& a_filename, MPI_Comm a_comm)
  : filename(a_filename), file(NULL), comm(a_comm), mpi_rank(0), mpi_size(1)   //save a copy of the filename
{
#ifdef HAVE_MPI
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
#endif
  Open();
}

C_ASCII_GReader::~C_ASCII_GReader()
{
  Close();
}

bool C_ASCII_GReader::Open()
{
  file = fopen(filename.c_str(), "r");
  return file;
}

bool C_ASCII_GReader::Open(const std::string& a_filename)
{
  filename = a_filename;
  return Open();
}

void C_ASCII_GReader::Close()
{
  fclose(file);
}

int C_ASCII_GReader::Read(long& i)
{
  char buffer[64];
  char c;
  long count=0, diff;
  char *last = buffer;
  while ( isspace( c = fgetc(file)) );
  buffer[count++] = c;
  while ( !isspace(buffer[count++] = fgetc(file)) && count < 64);
  buffer[count] = '\0';
  i = strtol(buffer, &last, 10);
  diff = last - buffer;
  fseek(file, diff - count, SEEK_CUR);
  return diff;
}

int C_ASCII_GReader::Read(double& d)
{
  char buffer[64];
  long count=0, diff;
  char *last = buffer;
  char c;
  while ( isspace( c = fgetc(file)) );
  buffer[count++] = c;
  while ( !isspace(buffer[count++] = fgetc(file)) && count < 64);
  buffer[count] = '\0';
  d = strtod(buffer, &last);
  diff = last - buffer;
  fseek(file, diff - count, SEEK_CUR);
  return diff;
}

int C_ASCII_GReader::Read(std::string& s)
{
  char buffer[256];
  long count=0;
  char c;
  while ( isspace( c = fgetc(file)) );
  buffer[count++] = c;
  while ( !isspace(buffer[count++] = fgetc(file)) && count < 256);
  ungetc(buffer[--count], file);
  buffer[count] = '\0';
  s = buffer;
  return count;
}

bool C_ASCII_GReader::operator!()
{
  return !file;
}

int C_ASCII_GReader::Skip(char comment)
{
  char c;
  char buffer[1024];
  long count = 0;
  while (isspace( c = fgetc(file)) );
  while (c == comment) {
      fgets(buffer, 1024, file);
      while (isspace( c = fgetc(file)) );
      ++count;
  }
  ungetc(c, file);
  return count;
}

bool C_ASCII_GReader::Select(const std::string& s)
{
  //nothing to do here
  return true;
}

bool C_ASCII_GReader::Seek(const std::string& s)
{
  if (marks.count(s) == 0)
    fgetpos(file, &marks[s]);
  return (fsetpos(file, &marks[s]) == 0);
}

//CPP functions
CPP_ASCII_GReader::CPP_ASCII_GReader(MPI_Comm a_comm)
    :comm(a_comm), mpi_rank(0), mpi_size(1)
{
#ifdef HAVE_MPI
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
#endif
}

CPP_ASCII_GReader::CPP_ASCII_GReader(const std::string& a_filename, MPI_Comm a_comm)
    : filename(a_filename), comm(a_comm), mpi_rank(0), mpi_size(1)
{
#ifdef HAVE_MPI
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
#endif
  Open();
}

CPP_ASCII_GReader::~CPP_ASCII_GReader()
{
  Close();
}

bool CPP_ASCII_GReader::Open()
{
  file.open(filename.c_str());
  return file.good();
}

bool CPP_ASCII_GReader::Open(const std::string& a_filename)
{
  filename = a_filename;
  return Open();
}

void CPP_ASCII_GReader::Close()
{
  file.close();
}

bool CPP_ASCII_GReader::operator!()
{
  return (!file);
}


int CPP_ASCII_GReader::Skip(char comment)
{
  std::string s;
  unsigned char c;
  long count = 0;

  file >> c;
  while (c == comment) { 
    std::getline( file, s );
    file >> c;
    ++count;
  }
  file.putback(c);

  return count;
}


bool CPP_ASCII_GReader::Select(const std::string& s)
{
  //nothing to do here
  return true;
}


bool CPP_ASCII_GReader::Seek(const std::string& s)
{
  if (marks.count(s) == 0)
    marks[s] = file.tellg();
  return file.seekg(marks[s]).good();
}


HDF5_GReader::HDF5_GReader(MPI_Comm a_comm)
    : file(0), plist(H5Pcreate(H5P_FILE_ACCESS)), index(0),  comm(a_comm), file_info(MPI_INFO_NULL), mpi_rank(0), mpi_size(1)
{
#ifdef HAVE_MPI
  MPI_Info_create(&file_info);
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
#endif 
//Initialize driver specific properties
#ifdef H5_HAVE_PARALLEL
  H5Pset_fapl_mpio(plist, comm, file_info);
#endif
  offset[0] = 0;
  offset[1] = 0;
  stride[0] = 1;
  stride[1] = 1;
}


HDF5_GReader::HDF5_GReader(const std::string& a_filename, MPI_Comm a_comm)
    : filename(a_filename), file(0), plist(H5Pcreate(H5P_FILE_ACCESS)), index(0), comm(a_comm), file_info(MPI_INFO_NULL), mpi_rank(0), mpi_size(1)
{
#ifdef HAVE_MPI
  MPI_Info_create(&file_info);
  MPI_Comm_size(comm, &mpi_size);
  MPI_Comm_rank(comm, &mpi_rank);
#ifdef __QK_USER__
  //Optimization parameters needed if running on the CRAY XT3
  MPI_Info_set(file_info, "cb_config_list", "*:*");
  MPI_Info_set(file_info, "cb_buffer_size", "1048576");
  MPI_Info_set(file_info, "romio_cb_read", "enable");
#endif
#endif
 //Initialize driver specific properties
#ifdef H5_HAVE_PARALLEL
  H5Pset_fapl_mpio(plist, comm, file_info);
#endif
  stride[0] = 1;
  stride[1] = 1;
  offset[0] = 0;
  offset[1] = 0;
  Open();
}


HDF5_GReader::~HDF5_GReader()
{
  Close();
#ifdef HAVE_MPI
  MPI_Info_free(&file_info);
#endif
}


bool HDF5_GReader::Open(const std::string& a_filename)
{
  filename = a_filename;
  return Open();
}


bool HDF5_GReader::Open()
{
  //Suppress error messges.
  //H5Eset_auto( NULL, NULL );
  file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist);
  H5Pclose(plist);
  group = H5Gopen(file, "/", H5P_DEFAULT);
#ifdef SEQUENTIAL_HDF5
  plist = H5P_DEFAULT;
#else
  //but collective I/O is more effective
  plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
#endif
        
  H5Eset_auto( H5E_DEFAULT, NULL, NULL );
  return (group >= 0);
}


void HDF5_GReader::Close()
{

  H5Gclose(group);
  H5Fclose(file);
}


bool HDF5_GReader::operator!()
{
  return (file<0);
}


int HDF5_GReader::Skip(char comment)
{
  //nothing to do here
  return 0;
}


bool HDF5_GReader::Select(std::string s)
{
  hid_t temp;
  //Test group
  temp = group;
  group = H5Gopen(temp, s.c_str(), H5P_DEFAULT);
  if (group < 0) {
    return false;
  }
  H5Gclose(temp);
  return true;
}


int HDF5_GReader::Read(const std::string& name, hid_t type, void* data, long num_global_elems, long num_my_elems, long elem_size, long my_offset)
{
  hsize_t len = num_my_elems*elem_size;
  hsize_t my_dims[2];
  my_dims[0] = num_my_elems;
  my_dims[1] = elem_size;
  
  hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
  if (dataset < 0)
    return -1;
  
  //Get dataspace of the dataset.
  hid_t dataspace = H5Dget_space(dataset);
  offset[0] = my_offset;
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, stride,  my_dims, NULL );

  hid_t memspace = H5Screate_simple( 1, &len, NULL );
  if (memspace <= 0) {
    hsize_t one=1;
    memspace = H5Screate_simple( 1, &one, NULL );
    H5Sselect_none(memspace);
  }
  
  H5Dread(dataset, type, memspace, dataspace, plist, data);
  H5Sclose(dataspace);
  H5Sclose(memspace);
  H5Dclose(dataset);
  return len;
}


int HDF5_GReader::Read(const std::string& name, hid_t type, void* data)
{
  hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
  if (dataset < 0)
    return -1;
  H5Dread(dataset, type, H5S_ALL, H5S_ALL, plist, data);
  H5Dclose(dataset);
  return 1;
}


int HDF5_GReader::Read(const std::string& name, std::string& data)
{
  const long buf_size = 80;
  int res = 0;
  char c_str[buf_size];
  if (mpi_rank == 0)
    res = Read(name, getNativeType(data), c_str);
  MPI_Bcast(c_str, buf_size, MPI_CHAR, 0, comm);
  data = c_str;
  return res;
}

int HDF5_GReader::Read(const std::string& name, hid_t type, void* data, hsize_t* my_offset, hsize_t *count, const int dims)
{
  herr_t status;
  hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
  if (dataset < 0)
    return -1;
  
  //Get dataspace of the dataset.
  hid_t dataspace = H5Dget_space(dataset);
  
// std::cout << "reading from " << my_offset[0] << " " << my_offset[1] << " " << my_offset[2] << "\n";
//  std::cout << "to " << count[0] << " " << count[1] << " " << count[2] << "\n";

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, my_offset, NULL,  count, NULL );
  
   if (status <0) {
      std::cout << "Error in Selection";
      return -1;
   }
   
  hid_t memspace = H5Screate_simple( dims, count, NULL );
  if (memspace <= 0) {
    std::cout << "Error in memspace";
    return -1;
  }
  
  H5Dread(dataset, type, memspace, dataspace, plist, data);

  H5Sclose(dataspace);
  H5Sclose(memspace);
  H5Dclose(dataset);
  return 1;
}

int HDF5_GReader::GetSizeOfDataset(const std::string& name, hsize_t* size, const int number_dims)
{
  if (mpi_rank == 0) {
    hid_t dataset = H5Dopen(group, name.c_str(), H5P_DEFAULT);
    if (dataset < 0)
      return -1;
  
    // Read the dimensions
    hsize_t* max_dims = new hsize_t[number_dims];
    hid_t dataspace = H5Dget_space(dataset);
    H5Sget_simple_extent_dims(dataspace, size, max_dims);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    delete[] max_dims;
  }
  MPI_Bcast(size, 3*sizeof(hsize_t), MPI_BYTE, 0, MPI_COMM_WORLD);
  return 1;
}

hid_t HDF5_GReader::getNativeType(int, hsize_t i) {
  if (i == 1)
    return H5T_NATIVE_INT;
  return H5Tarray_create(H5T_NATIVE_INT, 1, &i);
}

hid_t HDF5_GReader::getNativeType(double, hsize_t i) {
  if (i == 1)
    return H5T_NATIVE_DOUBLE;
  return H5Tarray_create(H5T_NATIVE_DOUBLE, 1, &i);
}

hid_t HDF5_GReader::getNativeType(float, hsize_t i) {
  if (i==1)
    return H5T_NATIVE_FLOAT; 
  return H5Tarray_create(H5T_NATIVE_FLOAT, 1, &i);
}
hid_t HDF5_GReader::getNativeType(short, hsize_t i) {
  if (i==1)
    return H5T_NATIVE_SHORT; 
  return H5Tarray_create(H5T_NATIVE_SHORT, 1, &i);
}

hid_t HDF5_GReader::getNativeType(const std::string&, hsize_t i) {
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, 32);
  if (i==1)
    return type; 
  return H5Tarray_create(type, 1, &i);
}

void HDF5_GReader::copy(void* dest, const void* src, long nblocks, long block_size, long stride)
{
  for(long i=0; i<nblocks; ++i) {
    memcpy(dest, src, block_size);
    dest = (void*)((unsigned long)(dest)+block_size);
    src = (void*)((unsigned long)(src)+stride);
  }
}
