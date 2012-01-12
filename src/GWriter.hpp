/*
 * ParFE: a micro-FE solver for trabecular bone modeling
 * Copyright (C) 2006, 2011, Uche Mennel and Marzio Sala, Cyril Flaig
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

#ifndef _GWRITER_HPP_
#define _GWRITER_HPP_

#define HAVE_MPI

#include <cstdlib> // GCC 4.3
#include<iostream>
#include<fstream>
#include<cstdio>
#ifdef HAVE_MPI
#include <mpi.h>
#else
#define MPI_Comm void*
#define MPI_Info void*
#define MPI_COMM_WORLD 0
#define MPI_INFO_NULL 0
#endif
#include <hdf5.h>
#include <map>

#ifndef H5_HAVE_PARALLEL
#ifdef HAVE_MPI
#define SEQUENTIAL_HDF5
#endif
#endif


//! A generic class to write scalars and (block) vectors sequentially into an ASCII file by using only C library functions.

/*! Writing into file concurrently requires synchronization of the processes involved into the writing process.
    The class C_ASCII_GWriter provides a safe and generic methods to allow multiple processes to write vector/scalar data into the
    same file.
    Only C library I/O functions are invoked for the case that they perform best on the target architecture.
    To use this class in parallel environment, an MPI implementation must be available.
*/

class C_ASCII_GWriter
{
 public:
  //! C_ASCII_GWriter constructor.
  /*! Creates a C_ASCII_GWriter object without performing any I/O operation.

    \param comm
    (In) optional MPI communicator containing information on the parallel environment.
	    
    \return Pointer to a C_ASCII_GWriter object.

  */ 
  C_ASCII_GWriter(MPI_Comm comm = MPI_COMM_WORLD);


  //! C_ASCII_GWriter constructor.
  /*! Creates a C_ASCII_GWriter object and opens the file with the given filename for writing.
      If the file already exits, it will be overwritten

    \param filename
    (In) Path name of the file to be written.
    The file will be overwritten, if already existent.  

    \param comm
    (In) Optional MPI communicator containing information on the parallel environment.
	    
    \return Pointer to a C_ASCII_GWriter object.

  */ 
  C_ASCII_GWriter(const std::string& filename, MPI_Comm comm = MPI_COMM_WORLD);

  //! C_ASCII_GWriter destructor.
  /*! Closes any open files.
  */ 
  ~C_ASCII_GWriter();

  //! Opens the file with the previously given filename.
   /*!
    \return true, only if successful.
  */
  bool Open();

  //! Opens the file with the given filename.
  /*!
    \param filename
    (In) Path name of the file to be written.
    The file will be overwritten, if already existent.

    \return true, only if successful.
  */
  bool Open(const std::string& filename);

  //! Closes any open files
  void Close();

  //! Generic function to write a scalar value.
   /*!
    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
 
    \param scalar
    (In) The scalar value to be written.
    Supported types, Integers, floating point numbers, strings. 

    \return number of records written.
  */
  template<typename T> int Write(const std::string& name, const T& scalar);

  //! Generic function to write a block vector.

  /*! A block will correspond to a line in the ASCII file.
    \param name
            (In) The name of a section in the file.
	    Each section in the file is associated with a file position.
	    If there exists already a section with this name, the function will jump to the corresponding file position.
 
    \param vec
            (In) Pointer to the block vector to be written.
	    Supported types of vector elements, Integers, floating point numbers, strings.

    \param num_global_elems 
            (In) Global number of vector blocks.

    \param num_my_elems
            (In) Local number of vector blocks.

    \param elem_size
            (In) Size of a vector block.
	    
    \param my_offset
            (In) Defines the block at which the writing operation should start.
	    In C_ASCII_GWriter, all processes must write consecutive blocks.

    \return number of records written.
  */
  template<typename T> int Write(const std::string& name, const T* vec, long num_global_elems, long num_my_elems, long elem_size, long my_offset);

  //! Generic function to write a block multityped vector.

  /*! A block will correspond to a line in the ASCII file. Each block consists of two subblocks of different types.
      The size and content of the subblocks is determined by the two given block vectors, i.e. each line is a
      concatenation of a block of vector one and a block of vector two.

    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
    
    \param name1
    (In) Description of the first vector.
    This argument has no effect in C_ASCII_GWriter.
	    
    \param vec1
    (In) Pointer to the block vector of the first type to be written.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    This argument has no effect in C_ASCII_GWriter.

    \param vec2
    (In) Pointer to the block vector of the second type to be written.
    Supported types of vector elements are Integers, floating point numbers, strings.	    
    
    \param num_global_elems
    (In) Global number of vector blocks.

    \param num_my_elems
    (In) Local number of vector blocks.

    \param elem_size1
    (In) Number of elements of type of the first vector in a vector block.

    \param elem_size2
    (In) Number of elements of type of the second vector in a vector block.
	    
    \param my_offset
    (In) Defines the block at which the writing operation should start.
    In C_ASCII_GWriter, all processes must write consecutive blocks.

    \return number of records written.
  */
  template<typename T1, typename T2> int Write(const std::string& name, const std::string& name1, const T1* vec1, const std::string& name2, const T2* vec2, long num_global_elems, long num_my_elems, long elem_size1, long elem_size2, long my_offset);
  
  //! Generic function to write a block multityped vector.
  /*! A block will correspond to a line in the ASCII file. Each block consists of three subblocks of different types.
      The size and content of the subblocks is determined by the three given block vectors, i.e. each line is a.
      concatenation of a block of vector one and a block of vector two and a block of vector three.

    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
    
    \param name1
    (In) Description of the first vector.
    This argument has no effect in C_ASCII_GWriter.
	    
    \param vec1
    (In) Pointer to the block vector of the first type to be written.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    This argument has no effect in C_ASCII_GWriter.

    \param vec2
    (In) Pointer to the block vector of the second type to be written.
    Supported types of vector elements are Integers, floating point numbers, strings.	    
    
    \param name3
    (In) Description of the third vector.
    This argument has no effect in C_ASCII_GWriter.
	    
    \param vec3
    (In) Pointer to the block vector of the third type to be written.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param num_global_elems
    (In) Global number of vector blocks.

    \param num_my_elems
    (In) Local number of vector blocks.

    \param elem_size1
    (In) Number of elements of type of the first vector in a vector block.

    \param elem_size2
    (In) Number of elements of type of the second vector in a vector block.

    \param elem_size3
    (In) Number of elements of type of the third vector in a vector block.	    

    \param my_offset
    (In) Defines the block at which the writing operation should start.
    In C_ASCII_GWriter, all processes must write consecutive blocks.
    
    \return number of records written.
  */
  template<typename T1, typename T2, typename T3> int Write(const std::string& name, const std::string& name1, const T1* vec1, const std::string& name2, const T2* vec2,const std::string& name3, const T3* vec3, long  num_global_elems, long num_my_elems, long elem_size1, long elem_size2, long elem_size3, long my_offset);

  //! Indicates wether the file has been properly opened
  /*!
      \return true, only if the file is NOT opened
  */
  bool operator!();

  //! This function does nothing. It is only implemented for convenience
  /*
      \return true
  */
  bool Select(const std::string&);

 private:
  bool Seek(const std::string&);
  int Write(const long& i);
  int Write(const double& d);
  int Write(const char* s);
  int Write(const std::string& s);
  bool Append();
  void Flush();
  FILE* file;
  std::string filename;
  MPI_Comm comm;
  int mpi_size;
  int mpi_rank;
  std::map<std::string, fpos_t> marks;

};

//! A generic class to write scalars and (block) vectors sequentially into an ASCII file by using only C++ library functions.

/*! Writing into file concurrently requires synchronization of the processes involved into the writing process.
    The class C_ASCII_GWriter provides a safe and generic methods to allow multiple processes to write vector/scalar data into the
    same file.
    Only C++ library I/O functions are invoked for the case that they perform best on the target architecture.
    To use this class in parallel environment, an MPI implementation must be available.
*/

class CPP_ASCII_GWriter
{
 public:
  //! CPP_ASCII_GWriter constructor.
  /*! Creates a CPP_ASCII_GWriter object without performing any I/O operation.

    \param comm
    (In) Optional MPI communicator containing information on the parallel environment.
	    
    \return Pointer to a CPP_ASCII_GWriter object.

  */ 
  CPP_ASCII_GWriter(MPI_Comm comm= MPI_COMM_WORLD);
  //! CPP_ASCII_GWriter constructor.
  /*! Creates a CPP_ASCII_GWriter object and opens the file with the given filename for writing.
      If the file already exits, it will be overwritten

    \param filename 
    (In) Path name of the file to be written.
    The file will be overwritten, if already existent.

    \param comm
    (In) Optional MPI communicator containing information on the parallel environment.
	    
    \return Pointer to a CPP_ASCII_GWriter object.

  */ 
  CPP_ASCII_GWriter(const std::string& filename, MPI_Comm comm= MPI_COMM_WORLD);


  //! CPP_ASCII_GWriter destructor.
  /*! Closes any open files.
  */ 
  ~CPP_ASCII_GWriter();

  //! Opens the file with the given filename.
  /*!
    \return true, only if successful.
  */
  bool Open();

  //! Opens the file with the given filename.
  /*!
    \param filename
    (In) Path name of the file to be written.
    The file will be overwritten, if already existent.

    \return true, only if successful
  */
  bool Open(const std::string& filename);

  //! Closes any open files
  void Close();

  //! Generic function to write a scalar value
   /*!
     \param name
     (In) The name of a section in the file.
     Each section in the file is associated with a file position.
     If there exists already a section with this name, the function will jump to the corresponding file position
 
    \param scalar
    (In) The scalar value to be written.
    Supported types, Integers, floating point numbers, strings. 

    \return number of records written.
  */
  template<typename T> int Write(const std::string& name, const T& scalar);

  //! Generic function to write a block vector

  /*! A block will correspond to a line in the ASCII file.
    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
 
    \param vec
    (In) Pointer to the block vector to be written.
    Supported types of vector elements, Integers, floating point numbers, strings.

    \param num_global_elems
    (In) Global number of vector blocks

    \param num_my_elems
    (In) Local number of vector blocks

    \param elem_size
    (In) Size of a vector block
	    
    \param my_offset
    (In) Defines the block at which the writing operation should start.
    In CPP_ASCII_GWriter, all processes must write consecutive blocks.

    \return number of records written.
  */
  template<typename T> int Write(const std::string& name, const T* vec, long num_global_elems, long num_my_elems, long elem_size, long my_offset);

  //! Generic function to write a block multityped vector.

  /*! A block will correspond to a line in the ASCII file. Each block consists of two subblocks of different types.
      The size and content of the subblocks is determined by the two given block vectors, i.e. each line is a
      concatenation of a block of vector one and a block of vector two.

    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
    
    \param name1
    (In) Description of the first vector.
    This argument has no effect in CPP_ASCII_GWriter.
	    
    \param vec1
    (In) Pointer to the block vector of the first type to be written.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    This argument has no effect in CPP_ASCII_GWriter.

    \param vec2
    (In) Pointer to the block vector of the second type to be written.
    Supported types of vector elements are Integers, floating point numbers, strings.	    
    
    \param num_global_elems
    (In) Global number of vector blocks.

    \param num_my_elems
    (In) Local number of vector blocks.

    \param elem_size1
    (In) Number of elements of type of the first vector in a vector block.

    \param elem_size2
    (In) Number of elements of type of the second vector in a vector block.
	    
    \param my_offset
    (In) Defines the block at which the writing operation should start.
    In CPP_ASCII_GWriter, all processes must write consecutive blocks.

    \return number of records written.
  */
  template<typename T1, typename T2> int Write(const std::string& name, const std::string& name1, const T1* vec1, const std::string& name2, const T2* vec2, long num_global_elems, long num_my_elems, long elem_size1, long elem_size2, long my_offset);

  //! Generic function to write a block multityped vector.
  /*! A block will correspond to a line in the ASCII file. Each block consists of three subblocks of different types.
      The size and content of the subblocks is determined by the three given block vectors, i.e. each line is a
      concatenation of a block of vector one and a block of vector two and a block of vector three.

    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
    
    \param name1
    (In) Description of the first vector.
    This argument has no effect in CPP_ASCII_GWriter.
	    
    \param vec1
    (In) Pointer to the block vector of the first type to be written.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    This argument has no effect in CPP_ASCII_GWriter.
    
    \param vec2
    (In) Pointer to the block vector of the second type to be written.
    Supported types of vector elements are Integers, floating point numbers, strings.    
    
    \param name3
    (In) Description of the third vector.
    This argument has no effect in CPP_ASCII_GWriter.
	    
    \param vec3
    (In) Pointer to the block vector of the third type to be written.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param num_global_elems
    (In) Global number of vector blocks

    \param num_my_elems
    (In) Local number of vector blocks

    \param elem_size1
    (In) Number of elements of type of the first vector in a vector block

    \param elem_size2
    (In) Number of elements of type of the second vector in a vector block

    \param elem_size3
    (In) Number of elements of type of the third vector in a vector block	    

    \param my_offset
    (In) Defines the block at which the writing operation should start
    In CPP_ASCII_GWriter, all processes must write consecutive blocks.

    \return Number of records written.
  */
  template<typename T1, typename T2, typename T3> int Write(const std::string& name, const std::string& name1, const T1* vec1, const std::string& name2, const T2* vec2,const std::string& name3, const T3* vec3, long num_global_elems, long num_my_elems, long elem_size1, long elem_size2, long elem_size3, long my_offset);

  //! Indicates wether the file has been properly opened
  /*!
      \return true, only if the file is NOT opened
  */
  bool operator!();

  //! This function does nothing. It is only implemented for convenience
  /*
      \return true
  */
  bool Select(const std::string&);

 private:
  bool Seek(const std::string&);
  template<typename T> int Write(const T& t) { file << t << ' '; return 0; }
  void Flush();
  bool Append();
  std::ofstream file;
  std::string filename;
  MPI_Comm comm;
  int mpi_size;
  int mpi_rank;
  std::map<std::string, std::fstream::pos_type> marks;
};

//! A generic class to write scalars and (block) vectors sequentially or in parallel into a HDF5 file.

/*! This class provides a generic interface to simplify the output of scalars and (block) vectors in the HDF5 file format.
    If parallel HDF5 was enabled, the functions will take advantage or parallel I/O.
    Otherwise, if the filesystem is not capable of parallel services, the class can be compiled to support only sequential I/O, 
    based on MPI synchronization.
*/

class HDF5_GWriter
{
 public:
  //! HDF5_GWriter constructor.
  /*! Creates a HDF5_GWriter object without performing any I/O operation.

    \param comm
    (In) Optional MPI communicator containing information on the parallel environment.
	    
    \return Pointer to a HDF5_GWriter object.

  */ 
  HDF5_GWriter(MPI_Comm comm= MPI_COMM_WORLD);

  //! HDF5_GWriter constructor.
  /*! Creates a HDF5_GWriter object and opens the file with the given filename for writing.
      If the file already exits, it will be overwritten.

    \param filename
    (In) Path name of the file to be written.
    The file will be overwritten, if already existent.

    \param comm
    (In) Optional MPI communicator containing information about the parallel environment.
	    
    \return Pointer to a HDF5_GWriter object.

  */ 
  HDF5_GWriter(const std::string& filename, MPI_Comm comm= MPI_COMM_WORLD);

  //! HDF5_GWriter destructor.
  /*! Closes any open files.
  */ 
  ~HDF5_GWriter();
  //! Opens the file with the given filename
  /*!
    \return true, only if successful
  */
  bool Open();

  //! Opens the file with the given filename
  /*!
    \param filename
    (In) Path name of the file to be written.
    The file will be overwritten, if already existent

    \return true, only if successful
  */
  bool Open(const std::string& filename);

  //! Closes any open files
  void Close();
  
  //! Generic function to write a scalar value
   /*!
    \param name
    (In) The name of the scalar HDF5 dataset.
 
    \param scalar
    (In) The scalar value to be written.
    Supported types, Integers, floating point numbers. 

    \return Number of records written.
  */
  template<typename T> int Write(const std::string& name, const T& scalar);
  
  //! Generic function to write a string value
   /*!
    \param name
    (In) the name of the scalar HDF5 dataset.
 
    \param scalar
    (In) The string value to be written.

    \return Number of records written.
  */
  int Write(const std::string& name, const std::string& scalar);

  //! Generic function to write a block vector

  /*! The block vector will be mapped into a dataspace of rank two.
      A block corresponds to a row while the block size equals the number of columns.

    \param name
    (In) The name of the HDF5 dataset.
 
    \param vec
    (In) Pointer to the block vector to be written.
    Supported types of vector elements, Integers, floating point numbers, strings.

    \param num_global_elems
    (In) Global number of vector blocks.

    \param num_my_elems
    (In) Local number of vector blocks.

    \param elem_size
    (In) Size of a vector block.
	    
    \param my_offset
    (In) Defines the block at which the writing operation should start.

    \return Number of records written.
  */
  template<typename T> int Write(const std::string& name, const T* vec, long num_global_elems, long num_my_elems, long elem_size, long my_offset);

  //! Generic function to write a block multityped vector.

  /*! A block will correspond to a row in a compound dataset. Each block consists of two subblocks of different types.
      The size and content of the subblocks is determined by the two given block vectors, i.e. each row is a
      concatenation of a block of vector one and a block of vector two.

    \param name
    (In) The name of the scalar HDF5 dataset.

    \param name1
    (In) Description of the first vector.
    Corresponds to the name of the first component of the compound dataset.
	    
    \param vec1
    (In) Pointer to the block vector of the first type to be written.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    Corresponds to the name of the second component of the compound dataset.

    \param vec2
    (In) Pointer to the block vector of the second type to be written.
    Supported types of vector elements are Integers, floating point numbers, strings.	    
    
    \param num_global_elems
    (In) Global number of vector blocks.

    \param num_my_elems
    (In) Local number of vector blocks.

    \param elem_size1
    (In) Number of elements of type of the first vector in a vector block.

    \param elem_size2
    (In) Number of elements of type of the second vector in a vector block.
	    
    \param my_offset
    (In) Defines the block at which the writing operation should start.

    \return Number of records written.
  */
  template<typename T1, typename T2> int Write(const std::string& name, const std::string& name1, const T1* vec1, const std::string& name2, const T2*vec2, long num_global_elems, long num_my_elems, long elem_size1, long elem_size2, long my_offset);

  //! Generic function to write a block multityped vector
  /*! A block will correspond to a row of a compound dataset. Each block consists of three subblocks of different types.
      The size and content of the subblocks is determined by the three given block vectors, i.e. each row is a
      concatenation of a block of vector one and a block of vector two and a block of vector three.

    \param name
    (In) The name of a section in the file.
	    

    \param name1
    (In) Description of the first vector.
    Corresponds to the name of the first component of the compound dataset.

    \param vec1
    (In) Pointer to the block vector of the first type to be written
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    Corresponds to the name of the second component of the compound dataset

    \param vec2
    (In) Pointer to the block vector of the second type to be written.
    Supported types of vector elements are Integers, floating point numbers, strings.	    
    
    \param name3
    (In) Description of the third vector.
    Corresponds to the name of the third component of the compound dataset
	    
    \param vec3
    (In) Pointer to the block vector of the third type to be written.
    Supported types of vector elements are Integers, floating point numbers, strings.

    \param num_global_elems
    (In) Global number of vector blocks.

    \param num_my_elems
    (In) Local number of vector blocks.

    \param elem_size1
    (In) Number of elements of type of the first vector in a vector block.

    \param elem_size2
    (In) Number of elements of type of the second vector in a vector block.

    \param elem_size3
    (In) Number of elements of type of the third vector in a vector block.	    

    \param my_offset
    (In) Defines the block at which the writing operation should start.

    \return Number of records written.
  */
  template<typename T1, typename T2, typename T3> int Write(const std::string& name, const std::string& name1, const T1* vec1, const std::string& name2, const T2* vec2, const std::string& name3, const T3* vec3, long num_global_elems, long num_my_elems, long elem_size1, long elem_size2, long elem_size3, long my_offset);
  
  //! Indicates wether the file has been properly opened
  /*!
      \return true, only if the file is NOT opened
  */
  bool operator!();

  //! Selects a group within a HDF5 file
  /*! If the group does not exist, it will be created
      
      \param name
      (In) Group name

      \return true, only if successful
  */ 
  bool Select(const std::string& name);

private:
  hid_t getNativeType(int, hsize_t=1);
  hid_t getNativeType(short, hsize_t=1);
  hid_t getNativeType(long, hsize_t=1);
  hid_t getNativeType(double, hsize_t=1);
  hid_t getNativeType(float, hsize_t=1);
  hid_t getNativeType(const std::string&, hsize_t=1);
  int Write(const std::string& name, hid_t type, const void* data, long num_global_elems, long num_my_elems, long elem_size, long my_offset);
  int Write(const std::string& name, hid_t type, const void* data);
  void copy(void* dest, const void* src, long nblocks, long block_size, long stride);
  std::string filename;
  std::string groupname;
  MPI_Comm comm;
  MPI_Info file_info;
  hid_t file;
  hid_t group;
  hid_t plist;
  int mpi_size;
  int mpi_rank;
  hsize_t offset[2];
  hsize_t stride[2];
  hsize_t index;
};


template<typename T> 
int CPP_ASCII_GWriter::Write(const std::string& name, const T& t) {
  int res = 0;
#ifdef HAVE_MPI
  //synchroniztion
  MPI_Barrier(comm);
  if (mpi_rank == 0) {
#endif
    Append();
    Seek(name); 
    res = Write(t); 
    file << '\n';
    Flush();
    Close();
#ifdef HAVE_MPI
  }
#endif
  return res; 
}

template<typename T>
int CPP_ASCII_GWriter::Write(const std::string& name, const T* data, long num_global_elems, long num_my_elems, long elem_size, long my_offset)
{
  long i=0;
#ifdef HAVE_MPI
  MPI_Status status;
  //synchroniztion
  MPI_Barrier(comm);
  if (mpi_rank != 0)
    MPI_Recv(&i, 1, MPI_INT, mpi_rank-1, 0, comm, &status);
#endif
  Append();
  Seek(name);
  long k=0;
  for (; i< num_my_elems + my_offset; ++i) {
    file << ' ';
    for (long j=0; j< elem_size; ++j)
      Write(data[k++]);
    file << '\n';
  }
  Close();
#ifdef HAVE_MPI
  //send my_offset+num_my_elems to next CPU
  if (mpi_rank != mpi_size-1)
    MPI_Send(&i, 1, MPI_INT, mpi_rank+1, 0, comm);
#endif
  return i;
}


template<typename T1, typename T2>
int CPP_ASCII_GWriter::Write(const std::string& name, 
			 const std::string& first_name, const T1* data1, 
			 const std::string& second_name, const T2* data2, 
			 long num_global_elems,
			 long num_my_elems,
			 long elem_size1, 
			 long elem_size2, 
			 long my_offset)
{
  //Wait for message
  long i=0;
#ifdef HAVE_MPI
  MPI_Status status;
  //synchroniztion
  MPI_Barrier(comm);
  if (mpi_rank != 0) 
    MPI_Recv(&i, 1, MPI_INT, mpi_rank-1, 0, comm, &status);
#endif
  Append();
  Seek(name);
  long k=0;
  long l=0;
  for (; i< num_my_elems+my_offset; ++i) {
    file << ' ';
    for (long j=0; j< elem_size1; ++j)
      Write(data1[k++]);
    for (long j=0; j< elem_size2; ++j)
      Write(data2[l++]);
    std::cout << '\n';
  }
  Close();
#ifdef HAVE_MPI
  //send my_offset+num_my_elems to next CPU
  if (mpi_rank != mpi_size-1)
    MPI_Send(&i, 1, MPI_INT, mpi_rank+1, 0, comm);
#endif
  return i;
}


template<typename T1, typename T2, typename T3>
int CPP_ASCII_GWriter::Write(const std::string& name, 
			 const std::string& first_name, const T1* data1,
			 const std::string& second_name, const T2* data2, 
			 const std::string& third_name, const T3* data3, 
			 long num_global_elems, 
			 long num_my_elems, 
			 long elem_size1,
			 long elem_size2,
			 long elem_size3,
			 long my_offset)
{
  //Wait for message
  long i=0;
#ifdef HAVE_MPI
  MPI_Status status;
  //synchroniztion
  MPI_Barrier(comm);
  if (mpi_rank != 0)
    MPI_Recv(&i, 1, MPI_INT, mpi_rank-1, 0, comm, &status);
#endif
  Append();
  Seek(name);
  long k=0;
  long l=0;
  long m=0;
  for (; i< num_my_elems+my_offset; ++i) {
    file << ' ';
    for (long j=0; j< elem_size1; ++j)
      Write(data1[k++]);
    for (long j=0; j< elem_size2; ++j)
      Write(data2[l++]);
    for (long j=0; j< elem_size3; ++j)
      Write(data3[m++]);
    std::cout << '\n';
  }
  Close();
#ifdef HAVE_MPI
  //send my_offset+num_my_elems to next CPU
  if (mpi_rank != mpi_size-1)
    MPI_Send(&i, 1, MPI_INT, mpi_rank+1, 0, comm);
#endif 
  return i;
}

template<typename T> 
int C_ASCII_GWriter::Write(const std::string& name, const T& t) { 
  long res = 0;
#ifdef HAVE_MPI
  //synchroniztion
  MPI_Barrier(comm);
  if (mpi_rank == 0) {
#endif
    Append();
    Seek(name); 
    res = Write(t); 
    fprintf(file, "\n");
    Flush();
    Close();
#ifdef HAVE_MPI
  }
#endif
  return res; 
}

template<typename T>
int C_ASCII_GWriter::Write(const std::string& name, const T* data, long num_global_elems, long num_my_elems, long elem_size, long my_offset)
{
  long i=0;
#ifdef HAVE_MPI
  MPI_Status status;
  //synchroniztion
  MPI_Barrier(comm);
  if (mpi_rank != 0) 
    MPI_Recv(&i, 1, MPI_INT, mpi_rank-1, 0, comm, &status);
#endif
  Append();
  Seek(name);
  long k=0;
  for (; i< num_my_elems+my_offset; ++i) {
    fprintf(file, " ");
    for (long j=0; j< elem_size; ++j)
      Write(data[k++]);
    fprintf(file, "\n");
  }
  Close();
#ifdef HAVE_MPI
  //send my_offset+num_my_elems to next CPU
  if (mpi_rank != mpi_size-1)
    MPI_Send(&i, 1, MPI_INT, mpi_rank+1, 0, comm);
#endif
  return i;
}


template<typename T1, typename T2>
int C_ASCII_GWriter::Write(const std::string& name, 
			 const std::string& first_name, const T1* data1, 
			 const std::string& second_name, const T2* data2, 
			 long num_global_elems,
			 long num_my_elems,
			 long elem_size1, 
			 long elem_size2, 
			 long my_offset)
{
  //Wait for message
  long i=0;
#ifdef HAVE_MPI
  MPI_Status status;
   //synchroniztion
  MPI_Barrier(comm);
  if (mpi_rank != 0) 
    MPI_Recv(&i, 1, MPI_INT, mpi_rank-1, 0, comm, &status);
#endif
  Append();
  Seek(name);
  long k=0;
  long l=0;
  for (; i< num_my_elems+my_offset; ++i) {
    fprintf(file, " ");
    for (long j=0; j< elem_size1; ++j)
      Write(data1[k++]);
    for (long j=0; j< elem_size2; ++j)
      Write(data2[l++]);
    fprintf(file, "\n");
  }
  Close();
#ifdef HAVE_MPI
  //send my_offset+num_my_elems to next CPU
  if (mpi_rank != mpi_size-1)
    MPI_Send(&i, 1, MPI_INT, mpi_rank+1, 0, comm);
#endif
  return i;
}


template<typename T1, typename T2, typename T3>
int C_ASCII_GWriter::Write(const std::string& name, 
			 const std::string& first_name, const T1* data1,
			 const std::string& second_name, const T2* data2, 
			 const std::string& third_name, const T3* data3, 
			 long num_global_elems, 
			 long num_my_elems, 
			 long elem_size1,
			 long elem_size2,
			 long elem_size3,
			 long my_offset)
{
  //Wait for message
  long i=0;
#ifdef HAVE_MPI
  MPI_Status status;
  //synchroniztion
  MPI_Barrier(comm);
  if (mpi_rank != 0) 
    MPI_Recv(&i, 1, MPI_INT, mpi_rank-1, 0, comm, &status);
#endif
  Append();
  Seek(name);
  long k=0;
  long l=0;
  long m=0;
  for (; i< num_my_elems+my_offset; ++i) {
    fprintf(file, " ");
    for (long j=0; j< elem_size1; ++j)
      Write(data1[k++]);
    for (long j=0; j< elem_size2; ++j)
      Write(data2[l++]);
    for (long j=0; j< elem_size3; ++j)
      Write(data3[m++]);
    fprintf(file, "\n");
  }
  Close();
#ifdef HAVE_MPI
  //send my_offset+num_my_elems to next CPU
  if (mpi_rank != mpi_size-1)
    MPI_Send(&i, 1, MPI_INT, mpi_rank+1, 0, comm);
#endif
  return i;
}


template<typename T>
int HDF5_GWriter::Write(const std::string& name, const T& data)
{
  return Write(name, getNativeType(data), &data);
}


template<typename T>
int HDF5_GWriter::Write(const std::string& name, const T* data, long num_global_elems, long num_my_elems, long elem_size, long my_offset)
{
  T dummy = 0;
  return Write(name, getNativeType(dummy), data,  num_global_elems, num_my_elems, elem_size, my_offset);
}


template<typename T1, typename T2>
int HDF5_GWriter::Write(const std::string& name, 
			 const std::string& first_name, const T1* data1, 
			 const std::string& second_name, const T2* data2, 
			 long num_global_elems,
			 long num_my_elems,
			 long elem_size1, 
			 long elem_size2, 
			 long my_offset)
{
  //Create datatype for first array
  T1 dummy1 = 0;
  hid_t t1 = getNativeType(dummy1, elem_size1);
  
  //Create datatype for second array
  T2 dummy2 = 0;
  hid_t t2 = getNativeType(dummy2, elem_size2);
  
  size_t sz = H5Tget_size(t1) + H5Tget_size(t2);
  hid_t s_tid = H5Tcreate (H5T_COMPOUND, sz);

  /* insert the component types at the appropriate offsets */
  
  H5Tinsert(s_tid, first_name.c_str(), 0, t1);
  H5Tinsert(s_tid, second_name.c_str(), H5Tget_size(t1), t2);

  //Create desired memory layout
  long nbytes1 = elem_size1*sizeof(T1);
  long nbytes2 = elem_size2*sizeof(T2);
  long stride = nbytes1+nbytes2;
  void* data = malloc(num_my_elems*(stride));
  copy(data, data1, num_my_elems, nbytes1, stride);
  copy((void*)((unsigned long)(data)+nbytes1), data2, num_my_elems, nbytes2, stride);
      
  long res = Write(name, s_tid, data, num_global_elems, num_my_elems, 1, my_offset);
  free(data);
  return res;
}


template<typename T1, typename T2, typename T3>
int HDF5_GWriter::Write(const std::string& name, 
			 const std::string& first_name, const T1* data1,
			 const std::string& second_name, const T2* data2, 
			 const std::string& third_name, const T3* data3, 
			 long num_global_elems, 
			 long num_my_elems, 
			 long elem_size1,
			 long elem_size2,
			 long elem_size3,
			 long my_offset)
{
  //Create datatype for first array
  T1 dummy1 = 0;
  hid_t t1 = getNativeType(dummy1, elem_size1);
  
  //Create datatype for second array
  T2 dummy2 = 0;
  hid_t t2 = getNativeType(dummy2, elem_size2);
  
   //Create datatype for third array
  T3 dummy3 = 0;
  hid_t t3 = getNativeType(dummy3, elem_size3);
  
  size_t sz = H5Tget_size(t1) + H5Tget_size(t2) + H5Tget_size(t3) ;
  hid_t s_tid = H5Tcreate (H5T_COMPOUND, sz);

  /* insert the component types at the appropriate offsets */
  
  H5Tinsert(s_tid, first_name.c_str(), 0, t1);
  H5Tinsert(s_tid, second_name.c_str(), H5Tget_size(t1), t2);
  H5Tinsert(s_tid, third_name.c_str(), H5Tget_size(t1) + H5Tget_size(t2), t3);

  //Create desired memory layout
  long nbytes1 = elem_size1*sizeof(T1);
  long nbytes2 = elem_size2*sizeof(T2);
  long nbytes3 = elem_size3*sizeof(T3);
  long stride = nbytes1+nbytes2+nbytes3;
  void* data = malloc(num_my_elems*(stride));
  copy(data, data1, num_my_elems, nbytes1,  stride);
  copy((void*)((unsigned long)(data)+nbytes1), data2, num_my_elems, nbytes2, stride);
  copy((void*)((unsigned long)(data)+nbytes1+nbytes2), data3, num_my_elems, nbytes3, stride);

  long res = Write(name, s_tid, data, num_global_elems, num_my_elems, 1, my_offset);
  free(data);
  return res;
}

#endif
