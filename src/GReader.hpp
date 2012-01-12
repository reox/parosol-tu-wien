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

//Generic Reader
#ifndef _GREADER_HPP_
#define _GREADER_HPP_

#define HAVE_MPI

#include <cstdlib> //GCC 4.3
#include<iostream>
#include<map>
#include<fstream>
#include<cstdio>
#include<hdf5.h>
#ifdef HAVE_MPI
#include <mpi.h>
#else
#define MPI_Comm void*
#define MPI_Info void*
#define MPI_COMM_WORLD 0
#define MPI_INFO_NULL 0
#endif

#ifndef H5_HAVE_PARALLEL
#ifdef HAVE_MPI
#define SEQUENTIAL_HDF5
#endif
#endif

//! A generic class to read scalars and (block) vectors from an ASCII file by using only C library functions.

/*! The class C_ASCII_GReader provides a generic methods to allow multiple processes to read vector/scalar data from the same file.
    Only C library I/O functions are invoked for the case that they perform best on the target architecture.
    To use this class in a parallel environment, an MPI implementation must be available.
*/

class C_ASCII_GReader
{
 public:
  //! C_ASCII_GReader constructor.
  /*! Creates a C_ASCII_GReader object without performing any I/O operation.

    \param comm
    (In) Optional MPI communicator containing information on the parallel environment.
	    
    \return Pointer to a C_ASCII_GReader object.

  */ 
  C_ASCII_GReader(MPI_Comm comm = MPI_COMM_WORLD);

  //! C_ASCII_GReader constructor.
  /*! Creates a C_ASCII_GReader object and opens the file with the given filename for reading.

    \param filename 
    (In) Path name of the file to be read.
    The file must already exist.  

    \param comm
    (In) Optional MPI communicator containing information on the parallel environment.
	    
    \return Pointer to a C_ASCII_GReader object.

  */ 
  C_ASCII_GReader(const std::string& filename, MPI_Comm comm = MPI_COMM_WORLD);

  //! C_ASCII_GReader destructor.
  /*! Closes any open files.
  */ 
  ~C_ASCII_GReader();
  
  //! Opens the file with the previously given filename.
   /*!
    \return true, only if successful.
  */
  bool Open();

  //! Opens the file with the given filename.
  /*!
    \param filename
    (In) Path name of the file to be read.
    The file must already exist.  

    \return true, only if successful.
  */
  bool Open(const std::string& filename);

  //! Closes any open files
  void Close();
  
  //! Generic function to read a scalar value.
  /*!
    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
 
    \param scalar
    (Out) The scalar value to be read.
    Supported types, Integers, floating point numbers, strings. 

    \return Number of records read.
  */
  template<typename T> int Read(const std::string& name, T& scalar) { Seek(name); return Read(scalar); }

   //! Generic function to read a block vector.

   /*! A block will correspond to a line in the ASCII file.
    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
 
    \param vec
    (Out) Pointer to the block vector to be read.
    Supported types of vector elements, Integers, floating point numbers, strings.

    \param num_global_elems
    (In) Global number of vector blocks.

    \param num_my_elems
    (In) Local number of vector blocks.

    \param elem_size
    (In) Size of a vector block.
	    
    \param my_offset
    (In) Defines the block at which the reading operation should start.
    In C_ASCII_GReader, all processes must read consecutive blocks.

    \return Number of records read.
  */
  template<typename T> int Read(const std::string& name, T* vec, long num_global_elems, long num_my_elems, long elem_size, long my_offset);

    //! Generic function to read an entire block vector.

    /*! A block will correspond to a line in the ASCII file.
    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.

    \param vec
    (Out) Pointer to the block vector to be read.
    Supported types of vector elements, Integers, floating point numbers, strings.

    \param num_elems
    (In) Global number of vector blocks.

    \param elem_size
    (In) Size of a vector block.

    \return Number of records read.
    */
    template<typename T> int Read(const std::string& name, T* vec, long num_elems, long elem_size) {
	return Read(name, vec, num_elems, num_elems, elem_size, 0);
    }

  //! Generic function to read a block multityped vector.

  /*! A block will correspond to a line in the ASCII file. Each block consists of two subblocks of different types.
      The size and content of the subblocks is determined by the two given block vectors, i.e. each line is a
      concatenation of a block of vector one and a block of vector two.

    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
    
    \param name1
    (In) Description of the first vector.
    This argument has no effect in C_ASCII_GReader.
	    
    \param vec1
    (Out) Pointer to the block vector of the first type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    This argument has no effect in C_ASCII_GReader.

    \param vec2
    (Out) Pointer to the block vector of the second type to be read.
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
    (In) Defines the block at which the reading operation should start.
    In C_ASCII_GReader, all processes must read consecutive blocks.

    \return number of records read.
  */
  template<typename T1, typename T2> int Read(const std::string& name, const std::string& name1, T1* vec1, const std::string& name2, T2* vec2, long num_global_elems, long num_my_elems, long elem_size1, long elem_size2, long my_offset);

  //! Generic function to read an entire block multityped vector.

  /*! A block will correspond to a line in the ASCII file. Each block consists of two subblocks of different types.
      The size and content of the subblocks is determined by the two given block vectors, i.e. each line is a
      concatenation of a block of vector one and a block of vector two.

    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
    
    \param name1
    (In) Description of the first vector.
    This argument has no effect in C_ASCII_GReader.
	    
    \param vec1
    (Out) Pointer to the block vector of the first type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    This argument has no effect in C_ASCII_GReader.

    \param vec2
    (Out) Pointer to the block vector of the second type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.	    
    
    \param num_elems
    (In) Global number of vector blocks.

    \param elem_size1
    (In) Number of elements of type of the first vector in a vector block.

    \param elem_size2
    (In) Number of elements of type of the second vector in a vector block.
	    
    \return number of records read.
  */
    template<typename T1, typename T2> int Read(const std::string& name, const std::string& name1, T1* vec1, const std::string& name2, T2* vec2, long num_elems, long elem_size1, long elem_size2) {
	return Read(name, name1, vec1, name2, vec2, num_elems, num_elems, elem_size1, elem_size2, 0);
    }
  
//! Generic function to read a block multityped vector.
  /*! A block will correspond to a line in the ASCII file. Each block consists of three subblocks of different types.
      The size and content of the subblocks is determined by the three given block vectors, i.e. each line is a.
      concatenation of a block of vector one and a block of vector two and a block of vector three.

    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
    
    \param name1
    (In) Description of the first vector.
    This argument has no effect in C_ASCII_GReader.
	    
    \param  vec1
    (Out) Pointer to the block vector of the first type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    This argument has no effect in C_ASCII_GReader.

    \param vec2
    (Out) Pointer to the block vector of the second type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.	    
    
    \param name3
    (In) Description of the third vector.
    This argument has no effect in C_ASCII_GReader.
	    
    \param vec3
    (Out) Pointer to the block vector of the third type to be read.
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
    (In) Defines the block at which the reading operation should start.
    In C_ASCII_GReader, all processes must read consecutive blocks.

    \return Number of records read.
  */
  template<typename T1, typename T2, typename T3> int Read(const std::string& name, const std::string& name1, T1* vec1, const std::string& name2, T2* vec2,const std::string& name3, T3* vec3, long num_global_elems, long num_my_elems, long elem_size1, long elem_size2, long elem_size3, long my_offset);


//! Generic function to read an entire block multityped vector.
  /*! A block will correspond to a line in the ASCII file. Each block consists of three subblocks of different types.
      The size and content of the subblocks is determined by the three given block vectors, i.e. each line is a.
      concatenation of a block of vector one and a block of vector two and a block of vector three.

    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
    
    \param name1
    (In) Description of the first vector.
    This argument has no effect in C_ASCII_GReader.
	    
    \param  vec1
    (Out) Pointer to the block vector of the first type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    This argument has no effect in C_ASCII_GReader.

    \param vec2
    (Out) Pointer to the block vector of the second type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.	    
    
    \param name3
    (In) Description of the third vector.
    This argument has no effect in C_ASCII_GReader.
	    
    \param vec3
    (Out) Pointer to the block vector of the third type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param num_elems
    (In) Global number of vector blocks.

    \param elem_size1 
    (In) Number of elements of type of the first vector in a vector block.

    \param elem_size2
    (In) Number of elements of type of the second vector in a vector block.

    \param elem_size3
    (In) Number of elements of type of the third vector in a vector block.	    

    \return Number of records read.
  */
    template<typename T1, typename T2, typename T3> int Read(const std::string& name, const std::string& name1, T1* vec1, const std::string& name2, T2* vec2, const std::string& name3, T3* vec3, long num_elems, long elem_size1, long elem_size2, long elem_size3){
	return Read(name, name1, vec1, name3, vec2, name3, vec3, num_elems, num_elems, elem_size1, elem_size2, elem_size3, 0);
    }
  
  //! Indicates wether the file has been properly opened
  /*!
      \return true, only if the file is NOT opened
  */
  bool operator!();

  //! This function does nothing. It is only implemented for convenience
  /*
      \return true
  */
  bool Select(const std::string& s);

 //! Skips lines beginning with the given character
  /*
      \param c
      (In) All consecutive lines beginning with c will be skipped

      \return Number of skipped lines
  */
  int Skip(char);

 private:
  bool Seek(const std::string& s);
  int Read(long& i);
  int Read(double& d);
  int Read(std::string& s);
  std::string filename;
  FILE* file;
  std::map<std::string, fpos_t> marks;
  MPI_Comm comm;
  int mpi_rank;
  int mpi_size;
};

//! A generic class to read scalars and (block) vectors from an ASCII file by using only C++ library functions.

/*! The class CPP_ASCII_GReader provides a generic methods to allow multiple processes to read vector/scalar data from the same file.
    Only C++ library I/O functions are invoked for the case that they perform best on the target architecture.
    To use this class in a parallel environment, an MPI implementation must be available.
*/
class CPP_ASCII_GReader
{
 public:
  //! CPP_ASCII_GReader constructor.
  /*! Creates a CPP_ASCII_GReader object without performing any I/O operation.

    \param comm
    (In) Optional MPI communicator containing information on the parallel environment.
	    
    \return Pointer to a CPP_ASCII_GReader object.

  */ 
  CPP_ASCII_GReader(MPI_Comm comm = MPI_COMM_WORLD);

  //! CPP_ASCII_GReader constructor.
  /*! Creates a CPP_ASCII_GReader object and opens the file with the given filename for reading.

    \param filename
    (In) Path name of the file to be read.
    The must already exist.

    \param comm
    (In) Optional MPI communicator containing information on the parallel environment.
	    
    \return Pointer to a CPP_ASCII_GReader object.

  */ 
  CPP_ASCII_GReader(const std::string& filename, MPI_Comm comm= MPI_COMM_WORLD);
  
  //! CPP_ASCII_GReader destructor.
  /*! Closes any open files.
  */ 
  ~CPP_ASCII_GReader();

  //! Opens the file with the previously given filename.
   /*!
    \return true, only if successful.
  */
  bool Open();

  //! Opens the file with the given filename.
  /*!
    \param filename
    (In) Path name of the file to be read.
    The file must already exist.

    \return true, only if successful.
  */
  bool Open(const std::string& filename);

  //! Closes any open files
  void Close();

  //! Generic function to read a scalar value.
   /*!
    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
    
    \param scalar
    (Out) The scalar value to be read.
    Supported types, Integers, floating point numbers, strings. 

    \return Number of records read.
  */
  template<typename T> int Read(const std::string& name, T& scalar) { Seek(name); return Read(scalar); }

  //! Generic function to read a block vector.

  /*! A block will correspond to a line in the ASCII file.
    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
 
    \param vec
    (Out) Pointer to the block vector to be read.
    Supported types of vector elements, Integers, floating point numbers, strings.

    \param num_global_elems
    (In) Global number of vector blocks.

    \param num_my_elems
    (In) Local number of vector blocks.

    \param elem_size
    (In) Size of a vector block.
	    
    \param my_offset
    (In) Defines the block at which the reading operation should start.
    In CPP_ASCII_GReader, all processes must read consecutive blocks.

    \return Number of records read.
  */
  template<typename T> int Read(const std::string& name, T* vec, long num_global_elems, long num_my_elems, long elem_size, long my_offset);

    //! Generic function to read an entire block vector.

    /*! A block will correspond to a line in the ASCII file.
    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.

    \param vec
    (Out) Pointer to the block vector to be read.
    Supported types of vector elements, Integers, floating point numbers, strings.

    \param num_elems
    (In) Global number of vector blocks.

    \param elem_size
    (In) Size of a vector block.

    \return Number of records read.
    */
    template<typename T> int Read(const std::string& name, T* vec, long num_elems, long elem_size) {
	return Read(name, vec, num_elems, num_elems, elem_size, 0);
    }

  //! Generic function to read a block multityped vector.

  /*! A block will correspond to a line in the ASCII file. Each block consists of two subblocks of different types.
      The size and content of the subblocks is determined by the two given block vectors, i.e. each line is a
      concatenation of a block of vector one and a block of vector two.

    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
    
    \param name1
    (In) Description of the first vector.
    This argument has no effect in CPP_ASCII_GReader.
	    
    \param vec1
    (Out) Pointer to the block vector of the first type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    This argument has no effect in CPP_ASCII_GReader.

    \param vec2
    (Out) Pointer to the block vector of the second type to be read.
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
    (In) Defines the block at which the reading operation should start.
    In CPP_ASCII_GReader, all processes must read consecutive blocks.

    \return Number of records read.
  */
  template<typename T1, typename T2> int Read(const std::string& name, const std::string& name1, T1* vec1, const std::string& name2, T2* vec2, long num_global_elems, long num_my_elems, long elem_size1, long elem_size2, long my_offset);

  //! Generic function to read an entire block multityped vector.

  /*! A block will correspond to a line in the ASCII file. Each block consists of two subblocks of different types.
      The size and content of the subblocks is determined by the two given block vectors, i.e. each line is a
      concatenation of a block of vector one and a block of vector two.

    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
    
    \param name1
    (In) Description of the first vector.
    This argument has no effect in CPP_ASCII_GReader.
	    
    \param vec1
    (Out) Pointer to the block vector of the first type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    This argument has no effect in CPP_ASCII_GReader.

    \param vec2
    (Out) Pointer to the block vector of the second type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.	    
    
    \param num_elems
    (In) Global number of vector blocks.

    \param elem_size1
    (In) Number of elements of type of the first vector in a vector block.

    \param elem_size2
    (In) Number of elements of type of the second vector in a vector block.
	    
    \return number of records read.
  */
    template<typename T1, typename T2> int Read(const std::string& name, const std::string& name1, T1* vec1, const std::string& name2, T2* vec2, long num_elems, long elem_size1, long elem_size2) {
	return Read(name, name1, vec1, name2, vec2, num_elems, num_elems, elem_size1, elem_size2, 0);
    }

  //! Generic function to read a block multityped vector.
  /*! A block will correspond to a line in the ASCII file. Each block consists of three subblocks of different types.
      The size and content of the subblocks is determined by the three given block vectors, i.e. each line is a.
      concatenation of a block of vector one and a block of vector two and a block of vector three.

    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
    
    \param name1
    (In) Description of the first vector.
    This argument has no effect in CPP_ASCII_GReader.
	    
    \param vec1
    (Out) Pointer to the block vector of the first type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    This argument has no effect in CPP_ASCII_GReader.

    \param vec2
    (Out) Pointer to the block vector of the second type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.	    
    
    \param name3
    (In) Description of the third vector.
    This argument has no effect in CPP_ASCII_GReader.
	    
    \param vec3
    (Out) Pointer to the block vector of the third type to be read.
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
    (In) Defines the block at which the reading operation should start.
    In CPP_ASCII_GReader, all processes must read consecutive blocks.

    \return Number of records read.
  */
  template<typename T1, typename T2, typename T3> int Read(const std::string& name, const std::string& name1, T1* vec1, const std::string& name2, T2* vec2,const std::string& name3, T3* vec3, long num_global_elems, long num_my_elems, long elem_size1, long elem_size2, long elem_size3, long my_offset);

//! Generic function to read an entire block multityped vector.
  /*! A block will correspond to a line in the ASCII file. Each block consists of three subblocks of different types.
      The size and content of the subblocks is determined by the three given block vectors, i.e. each line is a.
      concatenation of a block of vector one and a block of vector two and a block of vector three.

    \param name
    (In) The name of a section in the file.
    Each section in the file is associated with a file position.
    If there exists already a section with this name, the function will jump to the corresponding file position.
    
    \param name1
    (In) Description of the first vector.
    This argument has no effect in CPP_ASCII_GReader.
	    
    \param  vec1
    (Out) Pointer to the block vector of the first type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    This argument has no effect in CPP_ASCII_GReader.

    \param vec2
    (Out) Pointer to the block vector of the second type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.	    
    
    \param name3
    (In) Description of the third vector.
    This argument has no effect in CPP_ASCII_GReader.
	    
    \param vec3
    (Out) Pointer to the block vector of the third type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param num_elems
    (In) Global number of vector blocks.

    \param elem_size1 
    (In) Number of elements of type of the first vector in a vector block.

    \param elem_size2
    (In) Number of elements of type of the second vector in a vector block.

    \param elem_size3
    (In) Number of elements of type of the third vector in a vector block.	    

    \return Number of records read.
  */
    template<typename T1, typename T2, typename T3> int Read(const std::string& name, const std::string& name1, T1* vec1, const std::string& name2, T2* vec2, const std::string& name3, T3* vec3, long num_elems, long elem_size1, long elem_size2, long elem_size3){
	return Read(name, name1, vec1, name3, vec2, name3, vec3, num_elems, num_elems, elem_size1, elem_size2, elem_size3, 0);
    }

  //! Indicates wether the file has been properly opened
  /*!
      \return true, only if the file is NOT opened
  */
  bool operator!();

  //! This function does nothing. It is only implemented for convenience
  /*
      \return true
  */
  bool Select(const std::string& s);

  //! Skips lines beginning with the given character
  /*
      \param c
      (In) All consecutive lines beginning with c will be skipped

      \return Number of skipped lines
  */
  int Skip(char c);
  
 private:
  bool Seek(const std::string& s);
  template<typename T> int Read(T& t) { file >> t; return file.gcount();}
  std::ifstream file;
  std::string filename;
  std::map<std::string, std::fstream::pos_type> marks;
  MPI_Comm comm;
  int mpi_rank;
  int mpi_size;
};

//! A generic class to read scalars and (block) vectors from a HDF5 file.

/*! This class provides a generic interface to simplify the input of scalars and (block) vectors in the HDF5 file format.
    If parallel HDF5 was enabled, the functions will take advantage or parallel I/O.
*/

class HDF5_GReader
{
 public:
  //! HDF5_GReader constructor.
  /*! Creates a HDF5_GReader object without performing any I/O operation.

    \param comm
    (In) Optional MPI communicator containing information on the parallel environment.
	    
    \return Pointer to a HDF5_GReader object.

  */ 
  HDF5_GReader(MPI_Comm comm= MPI_COMM_WORLD);

  //! HDF5_GReader constructor.
  /*! Creates a HDF5_GReader object and opens the file with the given filename for reading.

    \param filename
    (In) Path name of the file to be read.
    The file must already exist.

    \param comm
    (In) Optional MPI communicator containing information on the parallel environment.
	    
    \return Pointer to a HDF5_GReader object.

  */ 
  HDF5_GReader(const std::string& filename, MPI_Comm comm= MPI_COMM_WORLD);

  //! HDF5_GReader destructor.
  /*! Closes any open files.
  */ 
  ~HDF5_GReader();

  //! Opens the file with the given filename
  /*!
    \return true, only if successful
  */
  bool Open();

  //! Opens the file with the given filename
  /*!
    \param filename
    (In) Path name of the file to be read.
    The file must already exist.

    \return true, only if successful
  */
  bool Open(const std::string& filename);

  //! Closes any open files
  void Close();

  //! Generic function to read a scalar value
   /*!
    \param name
    (In) The name of the scalar HDF5 dataset.
 
    \param scalar
    (Out) The scalar value to be read.
    Supported types, Integers, floating point numbers. 

    \return Number of records read.
  */
  template<typename T> int Read(const std::string& name , T& scalar);

    //! Generic function to read a string value
   /*!
    \param name
    (In) The name of the scalar HDF5 dataset.
 
    \param scalar
    (Out) The string value to be read.

    \return Number of records read.
  */
  int Read(const std::string& name, std::string& scalar);

  //! Generic function to partially read a block vector

  /*! The block vector will be mapped into a dataspace of rank two.
      A block corresponds to a row while the block size equals the number of columns.

    \param name
    (In) The name of the HDF5 dataset.
 
    \param vec
    (Out) Pointer to the block vector to be read.
    Supported types of vector elements, Integers, floating point numbers, strings.

    \param num_global_elems
    (In) Global number of vector blocks.

    \param num_my_elems
    (In) Local number of vector blocks.

    \param elem_size
    (In) Size of a vector block.
	    
    \param my_offset
    (In) Defines the block at which the reading operation should start.

    \return Number of records read.
  */
  template<typename T> int Read(const std::string& name, T* vec, long num_global_elems, long num_my_elems, long elem_size, long my_offset);

  //! Generic function to read an entire block vector and replicate it on every process.

  /*! The block vector will be mapped into a dataspace of rank two.
    A block corresponds to a row while the block size equals the number of columns.

  \param name
  (In) The name of the HDF5 dataset.

  \param vec
  (Out) Pointer to the block vector to be read.
  Supported types of vector elements, Integers, floating point numbers, strings.

  \param num_elems
  (In) Global number of vector blocks.

  \param elem_size
  (In) Size of a vector block.

  \return Number of records read.
    */
  template<typename T> int Read(const std::string& name, T* vec, long num_elems, long elem_size);


  //! Generic function to partially read a block multityped vector.

  /*! A block will correspond to a row in a compound dataset. Each block consists of two subblocks of different types.
      The size and content of the subblocks is determined by the two given block vectors, i.e. each row is a
      concatenation of a block of vector one and a block of vector two.

    \param name
    (In) The name of the scalar HDF5 dataset.

    \param name1
    (In) Description of the first vector.
    Corresponds to the name of the first component of the compound dataset.
	    
    \param vec1
    (Out) Pointer to the block vector of the first type to be read.
	    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    Corresponds to the name of the second component of the compound dataset.

    \param vec2
    (Out) Pointer to the block vector of the second type to be read.
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
    (In) Defines the block at which the reading operation should start.

    \return Number of records read.
  */
  template<typename T1, typename T2> int Read(const std::string& name, const std::string& name1, T1* vec1, const std::string& name2, T2* vec2, long num_global_elems, long num_my_elems, long elem_size1, long elem_size2, long my_offset);

  //! Generic function to read an entire block multityped vector an replicate it on every process.

  /*! A block will correspond to a row in a compound dataset. Each block consists of two subblocks of different types.
      The size and content of the subblocks is determined by the two given block vectors, i.e. each row is a
      concatenation of a block of vector one and a block of vector two.

  \param name
  (In) The name of the scalar HDF5 dataset.

  \param name1
  (In) Description of the first vector.
  Corresponds to the name of the first component of the compound dataset.

  \param vec1
  (Out) Pointer to the block vector of the first type to be read.
  Supported types of vector elements are Integers, floating point numbers, strings.

  \param name2
  (In) Description of the second vector.
  Corresponds to the name of the second component of the compound dataset.

  \param vec2
  (Out) Pointer to the block vector of the second type to be read.
  Supported types of vector elements are Integers, floating point numbers, strings.

  \param num_elems
  (In) Global number of vector blocks.

  \param elem_size1
  (In) Number of elements of type of the first vector in a vector block.

  \param elem_size2
  (In) Number of elements of type of the second vector in a vector block.

  \return Number of records read.
  */
  template<typename T1, typename T2> int Read(const std::string& name, const std::string& name1, T1* vec1, const std::string& name2, T2* vec2, long num_elems, long elem_size1, long elem_size2);

    //! Generic function to partially read a block multityped vector
  /*! A block will correspond to a row of a compound dataset. Each block consists of three subblocks of different types.
      The size and content of the subblocks is determined by the three given block vectors, i.e. each row is a
      concatenation of a block of vector one and a block of vector two and a block of vector three.

    \param name
    (In) The name of a section in the file.
	    

    \param name1
    (In) Description of the first vector.
    Corresponds to the name of the first component of the compound dataset.

    \param vec1
    (Out) Pointer to the block vector of the first type to be read
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    Corresponds to the name of the second component of the compound dataset

    \param vec2
    (Out) Pointer to the block vector of the second type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.	    
    
    \param name3
    (In) Description of the third vector.
    Corresponds to the name of the third component of the compound dataset
	    
    \param vec3
    (Out) Pointer to the block vector of the third type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.

    \param num_global_elems
    Global number of vector blocks.

    \param num_my_elems
    (In) Local number of vector blocks.

    \param elem_size1
    (In) Number of elements of type of the first vector in a vector block.

    \param elem_size2
    (In) Number of elements of type of the second vector in a vector block.

    \param elem_size3
    (In) Number of elements of type of the third vector in a vector block.	    

    \param my_offset
    (In) Defines the block at which the reading operation should start.

    \return Number of records read.
  */
  template<typename T1, typename T2, typename T3> int Read(const std::string& name, const std::string& name1, T1* vec1, const std::string& name2, T2* vec2,const std::string& name3, T3* vec3, long num_global_elems, long num_my_elems, long elem_size1, long elem_size2, long elem_size3, long my_offset);

    //! Generic function to read an entire block multityped vector and replicate it on every process.
  /*! A block will correspond to a row of a compound dataset. Each block consists of three subblocks of different types.
      The size and content of the subblocks is determined by the three given block vectors, i.e. each row is a
      concatenation of a block of vector one and a block of vector two and a block of vector three.

    \param name
    (In) The name of a section in the file.
	    

    \param name1
    (In) Description of the first vector.
    Corresponds to the name of the first component of the compound dataset.

    \param vec1
    (Out) Pointer to the block vector of the first type to be read
    Supported types of vector elements are Integers, floating point numbers, strings.
    
    \param name2
    (In) Description of the second vector.
    Corresponds to the name of the second component of the compound dataset

    \param vec2
    (Out) Pointer to the block vector of the second type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.	    
    
    \param name3
    (In) Description of the third vector.
    Corresponds to the name of the third component of the compound dataset
	    
    \param vec3
    (Out) Pointer to the block vector of the third type to be read.
    Supported types of vector elements are Integers, floating point numbers, strings.

    \param num_elems
    Global number of vector blocks.

    \param elem_size1
    (In) Number of elements of type of the first vector in a vector block.

    \param elem_size2
    (In) Number of elements of type of the second vector in a vector block.

    \param elem_size3
    (In) Number of elements of type of the third vector in a vector block.	    

    \return Number of records read.
  */
  template<typename T1, typename T2, typename T3> int Read(const std::string& name, const std::string& name1, T1* vec1, const std::string& name2, T2* vec2,const std::string& name3, T3* vec3, long num_elems, long elem_size1, long elem_size2, long elem_size3);
 
  
  //! Generic Readfunction for a hyperslap with the dimension count and a offset.
  /*!
    \param name
    (In) The name of a section in the file.
    
    \param field
    (Out) Array of of the dimension dims and the size of count
    
    \param my_offset
    (In) Starting offset of the hyperslap
    
    \param count
    (In) Size in each dimenson of the hyperslap

    \param number_dims
    Number of dimensions of the dataset
  */

  template<typename T> int Read(const std::string& name, T* field, hsize_t* my_offset, hsize_t *count, const int number_dims);

  //! Return the Size of the dataset
  /*!
    \param name
    (In) The name of a section in the file.
    
    \param size
    (Out) Array of the size of the dimension

    \param number_dims
    Number of dimensions of the dataset
  */

  int GetSizeOfDataset(const std::string& name, hsize_t *size, const int number_dims);
  
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
  bool Select(std::string name);

  //! This function does nothing. It is only implemented for convenience
  /*
      \return true
  */
  int Skip(char);
 
private:

  hid_t getNativeType(int, hsize_t=1);
  hid_t getNativeType(double, hsize_t=1);
  hid_t getNativeType(float, hsize_t=1);
  hid_t getNativeType(const std::string&, hsize_t=1);
  hid_t getNativeType(short, hsize_t t=1);
  int Read(const std::string& name, hid_t type, void* data, long num_global_elems, long num_my_elems, long elem_size, long my_offset);
  int Read(const std::string& name, hid_t type, void* data);
  int Read(const std::string& name, hid_t type, void* data, hsize_t* my_offset, hsize_t *count, const int dims);

  void copy(void* dest, const void* src, long nblocks, long block_size, long stride);
  std::string filename;
  hid_t file;
  hid_t group;
  hid_t plist;
  hsize_t offset[2];
  hsize_t stride[2];
  hsize_t index;
  MPI_Comm comm;
  MPI_Info file_info;
  int mpi_rank;
  int mpi_size;
};

template<typename T> 
int C_ASCII_GReader::Read(const std::string& name, T* data, long num_global_elems, long num_my_elems, long elem_size, long my_offset)
{
  Seek(name);
  T value;
  long k=0;
  long i;
  for (i=0; i< my_offset; ++i) 
    for (long j=0; j < elem_size; ++j)
      Read(value);
      
  for (; i< num_my_elems+my_offset; ++i) {
    for (long j=0; j < elem_size; ++j) {
      Read(value);
      data[k++] = value;
    }
  }

  for (; i< num_global_elems; ++i) 
    for (long j=0; j < elem_size; ++j)
      Read(value);

  return i;
}

template<typename T1, typename T2>
int C_ASCII_GReader::Read(const std::string& name, 
			 const std::string& first_name, T1* data1, 
			 const std::string& second_name, T2* data2, 
			 long num_global_elems,
			 long num_my_elems,
			 long elem_size1, 
			 long elem_size2, 
			 long my_offset)
{
  Seek(name);
  T1 value1;
  T2 value2;
  long k=0;
  long l=0;
  long i;
  for (i=0; i< my_offset; ++i) {
    for (long j=0; j < elem_size1; ++j)
      Read(value1);
    for (long j=0; j < elem_size2; ++j)
      Read(value2);
  }
  
  for (; i< num_my_elems+my_offset; ++i) {
    for (long j=0; j < elem_size1; ++j) {
      Read(value1);
      data1[k++] = value1;
    }
    for (long j=0; j < elem_size2; ++j) {
      Read(value2);
      data2[l++] = value2;
    }
  }

  for (; i< num_global_elems; ++i) {
    for (long j=0; j < elem_size1; ++j)
      Read(value1);
    for (long j=0; j < elem_size2; ++j)
      Read(value2);
  }
  return i;
}


template<typename T1, typename T2, typename T3>
int C_ASCII_GReader::Read(const std::string& name, 
			  const std::string& first_name, T1* data1,
			  const std::string& second_name, T2* data2, 
			  const std::string& third_name, T3* data3, 
			  long num_global_elems, 
			  long num_my_elems, 
			  long elem_size1,
			  long elem_size2,
			  long elem_size3,
			  long my_offset)
{
  Seek(name);
  T1 value1;
  T2 value2;
  T3 value3;
  long k=0;
  long l=0;
  long m=0;
  long i;
  for (i=0; i< my_offset; ++i) {
    for (long j=0; j < elem_size1; ++j)
      Read(value1);
    for (long j=0; j < elem_size2; ++j)
      Read(value2);
    for (long j=0; j < elem_size3; ++j)
      Read(value3);
  }
  
  for (; i< num_my_elems+my_offset; ++i) {
    for (long j=0; j < elem_size1; ++j) {
      Read(value1);
      data1[k++] = value1;
    }
    for (long j=0; j < elem_size2; ++j) {
      Read(value2);
      data2[l++] = value2;
    }
    for (long j=0; j < elem_size3; ++j) {
      Read(value3);
      data3[m++] = value3;
    }
  }

  for (; i< num_global_elems; ++i) {
    for (long j=0; j < elem_size1; ++j)
      Read(value1);
    for (long j=0; j < elem_size2; ++j)
      Read(value2);
    for (long j=0; j < elem_size3; ++j)
      Read(value3);
  }
  return i;
}


template<typename T> 
int CPP_ASCII_GReader::Read(const std::string& name, T* data, long num_global_elems, long num_my_elems, long elem_size, long my_offset)
{
  Seek(name);
  T value;
  long k=0;
  long i;
  for (i=0; i< my_offset; ++i)
    for (long j=0; j < elem_size; ++j)
      Read(value);

  for (; i< num_my_elems + my_offset; ++i) {
    for (long j=0; j < elem_size; ++j) {
      Read(value);
      data[k++] = value;
    }
  }

  for (; i< num_global_elems; ++i)
    for (long j=0; j < elem_size; ++j)
      Read(value);

  return i;
}


template<typename T1, typename T2>
int CPP_ASCII_GReader::Read(const std::string& name, 
			    const std::string& first_name, T1* data1, 
			    const std::string& second_name, T2* data2, 
			    long num_global_elems,
			    long num_my_elems,
			    long elem_size1, 
			    long elem_size2, 
			    long my_offset)
{
  Seek(name);
  T1 value1;
  T2 value2;
  long k=0;
  long l=0;
  long i;
  for (i=0; i< my_offset; ++i) {
    for (long j=0; j < elem_size1; ++j)
      Read(value1);
    for (long j=0; j < elem_size2; ++j)
      Read(value2);
  }
  
  for (; i< num_my_elems + my_offset; ++i) {
    for (long j=0; j < elem_size1; ++j) {
      Read(value1);
      data1[k++] = value1;
    }
    for (long j=0; j < elem_size2; ++j) {
      Read(value2);
      data2[l++] = value2;
    }
  }

  for (; i< num_global_elems; ++i) {
    for (long j=0; j < elem_size1; ++j)
      Read(value1);
    for (long j=0; j < elem_size2; ++j)
      Read(value2);
  }
  return i;
}


template<typename T1, typename T2, typename T3>
int CPP_ASCII_GReader::Read(const std::string& name, 
			    const std::string& first_name, T1* data1,
			    const std::string& second_name, T2* data2, 
			    const std::string& third_name, T3* data3, 
			    long num_global_elems, 
			    long num_my_elems, 
			    long elem_size1,
			    long elem_size2,
			    long elem_size3,
			    long my_offset)
{
  Seek(name);
  T1 value1;
  T2 value2;
  T3 value3;
  long k=0;
  long l=0;
  long m=0;
  long i;
  for (i=0; i< my_offset; ++i) {
    for (long j=0; j < elem_size1; ++j)
      Read(value1);
    for (long j=0; j < elem_size2; ++j)
      Read(value2);
    for (long j=0; j < elem_size3; ++j)
      Read(value3);
  }
  
  for (; i< num_my_elems+my_offset; ++i) {
    for (long j=0; j < elem_size1; ++j) {
      Read(value1);
      data1[k++] = value1;
    }
    for (long j=0; j < elem_size2; ++j) {
      Read(value2);
      data2[l++] = value2;
    }
    for (long j=0; j < elem_size3; ++j) {
      Read(value3);
      data3[m++] = value3;
    }
  }

  for (; i< num_global_elems; ++i) {
    for (long j=0; j < elem_size1; ++j)
      Read(value1);
    for (long j=0; j < elem_size2; ++j)
      Read(value2);
    for (long j=0; j < elem_size3; ++j)
      Read(value3);
  }
  return i;
}


template<typename T>
int HDF5_GReader::Read(const std::string& name, T* field, hsize_t* my_offset, hsize_t* count, const int dims)
{
  T dummy = 0;
  return Read(name, getNativeType(dummy), field, my_offset, count, dims);
  

}



template<typename T>
int HDF5_GReader::Read(const std::string& name, T& data)
{
  int res = 0;
  if (mpi_rank == 0)
    res = Read(name, getNativeType(data), &data);
#ifdef HAVE_MPI
  MPI_Bcast(&data, sizeof(T), MPI_BYTE, 0, comm);
#endif
  return res;
}


template<typename T>
int HDF5_GReader::Read(const std::string& name, T* data, long num_global_elems, long num_my_elems, long elem_size, long my_offset)
{
  T dummy = 0;
  return Read(name, getNativeType(dummy), data, num_global_elems, num_my_elems, elem_size, my_offset);
}

template<typename T>
int HDF5_GReader::Read(const std::string& name, T* data, long num_elems, long elem_size)
{
    int res;
    if (mpi_rank == 0)
	res = Read(name, data, num_elems, num_elems, elem_size, 0);
#ifdef HAVE_MPI
    MPI_Bcast(data, num_elems*elem_size*sizeof(T), MPI_BYTE, 0, comm);
#endif
    return res;
}

template<typename T1, typename T2>
int HDF5_GReader::Read(const std::string& name, 
			 const std::string& first_name, T1* data1, 
			 const std::string& second_name, T2* data2, 
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

  long nbytes1 = elem_size1*sizeof(T1);
  long nbytes2 = elem_size2*sizeof(T2);
  long stride = nbytes1+nbytes2;
  void* data = malloc(num_my_elems*(stride));      
  int res = Read(name, s_tid, data, num_global_elems, num_my_elems, 1, my_offset);

  //Create desired memory layout
  copy(data1, data, num_my_elems, nbytes1, stride);
  copy(data2, (void*)((unsigned long)(data)+nbytes1), num_my_elems, nbytes2, stride);

  free(data);
  return res;
}

template<typename T1, typename T2>
int HDF5_GReader::Read(const std::string& name,
                            const std::string& first_name, T1* data1,
                            const std::string& second_name, T2* data2,
                            long num_elems,
                            long elem_size1,
                            long elem_size2)
{
    long res;
    if (mpi_rank == 0)
	res = Read(name, first_name, data1, second_name, data2, num_elems, num_elems, elem_size1, elem_size2, 0);
#ifdef HAVE_MPI
    MPI_Bcast(data1, num_elems*elem_size1*sizeof(T1), MPI_BYTE, 0, comm);
    MPI_Bcast(data2, num_elems*elem_size2*sizeof(T2), MPI_BYTE, 0, comm);
#endif
    return res;
}


template<typename T1, typename T2, typename T3>
int HDF5_GReader::Read(const std::string& name, 
			 const std::string& first_name, T1* data1,
			 const std::string& second_name, T2* data2, 
			 const std::string& third_name, T3* data3, 
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

  int res = Read(name, s_tid, data, num_global_elems, num_my_elems, 1, my_offset);
  copy(data1, data, num_my_elems, nbytes1,  stride);
  copy(data2, (void*)((unsigned long)(data)+nbytes1), num_my_elems, nbytes2, stride);
  copy(data3, (void*)((unsigned long)(data)+nbytes1+nbytes2), num_my_elems, nbytes3, stride);

  free(data);
  return res;
}

template<typename T1, typename T2, typename T3>
int HDF5_GReader::Read(const std::string& name,
                            const std::string& first_name, T1* data1,
                            const std::string& second_name, T2* data2,
                            const std::string& third_name, T3* data3,
                            long num_elems,
                            long elem_size1,
                            long elem_size2,
                            long elem_size3)
{
    long res;
    if (mpi_rank == 0)
	res = Read(name, first_name, data1, second_name, data2, third_name, data3, num_elems, num_elems, elem_size1, elem_size2, elem_size3, 0);
#ifdef HAVE_MPI
    MPI_Bcast(data1, num_elems*elem_size1*sizeof(T1), MPI_BYTE, comm, 0);
    MPI_Bcast(data1, num_elems*elem_size2*sizeof(T2), MPI_BYTE, comm, 0);
    MPI_Bcast(data1, num_elems*elem_size3*sizeof(T3), MPI_BYTE, comm, 0);
#endif
    return res;
}

#endif
