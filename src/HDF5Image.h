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

#ifndef HDF5IMAGE_H
#define HDF5IMAGE_H

#include "CPULayout.h"
#include "GReader.hpp"
#include "ImageReader.h"
#include <string>

//! A class to generate a bone image.

/*! HDF5Image read in a Image from a ascii file.
 */

class HDF5Image : public ImageReader {

public:
  //! Constructor

  /**
   * @brief Constructor
   *
   * @param filename
   * @param layout Layout of the Grid on the cpu
   */
  HDF5Image(std::string &filename, CPULayout &layout);

  ~HDF5Image();

  //! Fills the image in to a grid

  /*!
    \param Grid
    (In) Pointer to the grid in which the Image is allocated and stored
  */
  virtual int Scan(BaseGrid *Grid);

private:
  void ReadBC(HDF5_GReader &reader, std::string s,
              std::vector<unsigned short> &coordinates,
              std::vector<float> &values);
  std::string _file;
  int _dim;
  CPULayout &_layout;
  union {
    struct {
      int _procx, _procy, _procz;
    };
    int _proc[3];
  };
  int MyPID, mpi_size;
};

#endif /* HDF5IMAGE_H */
