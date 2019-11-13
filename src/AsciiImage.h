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

#ifndef ASCIIIMAGE_H
#define ASCIIIMAGE_H

#include "CPULayout.h"
#include "ImageReader.h"

//! This call read an ascii image.

/*! AsciiImage reads an image from a ascii file.
 *
 *  The file must have the following format:<br>
 *  1. line x (int) y(int) z (int): x,y,z  dimension of the image <br>
 *  2. n (int) [E (float)]^n n indicates the number of E-moduli. E's E-moduli of
the materials. <br> *	3. i (int) x*y*z whitespace separated material
identifier for each voxel. <br>
 *
 *	All Values are just whitespace separated.
 *
 *  A minimal file:
\verbatim
2 2 2
2 0 1.1
1 1
1 1
0 1
1 1
\endverbatim
*/

class AsciiImage : public ImageReader {

public:
  //! Constructor

  /**
   * @brief Constructor
   *
   * @param filename
   * @param layout Layout of the Grid on the cpu
   */
  AsciiImage(const char *filename, CPULayout &layout);

  //! destructor
  ~AsciiImage();

  //! Fills the image in to a grid

  /*!
    \param Grid
    (In) Pointer to the grid in which the Image is allocated and stored
  */
  virtual int Scan(BaseGrid *Grid);

private:
  const char *_file;
  int _dim;
  CPULayout &_layout;
  union {
    struct {
      int _procx, _procy, _procz;
    };
    int _proc[3];
  };
  int MyPID;
};

#endif /* ASCIIIMAGE_H */
