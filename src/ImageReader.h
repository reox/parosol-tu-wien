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



#ifndef IMAGEREADER_H
#define IMAGEREADER_H

#include "BaseGrid.h"


//! An  abstract class for reading in different fileformat.

/*! This abstract class describe the interface to different fileformat.
*/
class ImageReader {
public:
  //! ImageReader Destructor
  virtual ~ImageReader(){};
  
  //! Print the mesh into a problem file

  /*!
    \param grid
    (In) Pointer to the grid in which the Image is allocated and stored 
  */
  virtual int Scan(BaseGrid* grid)=0;

private:
};


//some tools to read in
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

#endif /* IMAGEREADER_H */

