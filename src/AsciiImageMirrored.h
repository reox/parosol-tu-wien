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



#ifndef ASCIIIMAGEMIRROR_H
#define ASCIIIMAGEMIRROR_H

#include "ImageReader.h"
#include "CPULayout.h"

/** It reads in an AsciiImage and mirrored it on the CPU.
 *   
 *  This class is used to run weak scalability test. The fileformat is documented in AsciiImage.h
 */


class AsciiImageMirrored : public ImageReader
{

public:

  //! Constructor
  
	/** 
	 * @brief Constructor 
	 * 
	 * @param filename
	 * @param layout Layout of the Grid on the cpu
	 */
  AsciiImageMirrored(const char* filename, CPULayout &layout); 

  //! destructor
  ~AsciiImageMirrored();


  //! Fills the image in to a grid

  /*!
    \param Grid
    (In) Pointer to the grid in which the Image is allocated and stored 
  */
  virtual int Scan(BaseGrid* Grid);


private:
  const char* _file;
  int _dim;
  CPULayout &_layout;
  union{
    struct {int _procx, _procy, _procz; };
    int _proc[3];
  };
  int MyPID;
};
  
#endif /* ASCIIIMAGE_H */

