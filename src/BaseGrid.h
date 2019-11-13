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

#ifndef BASEGRID_H
#define BASEGRID_H

#include <map>
#include <vector>

#include "Config.h"
#include "ImageReader.h"

//! An interface to access implemented meshes.

/*! The BaseGrid is the abstraction of FE-meshes. This class should be inherited
 * protected, because all member are public accessable. The reason is, that a
 * reader cann fill in a image of a bone.
 */

class BaseGrid {
public:
  //! Constructor
  BaseGrid() : _grid(0){};

  //! Destructor
  virtual ~BaseGrid() {
    if (_grid)
      delete[] _grid;
  };

  //! Holds the Image in double format
  // Allocated by the reader
  double *_grid;

  //! Image dimension
  union {
    struct {
      t_coord gdim_x, gdim_y, gdim_z;
    };
    t_coord gdim[3];
  };

  //! Local image dimension
  union {
    struct {
      t_coord ldim_x, ldim_y, ldim_z;
    };
    t_coord ldim[3];
  };

  //! Cornercoordinate
  union {
    struct {
      t_coord corner_x, corner_y, corner_z;
    };
    t_coord corner[3];
  };

  //! Resolution in all directions. It is usually the same value for all
  //! dimension
  union {
    struct {
      double res_x, res_y, res_z;
    };
    double res[3];
  };

  //! Global Poison's ratio for the mesh
  double poisons_ratio;

  std::map<short, float> EmodMap;
  std::map<short, float> NuMap;
  std::map<float, short> invEmodMap;

  //! Index of the fix nodes. A node has 4 values (x,y,z,d) (d: direction)
  std::vector<unsigned short> fixed_nodes_coordinates;
  //! contains the displacements
  std::vector<float> fixed_nodes_values;
  //! Index of the loaded nodes. A node has 4 values (x,y,z,d) (d: direction)
  std::vector<unsigned short> loaded_nodes_coordinates;
  //! contains the loads
  std::vector<float> loaded_nodes_values;

private:
};

#endif /* BASEGRID_H */
