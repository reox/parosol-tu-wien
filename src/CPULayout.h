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

#ifndef CPULAYOUT_H
#define CPULAYOUT_H

#include <algorithm>
#include <mpi.h>
#include <ostream>
#include <vector>

/** This class generate and holds geometric cpu layout.
 * It tries to create an cuboid cpu grid by factoring the number
 * of cpu. If this number is prime, then the grid is linear
 */

class CPULayout {
public:
  CPULayout() {
    MPI_Comm_rank(MPI_COMM_WORLD, &_MyPID);
    MPI_Comm_size(MPI_COMM_WORLD, &_PSize);
    std::vector<unsigned> dims;
    computedims(_PSize, dims);
    for (int i = 0; i < 3; i++)
      _proc[i] = dims[i];

    _cpucord_x = _MyPID % _proc_x;
    _cpucord_y = (_MyPID / _proc_x) % _proc_y;
    _cpucord_z = (_MyPID / _proc_x) / _proc_y;
  }

  ~CPULayout() {}

  //! Debug function
  void Print();

  //! Returns the coordinate in the grid of the CPU
  const long *const CPUCoord() const { return _cpucord; }
  //! Returns the dimension of the grid.
  const long *const CPUGrid() const { return _proc; }

  int PID() { return _MyPID; }
  int PSIZE() { return _PSize; }

  friend std::ostream &operator<<(std::ostream &stream,
                                  const CPULayout &layout);

protected:
  //! Processeor id
  int _MyPID;
  int _PSize;
  //! CPU Coordinat
  union {
    struct {
      long _cpucord_x, _cpucord_y, _cpucord_z;
    };
    long _cpucord[3];
  };

  //! Number of Processor in each direction
  union {
    struct {
      long _proc_x, _proc_y, _proc_z;
    };
    long _proc[3];
  };

  // some print function
  std::ostream &print(std::ostream &stream) const {
    stream << "Layout: " << std::endl;
    return stream;
  }

private:
  // Compute the factorization of nr_proc and store it in the vector. This
  // is needed to compute the grid
  void factor(unsigned nr_proc, std::vector<unsigned> &factors) {
    while (nr_proc % 2 == 0) {
      factors.push_back(2);
      nr_proc = nr_proc / 2;
    }

    unsigned p = 3;
    while (p * p <= nr_proc) {
      if (nr_proc % p == 0) {
        factors.push_back(p);
        nr_proc = nr_proc / p;
      } else {
        p += 2;
      }
    }
    if (nr_proc > 1)
      factors.push_back(nr_proc);
    return;
  }

  // Tries do compute equally sized dimension by distribute the factor of nr_cpu
  // equally
  void computedims(unsigned nr_cpu, std::vector<unsigned> &dims) {
    std::vector<unsigned> factors;
    factors.reserve(20);
    factor(nr_cpu, factors);
    // take the 3 biggest factors
    for (int i = 0; i < 3; i++) {
      if (!factors.empty()) {
        dims.push_back(factors.back());
        factors.pop_back();
      } else {
        dims.push_back(1);
      }
    }
    // add an factor to the smallest dim
    sort(dims.begin(), dims.end());
    while (!factors.empty()) {
      dims[0] *= factors.back();
      factors.pop_back();
      sort(dims.begin(), dims.end());
    }
    sort(dims.begin(), dims.end());
    return;
  }
};

#endif /* CPULayoutH */
