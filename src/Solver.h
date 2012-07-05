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

#ifndef SOLVER_H
#define SOLVER_H

#include <string>
#include <eigen2/Eigen/Core>


//! An interface for solvers


class Solver 
{
public:

  virtual ~Solver() {};

  /** Solves \f[Ax=b\f]
   * @param [in] b rhs vector
   * @param [out] x solution vector
   */
  virtual int   Solve(Eigen::VectorXd &b, Eigen::VectorXd &x) = 0;

  /** Returns the name of the solver
   * 
   * @return return the name in an string
   */
  virtual const std::string  Label () const = 0;

};

#endif /* SOLVER_H */

