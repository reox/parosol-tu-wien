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


#ifndef JACOBIS_H
#define JACOBIS_H

#include "Solver.h"

#include "StiffnessMatrix.h"

//! This class implements a Jacobi Smoother.
/*! It uses external storage for the temporary vectors.
 *  The smoother can also be uses as a solver.
 */

class JacobiSmoother: public Solver
{
  public:
	/**
	 * @brief This is the constructor of the standard usage.
	 *
	 * @param M System matrix.
	 * @param Step Number of iteration.
	 * @param w Damping factor.
	 * @param ZeroStart If this parameter is set in assumes that it starts with x=0.
	 * @param m1 Vector, that is used to hold temporary results.
	 */

    JacobiSmoother(StiffnessMatrix &M, int Step, double w, bool ZeroStart, Eigen::VectorXd &m1):
      _mat(M),
      _idia(M.Diagonal()),
      _step(Step),
      _w(w),
      _zeroStart(ZeroStart),
      r(m1)
    {
    }


	/**
	 * @brief This constructor copies the configuration of the old Jacobi solver. ZeroStart must be extra set.
	 *
	 * @param old From this solver it reads the configuration.
	 * @param ZeroStart If this parameter is set in assumes that it starts with x=0.
	 */
    JacobiSmoother(JacobiSmoother &old, bool ZeroStart):
      _mat(old._mat),
      _idia(old._idia),
      _step(old._step),
      _w(old._w),
      _zeroStart(ZeroStart),
      r(old.r)
    {
    }

    virtual ~JacobiSmoother()
    {
    }

    int Solve(Eigen::VectorXd &b, Eigen::VectorXd &x);

    virtual const std::string  Label () const
    {
      return "Jacobi Smoother";
      }

  protected:
    StiffnessMatrix &_mat;
    Eigen::VectorXd &_idia;
    int _step;
    double _w;
    bool _zeroStart;
    Eigen::VectorXd &r;
};

#endif /* JACOBIS_H */
