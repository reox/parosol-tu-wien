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

#ifndef CHEBYSHEV_H
#define CHEBYSHEV_H

#include "Solver.h"
#include "StiffnessMatrix.h"

//! This class implements an Chebyshev solver.

/*! The Chebyshev solver is a symmetric procedure to solve a linear system.
 *  This implementation estimation the biggest and the smallest eigenvalue
 *  with the Lanczos algorithm. If it is used as smoother the lower eigenvalue
 *  can be set with a ration from the biggest estimate.
 */

class Chebyshev : public Solver {
public:
  /**
   * @brief Constructor of the Chebyshev solver.
   * If you use it to apply the inverse of the matrix choose zeroStart to be
   * true (example as Preeconditioner)
   *
   * @param M Stiffnessmatrix
   * @param degree Degree of the solver. This corresponds to the number of
   * iteration.
   * @param ZeroStart If this parameter is set in assumes that it starts with
   * x=0.
   * @param m1 temporary Vector, that is needed for the computation.
   * @param m2 temporary Vector, that is needed for the computation.
   * @param maxit The maximum number of iteration to estimate the eigenvalues.
   * (Default 10)
   * @param lmax Estimation of the biggest eigenvalue. (Default 0 -> means not
   * used)
   * @param lmin Estimation of the smallest eigenvalue. (Default 0 -> means not
   * used)
   * @param ratio Ratio of lmax/lmin. It is used to compute lmin.
   */
  Chebyshev(StiffnessMatrix &M, int degree, bool ZeroStart, Eigen::VectorXd &m1,
            Eigen::VectorXd &m2, int maxit = 10, double lmax = 0,
            double lmin = 0, double ratio = 16);
  /**
   * @brief This constructor copies the configuration of the old chebyshev
   * solver. ZeroStart must be extra set.
   *
   * @param old From this solver it reads the configuration.
   * @param ZeroStart If this parameter is set in assumes that it starts with
   * x=0.
   */
  Chebyshev(Chebyshev &old, bool ZeroStart);

  virtual ~Chebyshev() {}

  /**
   * @brief Solves the linear system \f[Ax=b\f] iteratively with the initial
   * guess \f[x\f]
   * @param[in] b RHS
   * @param[in,out] x initial guess and solution
   */
  virtual int Solve(Eigen::VectorXd &b, Eigen::VectorXd &x);

  virtual const std::string Label() const {
    return "Chebyshev Preconditioner/Smoother";
  }

protected:
  StiffnessMatrix &_mat;
  Eigen::VectorXd &_idia;
  long _ldofs;
  long _degree;
  double _lmax, _lmin;
  bool _zeroStart;

  Eigen::VectorXd &z, &p;
};

#endif /* CHEBYSHEV_H */
