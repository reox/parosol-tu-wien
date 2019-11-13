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

#ifndef EST_EV_H
#define EST_EV_H

#include "Solver.h"
#include "StiffnessMatrix.h"

struct eig {
  double large, small;
};

/** computes the eigenvalues of the preconditions system
 *
 * @param _mat the system matrix
 * @param _prec preconditioner
 * @param maxint number of iteration to compute the eigenvalues
 * @param m1 temporary Vector, that is needed for the computation.
 * @param m2 temporary Vector, that is needed for the computation.
 */

eig est_ev(StiffnessMatrix &_mat, Solver &_prec, int maxint,
           Eigen::VectorXd &m1, Eigen::VectorXd &m2);
#endif /* EST_EV_H */
