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




#include "est_ev.h"
#include "float.h"
#include "tri_eig.h"

#include <eigen2/Eigen/Core>
#include <eigen2/Eigen/Array>
#include "mpi.h"

/// Compute the EV of a tridiagonal symmetric matrix
int QL_Tridiagonal_Symmetric_Matrix( Eigen::VectorXd &diagonal, Eigen::VectorXd &p_off,
                                   int n, int max_iteration_count);

eig est_ev(StiffnessMatrix &_mat, Solver &_prec, int maxit, Eigen::VectorXd &r, Eigen::VectorXd &d)
{
  double alpha = 0, beta =0; // in first iteration beta_old is not used
  double alpha_old, res0, delta_new, delta_old;
  double _tol = 1e-7;
  Eigen::VectorXd dia_T(maxit+1);
  Eigen::VectorXd off_T(maxit+1);
  //_mat.SetVectorRandom(r);
  r.setConstant(_mat.GetNrDofs(), 1.0);

  int i = 0;
  _prec.Solve(r,d);
  Eigen::VectorXd s = r;
  res0 =  delta_new = _mat.dot(r,d);
  res0 = sqrt(res0);
  while ( i < maxit && sqrt(delta_new)/res0 > _tol) {
    _mat.Apply(d,s);
    alpha_old = alpha;
    alpha = delta_new /_mat.dot(d,s);

    //compute tridiagonal matrix
    dia_T[i] = 1/alpha;
    if (i > 0) {
      dia_T[i] += beta/alpha_old;
      off_T[i-1] = sqrt(beta)/alpha_old;
    }
    r = r - alpha*s;

    _prec.Solve(r,s);
    delta_old = delta_new;
    delta_new = _mat.dot(r,s);

    beta = delta_new/delta_old;
    d = s + beta*d;
    i++;
  }

  QL_Tridiagonal_Symmetric_Matrix(dia_T, off_T,i,30);

  double l_min=dia_T[0];
  double l_max=dia_T[0];
  for(int ii = 1; ii < i; ii++) {
    if (dia_T[ii] < l_min)
      l_min = dia_T[ii];
    if (dia_T[ii] > l_max)
      l_max = dia_T[ii];
  }


  eig ret;
  ret.large=l_max;
  ret.small=l_min;
  return ret;
}

