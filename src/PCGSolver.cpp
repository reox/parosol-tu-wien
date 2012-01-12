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


#include "PCGSolver.h"

PCGSolver::PCGSolver(StiffnessMatrix &M, Solver &precon, double tol, int maxIter, bool verbose, int freq):
  _mat(M),
  _prec(precon),
  _maxIter(maxIter),
  _tol(tol),
  _verbose(verbose),
  _freq(freq)
{
}

PCGSolver::~PCGSolver()
{
}

int PCGSolver::Solve(VectorXd &r, VectorXd &x)
{
  int MyPID = _mat.GetPID();
  std::ios_base::fmtflags origflag; 
  origflag = std::cout.flags();
  if (_verbose) {
    std::cout.setf(std::ios::scientific);
  }
  double alpha, res0, delta_new, delta_old, beta,norm;
  int i = 0;
  //VectorXd r = b;
  VectorXd d(_mat.GetNrDofs());  
  _mat.Apply(x,d);
  //residual
  r = r - d;
  _prec.Solve(r,d);
  VectorXd s = r;
  res0 =  delta_new = _mat.dot(r,d);
  res0 = sqrt(res0);
  double resreal0 = sqrt(_mat.dot(r,r));
  norm = resreal0;
  if(_verbose) {
	  PCOUT(MyPID, "Initial Residuum: "<< resreal0 << std::endl)
  }
  while ( i < _maxIter && norm/resreal0 > _tol) {
    if (_verbose && i % _freq == 0) {
      PCOUT(MyPID, "Iteration: "<< i << "\t Residuum precondioned: " << sqrt(delta_new)/res0)
      PCOUT(MyPID, "\t Residuum: " << norm/resreal0 << std::endl)
    }
    _mat.Apply(d,s);
    alpha = delta_new / _mat.dot(d,s);

    x = x + alpha*d;
    //r = r - alpha*q;
    r = r - alpha*s;
    _prec.Solve(r,s);
    delta_old = delta_new;
    delta_new = _mat.dot(r,s);

    beta = delta_new/delta_old;
    d = s + beta*d;
    i++;
    norm =  sqrt(_mat.dot(r,r));
  }
  if (_verbose) {
    PCOUT(MyPID, "Iteration: "<< i << "\t Residuum precondioned: " << sqrt(delta_new)/res0 << std::endl)
    PCOUT(MyPID, "Iteration: "<< i << "\t Residuum:              " << norm/resreal0 << std::endl);
    std::cout.flags(origflag); 
    }

  return i;
}

const std::string PCGSolver::Label() const
{
  return("Preconditioned CG-Solver\n");
}
