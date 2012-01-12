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

#ifndef PCGSOLVER_H
#define PCGSOLVER_H

#include "Solver.h"

#include "StiffnessMatrix.h"

//! This class implements an PCG-method. 

/*! The preconditioned conjugate gradient method is a solver for symmetric
 *  positive definite linear system. The implementation is optimized in respect
 *  to low memory usage. 
*/



class PCGSolver: public Solver
{
  public:

      /** 
       * @brief Constructor of the PCG solver.
       * If you use it to apply the inverse of the matrix choose zeroStart to be true (example as Preeconditioner)
       * 
       * @param M Stiffnessmatrix
	   * @param precon Preconditioner that has to implement the solver interface.
	   * @param tol Relative tolerance for solving.
	   * @param maxIter Maximal number of iteration of the iterative solving process.
	   * @param output If this parameter is set to true, the solver outputs the residual
	   * @param freq If output set it writes freq-th iteration the residual 
	   *
	   */
    PCGSolver(StiffnessMatrix &M, Solver &precon, double tol=1e-10, int maxIter=500, bool output=true, int freq = 4);
    virtual ~PCGSolver();

    /**
     * @brief Solves the linear system \f[Ax=b\f] iteratively with the initial guess \f[x\f]
     * @param[in] b RHS
     * @param[in,out] x initial guess and solution
     */
    virtual int Solve(VectorXd &b, VectorXd &x);
    virtual const std::string  Label () const;

  protected:
    StiffnessMatrix &_mat;
    Solver &_prec;
    long _ldofs;
    long _maxIter;
    double _tol;
    bool _verbose;
	int _freq;
};

#endif /* PCGSOLVER_H */
