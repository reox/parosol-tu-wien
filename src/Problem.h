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

#ifndef PROBLEM_H
#define PROBLEM_H

#include "GenericMatrix.h" 
#include "BoundaryCondition.h"
#include "Solver.h"



//! This class holds the linear system 

/*! Problem stores the Matrix, RHS and LHS. It also computes the the
 *  LHS out of the boundary conditions
 */

template <class Grid>
class Problem
{
	public:
		//!Problem constructor

		Problem(GenericMatrix<Grid> &mat):
			 _ldofs(mat.GetNrDofs()), _mat(mat), _bcond(mat.GetBC())
		{

			//allocate the LHS and RHS. Value are set through the BC
			_x = new Eigen::VectorXd(_ldofs);
			_b = new Eigen::VectorXd(_ldofs);
		}

		//! Problem destructor
		~Problem()
		{
			delete _x;
			delete _b;
		}

		//! Impose boundary condition
		int Impose();
		
		/** Set the solver
		 * @param s Implemented Solver
		 */
		void SetSolver(Solver &s)
		{
			solver = &s;
		}

		//! Start the solver
		int Solve() {
            Eigen::VectorXd tmp = *_b; // The solver may change b
			return solver->Solve(tmp,*_x);
		}
		long _ldofs;

		/** Gets the solution
		 * @return a refence to the solution vector
		 */
        Eigen::VectorXd & GetSol() {
			return *_x;
		}
		
		/** Gets the residual vector
		 * @return a refence to a the residual vector
		 */
        Eigen::VectorXd & GetRes() {
            Eigen::VectorXd *tmp = new Eigen::VectorXd(*_x);
			_mat.Apply(*_x, *tmp);
			*tmp = *_b -*tmp;
			return *tmp;
		}
		
		/** computes the true residuum
		 * @return the residuum
		 */
		double Res() {
            Eigen::VectorXd tmp =*_x;
			_mat.Apply(*_x, tmp);
			tmp = *_b -tmp;
			return sqrt(_mat.dot(tmp,tmp));
		}

		
				
		/** @name Some Debug function */
		//@{
  		/** Compares a vector with the solution
		 * @param x vector that is compared to the solution of the problem
		 * @return the norm of the difference
		 */
		double CompareSol(Eigen::VectorXd &x) {
            Eigen::VectorXd tmp = x - *_x;
			return sqrt(_mat.dot(tmp,tmp));
		}
		//! Print solutioin vector
		void PrintSol() {
			std::cout << "\nSolution\n";
			_mat.GetGrid().PrintVector(*_x);
		}

		//! Print the RHS
		void PrintRHS() {
			std::cout << "\n RHS \n";
			_mat.GetGrid().PrintVector(*_b);
		} //@}

	private:
		//! Helper function which set BC on the RHS
		int SetBoundaryConditions(Eigen::VectorXd &rhs);
		GenericMatrix<Grid> &_mat;
		BoundaryCondition &_bcond;
        Eigen::VectorXd *_x, *_b;
		Solver *solver;
};

template <class Grid>
int Problem<Grid>::Impose() {

	//renaming
    Eigen::VectorXd &x = *_x;
    Eigen::VectorXd &b = *_b;
	x.setZero(_ldofs);
	//Store Loads in B
	//At the moment none load supported
	b.setZero(_ldofs);

	//apply fixed nodes on Matrix
	long num_ind = _bcond.FixedNodes_Ind.size();
	//write displacement in x
	//corr = A*fixed_disp -> load which is cause by the fixed nodes 
	for (long i=0; i<num_ind; ++i) {
		x[_bcond.FixedNodes_Ind[i]] = _bcond.FixedNodes[i];
	}

	//b = b - A*fixed
	_mat.Apply_NoResetBoundaries(x, b);
	b = -1*b;

    //add the loaded nodes
	num_ind = _bcond.LoadedNodes_Ind.size();
	for (long i=0; i<num_ind; ++i) {
		b[_bcond.LoadedNodes_Ind[i]] += _bcond.LoadedNodes[i];
	}

	SetBoundaryConditions(b);
	//x.setRandom(_ldofs);
	//x.setConstant(_ldofs, 0.02);
	x.setZero(_ldofs);
	num_ind = _bcond.FixedNodes_Ind.size();
	for (long i=0; i<num_ind; ++i)
		x[_bcond.FixedNodes_Ind[i]] = _bcond.FixedNodes[i];
	return 0;
}

template <class Grid>
int Problem<Grid>::SetBoundaryConditions(Eigen::VectorXd &rhs)
{
	long num_ind = _bcond.FixedNodes_Ind.size();
	// write the disp in the vector 
	for (long i=0; i<num_ind; ++i)
		rhs[_bcond.FixedNodes_Ind[i]] = _bcond.FixedNodes[i];
	return 0;
}

#endif /* PROBLEM_H */

