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

#ifndef MLCYCLE_H
#define MLCYCLE_H

#include "Solver.h"
#include "JacobiSmoother.h"
#include "Chebyshev.h"
#include "OctreeGrid.h"
#include "MlOctreeGrid.h"
#include "GenericMatrix.h"
#include "MlLevelCG.h"
#include "MlLevelCheb.h"

#include "Testtools.h"

#include <eigen2/Eigen/Core>

//! A class to perform an MG-Cycle 

/*! This is the MG-Algorithm. Attention: The matrix must be prepared (Must have
 * included the boundaries condition)

*/

//#define DEBOUT

template <typename T>
class MlCycle : public Solver
{
	public:
       /** 
       * @brief Constructor of the MG-Algorithm.
       * The matrix must be prepared. This means it must have included the boundaries.
       * 
       * @param finegrid The grid of the finest level.
	   * @param mat The matrix of the finest level.
	   * @param degree Degree of the smoother
	   * @param level Number of max level.
	   * @param stype Smoother type (0: Chebyshev, 1: Jacobismoother)
	   * @param ratio Ration of the Eigenvlaue. See Chebyshev smoother
       * @param w weight parameter of the Jacobismoother
	   */
		MlCycle(OctreeGrid<T> &finegrid, GenericMatrix<OctreeGrid<T> > &mat,
			    int degree =8, int level = -1, int stype =0, double ratio=16, double w=0.61);

		~MlCycle()
		{
			delete _coarsegrid;
			delete _coarsemat;
			delete _solver;
			delete presmoother;
			delete postsmoother;
		}

		int Solve(Eigen::VectorXd& b, Eigen::VectorXd& x);

		const std::string Label() const
		{
			return "Multilevel solver";
		}


	private:

		OctreeGrid<T> &_finegrid;
		GenericMatrix<OctreeGrid<T> > &_finemat;

		int _smoothingsteps;
		int _level;
		double _w;

        Eigen::VectorXd r_fine, e_fine, r_coarse, e_coarse;

		MlOctreeGrid<T> * _coarsegrid;
		GenericMatrix<OctreeGrid<T> > * _coarsemat;
		Solver *_solver;
		Solver *presmoother;
		Solver *postsmoother;
};

	template <class T>
MlCycle<T>::MlCycle(OctreeGrid<T> &finegrid, GenericMatrix<OctreeGrid<T> > &mat, int degree, int level, int stype, double ratio, double w): _finegrid(finegrid), _finemat(mat), _level(level), _w(w)
{
//Debug output
//	if (mat.GetPID() ==0) {
//		COUT << "ML Level: " << level << " degree: " << degree << " w: " << w << "\n";
//	}
	//Create coarser matrix and grid
	_coarsegrid = new MlOctreeGrid<T>(&finegrid);
	_coarsegrid->GenerateBC();
	_coarsemat = new GenericMatrix<OctreeGrid<T> >(*_coarsegrid);


	//Get local DOFS
	t_index nr_fine_dofs = _finemat.GetNrDofs();
	t_index nr_coarse_dofs = _coarsemat->GetNrDofs();
	
	r_fine.setZero(nr_fine_dofs);
	e_fine.setZero(nr_fine_dofs);
	r_coarse.setZero(nr_coarse_dofs);
	e_coarse.setZero(nr_coarse_dofs);

	// setting up cheb
	int MaxIters = 10;

	if (stype ==0) {
		Chebyshev *tmpsmoother = new Chebyshev(_finemat,degree ,true, r_fine, e_fine, MaxIters, 0, 0, ratio);
		presmoother = tmpsmoother;
		postsmoother = new Chebyshev(*tmpsmoother, false);
	} else {
		JacobiSmoother *tmpsmoother = new JacobiSmoother(_finemat,degree+level ,w, true, r_fine);
		presmoother = tmpsmoother;
		postsmoother = new JacobiSmoother(*tmpsmoother, false);
	}
	if (mat.GetPID() ==0) {
		std::cout << *_coarsemat;
	}

	t_coord ld[3];
	_coarsegrid->GetLocalDim(ld);
	_coarsegrid->GetGlobalDim(ld);
	if (level == 0|| (ld[0] == 2) || (ld[1] == 2)|| (ld[2] == 2)) 
		_solver = new MlLevelCG(*_coarsemat);
	else
		_solver = new MlCycle(*_coarsegrid,  *_coarsemat, degree, level-1, stype, ratio, w);

	return;
}


	template <class T>
int MlCycle<T>::Solve(Eigen::VectorXd& b, Eigen::VectorXd& x)
{
	long nrdofs = _finemat.GetNrDofs();
#ifdef DEBOUT
	PrintVector<OctreeKey_Lookup> fine_converter(_finegrid);
	PrintVector<OctreeKey_Lookup> coarse_converter(*_coarsegrid);
    Eigen::VectorXd debr(nrdofs);
	x.setZero(nrdofs);
	_finemat.Apply(x, debr);
	debr = b - debr;
	COUT << "(level " << _level << ") " << "Res before smoothing " << debr.norm() << std::endl;
	COUT << "(level " << _level << ") " << "b  before smoothing " << debr.norm() << std::endl;
#endif

	presmoother->Solve(b,x);
#ifdef DEBOUT
	_finemat.Apply(x, debr);
	debr = b - debr;
	COUT << "(level " << _level << ") " << "Res after smoothing " << debr.norm() << std::endl;
#endif

	// W cycle ( cy == 2)
	for (int cy = 0; cy < 2; cy++) {
		//r = b - Ax
		_finemat.Apply(x,r_fine);
		r_fine = b - r_fine;
#ifdef DEBOUT
	COUT << "(level " << _level << ") " << "Res in solve after 1. smoothing " << r_fine.norm() << std::endl;
#endif

		//r_k = restrict(r)
		r_coarse.setZero(_coarsemat->GetNrDofs());
		_coarsegrid->Restrict(r_fine,r_coarse);
#ifdef DEBOUTV
		COUT << "r_fine: \n";
		fine_converter.PrintVectorOct(r_fine,3);
		COUT << "\n\n\nr_coarse: \n";
		coarse_converter.PrintVectorOct(r_coarse,3);
#endif
		e_coarse.setZero(_coarsemat->GetNrDofs());

		//Solve the coarse size problem
		_solver->Solve(r_coarse, e_coarse);

#ifdef DEBOUT
	COUT << "(level " << _level << ") " << "e_coarse before prolongate " << e_coarse.norm() << std::endl;
#endif 
		//x = x + Prolongate(e)
		//Prolongate_Interpolation(*ecoarse,*rfine);
		e_fine.setZero(nrdofs);
		_coarsegrid->Prolongate(e_coarse,e_fine);

#ifdef DEBOUT
	COUT << "(level " << _level << ") " << "e_fine before prolongate " << e_fine.norm() << std::endl;
#endif 
#ifdef DEBOUTV
		COUT << "XXXXXX: \n";
		fine_converter.PrintVectorOct(x,3);
		COUT << "\n\n\ne_coarse: \n";
		coarse_converter.PrintVectorOct(e_coarse,3);
		COUT << "\n\n\ne_fine: \n";
		fine_converter.PrintVectorOct(e_fine,3);
		COUT << "\n\n\n";
#endif
		
		x= x + e_fine;

#ifdef DEBOUT
	_finemat.Apply(x, debr);
	debr = b - debr;
	COUT << "(level " << _level << ") " << "Res before post " << debr.norm() << std::endl;
#endif

		postsmoother->Solve(b,x);

#ifdef DEBOUT
	_finemat.Apply(x, debr);
	debr = b - debr;
	COUT << "(level " << _level << ") " << "Res after post " << debr.norm() << std::endl;
#endif
	}
	for(long i=0; i < _finemat.GetNrDofs(); i++) {
		if (_finemat.Diagonal()[i] == 0)
			x[i] = 0;
	}

	return 1;
} 
#endif /* FULLGRIDMLCYCLE_H */

