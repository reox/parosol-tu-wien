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


#ifndef JACOBI_H
#define JACOBI_H

#include "Solver.h"

#include "StiffnessMatrix.h"


//! This class provied a Jacobi preconditioner
class Jacobi: public Solver
{
	public:
		Jacobi(StiffnessMatrix &M): _idia(M.Diagonal())
		{
		}

		virtual ~Jacobi()
		{
		}

		virtual int Solve(VectorXd &b, VectorXd &x)
		{
			x = (_idia.cwise() * b).lazy();  
			return 0;
		}
		virtual const std::string  Label () const
		{
			return "Jacobi Preconditioner";
		}

	protected:
		VectorXd &_idia;
		long _ldofs;
		long _maxIter;
};

#endif /* JACOBI_H */
