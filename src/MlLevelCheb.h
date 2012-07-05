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

#ifndef MLLEVELCHEB_H 
#define MLLEVELCHEB_H 

#include "Solver.h"
#include "StiffnessMatrix.h"

#include "Chebyshev.h"

/*! This class wraps a Chebishev solver for coarse level solver */

class MlLevelCheb: public Solver
{
  public:
	/**
	 * @brief Construcutor of the Chebyshev solver for the corse size mesh. 
	 *
	 * @param M StiffnessMatrix that is solvd with the chebyshev solver.
	 * @param approx_it Number of iteration to estimate min/max eigenvalue.
	 */
    MlLevelCheb(StiffnessMatrix &M, int approx_it):a(M.GetNrDofs()),b(M.GetNrDofs()), solver(M,200, true, a,b,approx_it,0,0,-1)
    {
    }
    virtual ~MlLevelCheb()
    {
    }

    virtual int Solve(Eigen::VectorXd &b, Eigen::VectorXd &x) {
      return solver.Solve(b,x);
    }

    const std::string Label() const
    {
      return "Chebishev\n";
    }
  private:
    Eigen::VectorXd a,b;
    Chebyshev solver;

};
#endif /* MLLEVELSOLVERCHEB_H  */
