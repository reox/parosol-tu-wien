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


#ifndef MLLEVELCG_H 
#define MLLEVELCG_H 

#include "Solver.h"
#include "StiffnessMatrix.h"

#include "Jacobi.h"
#include "PCGSolver.h"

/*! This class wraps a PCG solver for coarse level solver */
class MlLevelCG: public Solver
{
  public:
    MlLevelCG(StiffnessMatrix &M):prec(M),solver(M,prec,1e-7,20,false)
    {
    }
    virtual ~MlLevelCG()
    {
    }

    virtual int Solve(VectorXd &b, VectorXd &x) {
      return solver.Solve(b,x);
    }

    const std::string Label() const
    {
      return "Coarsegrid PCG solver\n";
    }
  private:
    Jacobi prec;
    PCGSolver solver;

};
#endif /* MLLEVELSOLVERCG_H  */
