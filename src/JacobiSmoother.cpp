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


#include "JacobiSmoother.h"

int JacobiSmoother::Solve(Eigen::VectorXd &b, Eigen::VectorXd &x)
{
  int it = 0;
  if (_zeroStart) {
    x = _w*_idia.cwiseProduct(b);
    it = 1;
   }
 for(; it<_step;it++) {
    _mat.Apply(x,r);
    // x = x + w*M^-1 *r;
    // x = x + w*M^-1 (b - A*x);
    x += _w*_idia.cwiseProduct(b-r);
  }
 return 1;
}
