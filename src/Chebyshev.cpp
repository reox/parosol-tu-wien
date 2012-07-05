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



#include "Chebyshev.h"
#include "est_ev.h"
#include "Jacobi.h"

Chebyshev::Chebyshev(StiffnessMatrix &M, int degree, bool ZeroStart, Eigen::VectorXd &m1, Eigen::VectorXd &m2, int maxit, double lmax, double lmin, double ratio):
  _mat(M),
  _idia(M.Diagonal()),
  _ldofs(M.GetNrDofs()),
  _degree(degree),
  _lmax(lmax),
  _lmin(lmin),
  _zeroStart(ZeroStart),
  z(m1),
  p(m2)
{
 if (maxit != 0) {
  Jacobi prec(_mat);
  eig e = est_ev(_mat, prec, maxit,z,p);
  _lmax = e.large*1.05;
  _lmin = e.small;
  //use ration if set
  if (ratio >0)
	  _lmin = e.large/ratio;
  if (0 == _mat.GetPID())
      std::cout << "lmax = " << _lmax << "   lmin = " << e.small << "   lmax used = " << _lmin  << std::endl;
  }
}

Chebyshev::Chebyshev(Chebyshev &old, bool ZeroStart):
  _mat(old._mat),
  _idia(old._idia),
  _ldofs(old._ldofs),
  _degree(old._degree),
  _lmax(old._lmax),
  _lmin(old._lmin),
  _zeroStart(ZeroStart),
  z(old.z),
  p(old.p)
{
}

int Chebyshev::Solve(Eigen::VectorXd &b, Eigen::VectorXd &x)
{
//  std:: cout << "solve with:\n";
//  std::cout << "lmax = " << _lmax << "   lmin = " << _lmin << endl;
  
  double alpha = _lmin;
  double beta = _lmax;
  double delta=2.0/ (beta - alpha);
  double theta=(beta + alpha)/2;
  double s1 = theta *delta;
  double invTheta = 1.0 / theta;

  if (_zeroStart) {
    p = _idia.cwise()*b*invTheta;
    x = p;
  } else {
    _mat.Apply(x,z);
    p = _idia.cwise() *( b -z)*invTheta;
    x = x + p;
  }

  double rhok = 1.0/s1, rhokp1;
  double dtmp1, dtmp2;
  for (int i = 1; i < _degree; i++)
  {
    _mat.Apply(x,z);
    rhokp1 = 1.0/ (2.0*s1 -rhok);
    dtmp1 = rhokp1 * rhok;
    dtmp2 = 2.0 * rhokp1 * delta;
    rhok = rhokp1;
    p = dtmp1 * p;
    p = p  + _idia.cwise() *  (b - z) * dtmp2;
    x = x + p;
  }

 return 1;
}
