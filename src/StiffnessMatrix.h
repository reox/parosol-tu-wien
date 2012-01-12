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



#ifndef STIFFNESSMATRIX_H
#define STIFFNESSMATRIX_H

#include <string>
#include <eigen2/Eigen/Core>
#include "Config.h"

USING_PART_OF_NAMESPACE_EIGEN



//! An Interface which describes a MatrixVector Product.


class StiffnessMatrix 
{
public:

  //! ElementByElementMatrix destructor
  virtual ~StiffnessMatrix() {};

   /** Compute the matrix vetctor product.
   * 
   * @param x input vector
   * @param y result vector
   * @return error value
   */
  virtual int   Apply (VectorXd &x, VectorXd &y) = 0;

  /** Name of the Operator
   * 
   * @return Name in an string
   */
  virtual const std::string  Label () const = 0;

   /** Gets the Diagonal of the operator.
   * 
   * @return diagonal
   */
  virtual VectorXd& Diagonal() = 0;

   /** Gets the number of degrees of freedom.
   * 
   * @return number of dofs
   */
  virtual t_index GetNrDofs() =0;

   /** Computes the dot product of two vectors according the communicator is used in the matrix operator.
   * 
   * @param a first vector
   * @param b second vector
   * @return dot product
   */
  virtual double dot(VectorXd &a, VectorXd &b)=0;
  
  /** Gets the process id
   * @return PID
   */
  virtual int GetPID()=0;
  
  /** initializes an vector with random values according the distribution of the mesh on the computing nodes.
   * @param x vector that is filled with random numbers.
   */
  virtual void SetVectorRandom(VectorXd &x)=0;
};

#endif /* STIFFNESSMATRIX_H */


