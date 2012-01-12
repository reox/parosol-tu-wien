/*
 * ParFE: a micro-FE solver for trabecular bone modeling
 * Copyright (C) 2006, Cyril Flaig, Uche Mennel and Marzio Sala
 *
 * ParOSol: a parallel FE solver for trabecular bone modeling
 *  Copyright (C) 2011, Cyril Flaig
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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  
 * 02110-1301, USA.
 */


  //! Element stiffness matrix integration routine.

  /*! This subroutine computes the element stiffness matrix by numerical integration (Gauss).
      
  \param mat_prop
  (In) An array of material properties

  \param nprops
  (In) Number of material properties.

  \param nod
  (In) Number of nodes per element.

  \param ndof
  (In) Number of DOFs per element

  \param nip
  (In) Number of integration points for the numerical integration

  \param ndim
  (In) Number of dimensions

  \param nst
  (In) Size of the stress-strain matrix

  \param coord
  (In) An array containing the coordinates of the mesh vertices

  \param km
  (Out) The element stiffness matrix
  */
  
  void Stiffness_Matrix(const double* mat_prop, const int& nprops,
                        const int& nod, const int& ndof, const int& ndim,
                        const int& nip, const int& nst, const double* coord, double* km);

//! Element stress and strain computation routine.

  /*! This subroutine computes the element stress and strain at a given number of Gauss points.
      
  \param mat_prop
  (In) An array of material properties

  \param nprops
  (In) Number of material properties.

  \param nod
  (In) Number of nodes per element.

  \param ndof
  (In) Number of DOFs per element.

  \param nip
  (In) Number of Gauss points.

  \param ndim
  (In) Number of dimensions.

  \param nst
  (In) Size of the stress-strain matrix.

  \param coord
  (In) An array containing the coordinates of the mesh vertices.

  \param eld
  (In) An array containing the displacements of the element nodes.

  \param stress
  (Out) An array containing the element stresses at the Gauss points.
        One column per Gauss point.

  \param strain
  (Out) An array containing the element strains at the Gauss points.
        One column per Gauss point.

  \param sigma
  (Out) An array containing the sigma values at each Gauss point.

  \param theta
  (Out) An array containing the von theta values at each Gauss point.

  */
  
  void Element_Stress(const double* mat_prop, const int& nprops,
                      const int& nod, const int& ndof, const int& ndim,
                      const int& nip, const int& nst, const double* coord,
                      double* eld, double* strain, double* stress,
                      double* sigma, double* theta);

