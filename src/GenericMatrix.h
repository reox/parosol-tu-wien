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


#ifndef GENERICMATRIX_H
#define GENERICMATRIX_H

#include <iomanip>
#include "Config.h"
#include "Timing.h"
#include "StiffnessMatrix.h"
#include "BoundaryCondition.h"
#include "Toolbox.h"
#include "fem.h"
#include "OptLocalMatrix.h"

void decarbenz(double  *Y, double  *X);

USING_PART_OF_NAMESPACE_EIGEN

//! A class to perform an matrix-vector product for elasticity problem based on grid

/*! The GenericMatrix class is used to multiply a vector by the global stiffness
 * matrix, whithout ever assembling the matrix. A local stiffness matrix is
 * computed. This is applied element wise to the nodes corresponding with the
 * weight of the bone. The base grid is given by a template parameter
 */


template <class Grid>
class GenericMatrix: public StiffnessMatrix 
{
	public:
		GenericMatrix(Grid &grid);

		~GenericMatrix();

		/** The Application is done with respect to the boundary condition. The matrix is not assembled
		 *  This function needs communication.
		 *
		 *  @param [in] x input vector
		 *  @param [out] y output vector
		 */
		int  Apply (VectorXd &x, VectorXd &y);

		/** This application of the matrix will not apply the boundary condition. It uses the indefinite matrix.
		 *  This function needs communication.
		 *
		 *  @param [in] x input vector
		 *  @param [out] y output vector
		 */
		int Apply_NoResetBoundaries(VectorXd &x, VectorXd &y);

		inline virtual const std::string  Label () const
		{
			return("Matrix");
		}

		/// Prepares the data structure that is needed by apply
		int PrepareApply();

		//!Extract the Diagonal
		VectorXd& Diagonal();

		Grid& GetGrid() {
			return _grid;
		}

		BoundaryCondition & GetBC()
		{
			return _grid.bc;
		}

		
		//!Get the local Number of Dofs  
		t_index GetNrDofs()
		{
			return _mydofs;
		}

		/** Computes the dot product of the vectors based on the grid
		 *  @param a first input vector
		 *  @param b second input vector
		 */
		double dot(VectorXd &a, VectorXd &b)
		{
			return _grid.dot(a,b);
		}


		friend std::ostream& operator<<(std::ostream& stream, const GenericMatrix &mat)
		{
			return mat.print_prop(stream);
		}

		int GetPID() {
			return MyPID;
		}
		void SetVectorRandom(VectorXd &x)
		{
			x.setRandom(GetNrDofs());
			x.setRandom(GetNrDofs());
			_grid.Recv_import_Ghost(x);
			MPI_Barrier(MPI_COMM_WORLD);
			_grid.Send_import_Ghost(x);
			_grid.Wait_import_Ghost();
		}

		void PrintTimings() {
			t_timing ela = Timing.ElapsedTime("Apply");
			PCOUT(MyPID, "Applytime " << COUTTIME(ela) << std::endl)
		}


	private:
		void ResetBoundaries(VectorXd &y);
		void SetBoundaries(VectorXd &, VectorXd &) const;

		Matrix<double, 24, 24>  *_stiffnessmatrix;
		double* _matprop;
		t_index _mydofs;
		Grid& _grid;
		VectorXd *_dia;

		std::vector<double> _stored_disp;

		double _factor;

		//some print function
		std::ostream& print_mat(std::ostream& stream) const;
		std::ostream& print_prop(std::ostream& stream) const;

		int MyPID;
		Timer Timing;

};

	template <class Grid>
	GenericMatrix<Grid>::GenericMatrix(Grid &grid)
:_grid(grid), _dia(0), Timing(MPI_COMM_WORLD)
{
	MyPID = _grid.GetPID();

	_mydofs = grid.GetNrDofs();

	//set up the local stiffness matrix
	double GridDim[3];
	int Dimension = 3;
	const int NumMaterialProps = 2;
	_matprop = new double[2];
	_matprop[0] = 1000; //refernece value Emodule is linear
	_matprop[1] = 0.3;
	int NumNodesPerElement = 8;
	int NumDofsPerElement = 24;
	int NumIntegrationPoints = 8;
	int SSMatrixSize = 6;
	double *coord = new double[Dimension * NumNodesPerElement];
	grid.GetRes(GridDim);
	setcoord(GridDim,coord);

	_stiffnessmatrix = new Matrix<double, 24, 24>;
	double *tmpstiff = new double[24*24];

	Stiffness_Matrix(_matprop, NumMaterialProps,
			NumNodesPerElement, NumDofsPerElement,
			Dimension, NumIntegrationPoints,
			SSMatrixSize, coord, tmpstiff); //TODO: not sure if this works
	for( int i = 0; i < 24*24;i++)
		(*_stiffnessmatrix)[i]=tmpstiff[i];

	delete[] coord;
	delete[] tmpstiff;
  	double nu = 0.3;
	//if the optimized code is used, following factor has to be used. Else 1
	//_factor = (1-nu)/((1+nu)*(1-2*nu))*1000/(144*(1-nu))*GridDim[0];
	_factor = 1;

	PrepareApply();

}

	template <class Grid>
GenericMatrix<Grid>::~GenericMatrix()
{
	delete _stiffnessmatrix;
	delete[] _matprop;
	if (_dia)
		delete _dia;

}

	template <class Grid>
int  GenericMatrix<Grid>::Apply(VectorXd &x, VectorXd &y)
{
	//Store the values on boundary... replace it with 0
	ResetBoundaries(x);
	//Import needed degree of freedoms....

	//Apply the singular matrices
	Apply_NoResetBoundaries(x,y);

	//Write back the stored disp
	SetBoundaries(x,y);
	return 0;
}

	template <class Grid>
int GenericMatrix<Grid>::Apply_NoResetBoundaries(VectorXd &x, VectorXd &y)
{
	Timing.Restart("Apply");
	//fetch the nodes of the neighbours
	_grid.Recv_import_Ghost(x);
	_grid.Send_import_Ghost(x);


	//fetch the 24 values in pref and store store
	// res = K_e * xpref
	Matrix<double,24,1> xpref;
	Matrix<double,24,1> res;
	int nr_ele = 0;

	y.setZero(_mydofs);
	_grid.Wait_import_Ghost();
	_grid.Recv_export_Ghost();

	for(_grid.initIterateOverElements(); _grid.TestIterateOverElements(); _grid.IncIterateOverElements()){
		_grid.GetNodalDisplacementsOfElement(x, xpref);

		nr_ele++;
		//decarbenz(res.data(), xpref.data());
        res = (*_stiffnessmatrix)*xpref;
		_grid.SumInToNodalDisplacementsOfElement(y, res, _factor*_grid.GetElementWeight());
	}

	_grid.Send_export_Ghost(y);
	_grid.WaitAndCopy_export_Ghost(y);

    
	Timing.Stop("Apply");
	return 0;
}

template <class Grid>
VectorXd& GenericMatrix<Grid>::Diagonal()
{
	if (_dia == 0) {
		_grid.Recv_export_Ghost();
		VectorXd* diagonal = new VectorXd(_mydofs);
		diagonal->setZero(_mydofs);

		Matrix<double,24,1> res;
		for(_grid.initIterateOverElements(); _grid.TestIterateOverElements(); _grid.IncIterateOverElements()){
			for(int i = 0; i <24; i++) {
				res[i] = _stiffnessmatrix->operator()(i,i);
			}
			_grid.SearchIndexes();
			_grid.SumInToNodalDisplacementsOfElement(*diagonal,res, _grid.GetElementWeight());
		}
		_grid.Send_export_Ghost(*diagonal);
		_grid.WaitAndCopy_export_Ghost(*diagonal);
		_grid.Recv_import_Ghost(*diagonal);

		// apply fixed nodes on vec
		for (t_index i=0; i < _grid.bc.FixedNodes_Ind.size(); ++i) {
			//1.0 for diagonal, 0.0 for rhs
			diagonal->operator[](_grid.bc.FixedNodes_Ind[i]) = 1.0;
		}

		double d;
		double dmax = 0;
		for (t_index i=0; i < _mydofs; ++i) {
			if (dmax < diagonal->operator[](i));
			dmax = diagonal->operator[](i);
		}
		//compute the pseudo invers of the diagonal. If the entry is zero than set it also in 
		//the inverse to zero
		for (t_index i=0; i < _mydofs; ++i) {
			d = diagonal->operator[](i);
			diagonal->operator[](i) = 1.0 / d;

			if (d == 0)
				diagonal->operator[](i) = 0;
	}

	_grid.Send_import_Ghost(*diagonal);
	_grid.Wait_import_Ghost();
	//*diagonal = diagonal->cwise()*(*diagonal);

	_dia = diagonal;
}
return *_dia;

}

//Compute the indicies of the nodes affactet by the BC
	template <class Grid>
int GenericMatrix<Grid>::PrepareApply()
{
	//local indices mapping from fixed dofs to domain map is down by bc
	_stored_disp.resize(_grid.bc.FixedNodes_Ind.size());

	return 0;
}
template <class Grid>
void GenericMatrix<Grid>::SetBoundaries(VectorXd &y, VectorXd &x) const
{
	t_index i = 0;
	std::vector<t_index>::const_iterator it;
	
			for (it = _grid.bc.FixedNodes_Ind.begin(); it != _grid.bc.FixedNodes_Ind.end(); ++it) {
				x[*it] = _stored_disp[i++];
				y[*it] =  x[*it];
			}
}

//! Sets all boundary nodes (both fixed and restrained) to zero.
template <class Grid>
void GenericMatrix<Grid>::ResetBoundaries(VectorXd &y)
{
	t_index i =0;
	std::vector<t_index>::const_iterator it;
	for (it = _grid.bc.FixedNodes_Ind.begin(); it != _grid.bc.FixedNodes_Ind.end(); ++it) {
		_stored_disp[i++] = y[*it];
		y[*it] = 0.0;
	}
}

template <class Grid>
std::ostream& GenericMatrix<Grid>::print_mat(std::ostream& stream) const
{
	for (int i = 0;i < 24; i++) {
		for (int j = 0;j < 24; j++)
			stream << std::setiosflags(std::ios::fixed)<< std::setw( 7) << std::setprecision( 4 )<< (*_stiffnessmatrix)(i,j) << " ";
		stream << std::endl;
	}
	return stream;
}

template <class Grid>
std::ostream& GenericMatrix<Grid>::print_prop(std::ostream& stream) const
{
	double GridDim[3];
	_grid.GetRes(GridDim);
	t_coord GlobalDim[3];
	_grid.GetGlobalDim(GlobalDim);
	stream << "Matrix:\n";
	stream << "   used Grid Size: " << GlobalDim[0] << ", " << GlobalDim[1] << ", " << GlobalDim[2]  << "\n";
	stream << "   used Grid Resolution: " << GridDim[0] << ", " << GridDim[1] << ", " << GridDim[2]  << "\n";
	stream << "   used Material properties:   E: " << _matprop[0] << " Poisson: "<< _matprop[1] << std::endl;
	return stream;
}

#endif /* GENERICMATRIX_H */


