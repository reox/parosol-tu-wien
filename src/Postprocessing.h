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

#ifndef POSTPROCESSING_H
#define POSTPROCESSING_H

#include "Config.h"
#include "StiffnessMatrix.h"
#include "Toolbox.h"
#include "fem.h"


/*! The PostProcess Class computes the von mises stress and the strain energy density
 *  from the resulting displacement of the computation.
 */


template <class Grid>
class PostProcess 
{
	public:
		PostProcess(Grid &grid) :_grid(grid)
	{
	}

		//! ElementByElementMatrix destructor
		~PostProcess() { }

		void ComputeStressAndStrain(Eigen::VectorXd &disp, Eigen::VectorXd &VonMises, Eigen::VectorXd &SED, Eigen::VectorXd &eff) {
			//fetch the nodes of the neighbours
			_grid.Recv_import_Ghost(disp);
			_grid.Send_import_Ghost(disp);

            // Set the reference element up
			double GridDim[3];
			int Dimension = 3;
			const int NumMaterialProps = 2;
			double _matprop[2];
			_matprop[0] = 1000; //reference value Emodule is linear
			_matprop[1] = 0.3;
			int NumNodesPerElement = 8;
			int NumDofsPerElement = 24;
			int NumGaussPoints = 1;
			int SSMatrixSize = 6;
			double *coord = new double[Dimension * NumNodesPerElement];
			_grid.GetRes(GridDim);
			setcoord(GridDim,coord);

			double* strainbuf = new double[(SSMatrixSize +1) * NumGaussPoints ];
			double* stressbuf = new double[(SSMatrixSize +1) * NumGaussPoints ];
			double sigma, theta;
			t_index nr_elem =_grid.GetNrElem();
			VonMises.setZero(nr_elem);
			SED.setZero(nr_elem);
			eff.setZero(nr_elem);


			//fetch the 24 values in pref and store store
			// res = K_e * xpref
            Eigen::Matrix<double,24,1> xpref;
			int nr_ele = 0;

			_grid.Wait_import_Ghost();
			MPI_Barrier(MPI_COMM_WORLD);

			for(_grid.initIterateOverElements(); _grid.TestIterateOverElements(); _grid.IncIterateOverElements()){

				_grid.GetNodalDisplacementsOfElement(disp, xpref);

				_matprop[0] = 1000*_grid.GetElementWeight();
                double emoduli = 1000*_grid.GetElementWeight();
				Element_Stress(_matprop, NumMaterialProps,
						NumNodesPerElement, NumDofsPerElement,
	                    Dimension, NumGaussPoints, SSMatrixSize,
						coord, &xpref[0], 
						strainbuf, stressbuf, 
						&sigma, &theta);
				VonMises[nr_ele] = stressbuf[6];
				SED[nr_ele] = strainbuf[6];
                eff[nr_ele] = sqrt(2*SED[nr_ele]/emoduli); 

				nr_ele++;

			}
			delete[] strainbuf;
			delete[] stressbuf;
		}

	private:
		Grid& _grid;
};

#endif /* POSTPROCESSING_H */
