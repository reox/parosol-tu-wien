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

#include <iostream>

#include <mpi.h>
#include <eigen3/Eigen/Core>


#include "Timing.h"
#include "Config.h"

//is chosen by the templates
#include "KeyGenerator.h"

#include "HDF5Image.h"
#include "OctreeGrid.h"
#include "GenericMatrix.h"
#include "PCGSolver.h"
#include "Problem.h"
#include "MlCycle.h"

#include "VTKPrinter.h"
#include "HDF5Printer.h"
#include "HDF5ParfePrinter.h"
#include "PfePrinter.h"
 
#define EXIT(X) MPI_Finalize(); exit(X);

//measure mflops with a matrix
template <class T>
void mflops(GenericMatrix<T> &matr, int MyPID);

template <class T>
void print(GenericMatrix<T> &matr, Problem<T> &problem, std::string file, int MyPID);

double global_poisson_ratio;

int main( int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int pid, psize;
	MPI_Comm_rank (MPI_COMM_WORLD, &pid);
	int MyPID = pid;
	MPI_Comm_size (MPI_COMM_WORLD, &psize);

	std::string file;
    int level = 6, degree = 10;
    double tol = 1e-6;

    std::string usage = std::string("usage: ") + argv[0]
                      + std::string(" [--level arg (6)]")
                      + std::string(" [--tol arg (1e-6)]")
                      + std::string(" filename\n"); 
	
    if (argc < 2) {
		PCOUT(MyPID, usage)
		MPI_Finalize();
		return -1;
	}

    int i = 1;
    std::string param;
    while(i < argc) {
        param = std::string(argv[i]);
        if (param.compare("--level") == 0) {
            i++;
            level = atoi(argv[i]);
        }
        else if(param.compare("--tol") == 0) {
            i++;
            tol = atof(argv[i]);
        } else {
            file = param;
        }
        i++;
    }

	PCOUT(MyPID, "file: " << file << " tolerance: " << tol << " max. num. level: "<< level << std::endl);
    if (file.compare("") == 0) {
		PCOUT(MyPID, "No file given\n");
		MPI_Finalize();
		return -10;
    }

	Timer timer(MPI_COMM_WORLD);
	timer.Start("All");

	t_timing elapsed_time;

	CPULayout layout;

    //Read in the Problem
	timer.Start("Read");
	HDF5Image ir(file, layout);
	typedef OctreeGrid<OctreeKey_Lookup> t_Ogrid;
	t_Ogrid grid(ir);
	timer.Stop("Read");
	elapsed_time = timer.ElapsedTime("Read");
	PCOUT(MyPID, "Time for Read in: "  << COUTTIME(elapsed_time) << "s\n");

	PCOUT(MyPID, "Nr. of global elements: " << grid.GetNrElemGlobal() << " global nodes: " << grid.GetNrNodesGlobal() << "\n");

    //Generate the BC
	timer.Start("BC");
	grid.GenerateBC();
	timer.Stop("BC");
	elapsed_time = timer.ElapsedTime("BC");
	PCOUT(MyPID, "Time for Construction of Boundary Condition: " << COUTTIME(elapsed_time) << "s\n");

	PCOUT(MyPID, grid);

    //Construct the matrix
	timer.Start("Mat");
	GenericMatrix<t_Ogrid> matr(grid);
	timer.Stop("Mat");
	elapsed_time = timer.ElapsedTime("Mat");
	PCOUT(MyPID, "Time for Construction of Matrix: " << COUTTIME(elapsed_time) << "s\n");
	PCOUT(MyPID, matr);

	std::cout << std::setprecision (13); 
	MPI_Barrier(MPI_COMM_WORLD);

    //Construct the preconditioner
	timer.Start("Prec");
	MlCycle<OctreeKey_Lookup>  prec(grid, matr, degree, level-1, 0, 16,0);
	PCOUT(MyPID, "Prec done\n");
	timer.Stop("Prec");
	elapsed_time = timer.ElapsedTime("Prec");
	PCOUT(MyPID, "Time for Building Preconditioner: " << COUTTIME(elapsed_time) << "s\n");

    //Setting up the problem
	timer.Start("Setup");
    int maxiter = 1000;
    int outputfreq = 4;
	PCGSolver solver(matr, prec, tol, maxiter, true, outputfreq);
	Problem<t_Ogrid> problem(matr);
	if (problem.Impose())
		std::cout << "Error in impose\n";
	problem.SetSolver(solver);
	timer.Stop("Setup");
	elapsed_time =  timer.ElapsedTime("Setup");
	PCOUT(MyPID, "Time for Setting up the Problem: " << COUTTIME(elapsed_time) << "s\n");

    //solving the system
	double resid = problem.Res();
	PCOUT(MyPID, "Residuum before solving: " << resid << std::endl);
	timer.Start("Solve");
	problem.Solve();
	resid = problem.Res();
	PCOUT(MyPID, "Residuum after solving: " << resid << std::endl);
	timer.Stop("Solve");
	t_timing time = timer.ElapsedTime("Solve");
	PCOUT(MyPID, "Solving time: " << COUTTIME(time) << "s\n");

    //print out the solution
    print(matr, problem, file, MyPID);

	timer.Stop("All");
	elapsed_time = timer.ElapsedTime("All");

	PCOUT(MyPID, "Overall time: " << COUTTIME(elapsed_time) << "s\n");

	//mflops(matr, MyPID);

	MPI_Finalize();
	return 0;
}

template <class T>
void print(GenericMatrix<T> &matr, Problem<T> &problem, std::string file, int MyPID) {
	Timer timer(MPI_COMM_WORLD);
	//VTKPrinter<OctreeKey_Lookup> print("out", matr.GetGrid());
	//HDF5Printer<OctreeKey_Lookup> print(file, matr.GetGrid());
	//HDF5ParfePrinter<OctreeKey_Lookup> print2(file, matr.GetGrid());
	PfePrinter<OctreeKey_Lookup> print3(file, matr.GetGrid());
	matr.PrintTimings();
	timer.Start("Write");
    int dofs = matr.GetGrid().GetNrDofs();
    Eigen::VectorXd force(dofs);
    matr.Apply_NoResetBoundaries(problem.GetSol(), force);
	//print.PrintAll(problem.GetSol(), force, problem.GetRes());
	//print2.PrintAll(problem.GetSol(), force, problem.GetRes());
	print3.PrintAll(problem.GetSol(), force, problem.GetRes());
	timer.Stop("Write");
	t_timing time = timer.ElapsedTime("Write");
	PCOUT(MyPID, "Outputtime: " << COUTTIME(time) << "s\n");
}

template <class T>
void mflops(GenericMatrix<T> &matr, int MyPID) {
	Timer timer(MPI_COMM_WORLD);
	int dofs = matr.GetGrid().GetNrDofs(); 
    Eigen::VectorXd a,b;
	a.setRandom(dofs);
	b.setZero(dofs);
	t_octree_key globalelem = matr.GetGrid().GetNrElemGlobal();
	int size = matr.GetGrid().GetNrCPU();
	int loops = 16e6/globalelem*size+1;
	double flop = globalelem*(24*(24+24)+48)*loops; //+48 the factor

	PCOUT(MyPID, "Meassure with " << loops << "loops" << std::endl);
	MPI_Barrier(MPI_COMM_WORLD);
	timer.Start("gflop");
	for(int i = 0; i <loops; i++)
		matr.Apply_NoResetBoundaries(a, b);
	timer.Stop("gflop");
	t_timing ela = timer.ElapsedTime("gflop");
	double elaps = ela.max; 
	PCOUT(MyPID, "time needed: " << elaps << std::endl)
	PCOUT(MyPID, "MFlop/s: " << flop/elaps/1e6<< std::endl)
}
