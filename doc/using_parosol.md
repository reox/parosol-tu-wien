# Synopsis #

    parosol [options] input-file

# Description #

ParOSol is a fully-parallel micro-FE analysis code based on octrees to
solve linear elasticity problems. ParOSol solves the problem directly on
3D images, which can be obtained by CT-scans. It can be run on a desktop
machine as well as on a cluster using [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface).

# Using ParOSol #

If no options are given ParOSol solves up to a tolerance of 1e–6 and uses a
maximum of seven levels. An input file has to be provided. ParOSol writes
the result back into the input file once the solution is complete.

A sample command is given below:

    parosol --level 3 --tol 1e-4 simple.h5

This solves the elasticity equation on mesh given in the file
`simple.h5`. It uses maximal three multigrid levels and solves
with a relative residual of 1e–4.

If ParOSol was compiled using an MPI C++ compiler, you can use `mpirun` to let
ParOSol run in parallel:

    mpirun -n 4 parosol --level 3 --tol 1e-4 simple.h5

will run ParOSol using four threads.
