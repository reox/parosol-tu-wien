## Visualizing the Results ##

ParOSol ships with some helper tools. These are located in the directory `tools`.
There is also a python script, `createxmf.py`, that generates a xmf file. The xmf
file describes how to read the h5 mesh file.

    createxmf.py sphere.h5

A file `sphere.xmf` will be placed next to the `.h5` file.
The newly generated file can be opened with [Paraview](http://paraview.org).

Paraview has many features and filters for post-processing of the data.

Please read the [Paraview documentation](https://www.paraview.org/documentation/) for more information!

