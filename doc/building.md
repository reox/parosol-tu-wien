# Getting the Source Code #

The source code of this fork is available at
[GitHub](https://github.com/reox/parosol-tu-wien).
The original source code is available at
[bitbucket](https://bitbucket.org/cflaig/parosol/overview).

To clone the repository use the following command:

    git clone https://github.com/reox/parosol-tu-wien

The source tree consists of the following directories:

* `build` directory to build ParOSol. It includes a sample configuration file
* `doc` holds the user's guide
* `src` holds all source code
* `tools` includes some non documented helper tools for generating meshes.

When extracting the file `mesh.tar.bz2` a directory will be created
that includes some example meshes.

# Building and running ParOSol #

This section describes how to build ParOSol from source and to run an
example mesh.


## Build Debian Package ##

If you are using a Debian based distribution, the easiest way to build ParOSol
is by creating a .deb package.
The description to build a debian package is available in the source code.

You need at least these packages installed:

    apt install build-essential devscripts

Then you can simply run

    mk-build-deps --install debian/control

to gather all other dependencies.
Now start the build by calling

    dpkg-buildpackage -us -uc

to create an unsigned package.

## Dependencies ##

ParOSol has minimal dependencies on other libraries. This simplifies
the installation process. To build ParOSol 
[HDF5 library 1.8](http://www.hdfgroup.org/HDF5/) and [eigen
v3](http://eigen.tuxfamily.org/index.php?title=Main_Page) has to be
installed. An MPI-compiler has to be used. OpenMPI and MPICH2 are
supported.
If you are using an MPI compiler, please use the HDF5 library that is supplied
by the compiler.

## Compiling ParOSol ##

If all libraries are installed from the Linux distribution, ParOSol can
be built in following way:

    cd build/
    ./do_configure
    make

The MPI compiler can be set with the following variable:

    -DCMAKE_CXX_COMPILER:FILEPATH=mpic++

If HDF5 is not installed in the standard library path then set the following
variables:

    HDF5_LIBRARIES:PATH=
    HDF5_INCLUDE_DIRS:PATH=

Eigen is a pure template library. This means the library consists of header
files only. This can be easily included by the header path switch (-I
for _GCC_ ). If eigen3 is installed set following variable:

    CMAKE_CXX_FLAGS:STRING=-I/path/toeigen

In the directory `build/` there is a sample configure script called
`do_configure.sh`. The script can be executed with following commands:

    cd build
    ./do_configure.sh

Instead of using the `do_configure.sh` script, **cmake** can be called
directly. If the **HDF5** library is installed at `/path/to/hdf5/` the
following command should be used to configure ParOSol:

    cmake \
        -D CMAKE_CXX_COMPILER:FILEPATH=mpic++ \
        -D CMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
        -D HDF5_LIBRARIES:PATH=/path/to/hdf5/lib/libhdf5.so\
        -D HDF5_INCLUDE_DIRS:PATH=/path/to/hdf5/include \
        ../src

After running the configure script, ParOSol can be built by simply running
the following command:

    make

Now the software is built in the `build/` directory.

## Building the Source Code Documentation ##

The source code of ParOSol is documented with
[Doxygen](http://www.doxygen.org). The creation of the documentation is
included in the build system. If Doxygen is installed and found by
cmake you can run following command to build the documentation:

    make docu


