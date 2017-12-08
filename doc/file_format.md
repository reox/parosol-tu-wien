# File Format #

ParOSol reads and writes into the same file. Thus the file consists of
an input and output part. Because HDF5 is used the file is organized in
datasets and group of datasets.

There are several different file formats and adjustments for special purposes.
This document describes the file format which is used by this version of
ParOSol!

## Data Types Involved ##

The HDF5 library handles the correct endianness of the hardware. The Table
below shows the data types used in the I/O operation. Since most
of the computers use the little endian format, a `double` is interpreted
in little endian.

`H5T_IEEE_F32LE` has to be read as IEEE floating point number in single
precision and little endian.

The table gives an overview on the used formats

  native      hdf5 type
  ---------   -------------------
  `float`     `H5T_IEEE_F32LE`
  `double`    `H5T_IEEE_F64LE`
  `short`     `H5T_STD_U16LE`
  `long`      `H5T_STD_I64LE`


## Input File Format ##

ParOSol reads all required data from the dataset group `Image_Data`.
This group must contain following datasets:

-   `/Image_Data/Image`
    -   Type: `H5T_IEEE_F32LE`
    -   Size: `(z,y,x)`
    -   Description: Holds the 3D image of size z,y,x. The values are the
        Young's modulus in MPa.

-   `/Image_Data/Voxelsize`
    -   Type: `H5T_IEEE_F64LE`
    -   Size: `(1, )`
    -   Description: Size of the voxels in mm.

-   `/Image_Data/Poison_ratio`
    -   Type: `H5T_IEEE_F64LE`
    -   Size: `(1, )`
    -   Description: The Poison’s ratio of the material `[0,0.5)`.
        Remark: Poison’s ratio is constant over the whole image.

Fixation and displacements of nodes are possible using the following two
datasets:

-   `/Image_Data/Fixed_Displacement_Coordinates`
    -   Type: `H5T_STD_U16LE`
    -   Size: `(k, 4)`
    -   Description: This datasets describes the boundary conditions. The
        first 3 values are the index of the boundary nodes. The
        fourth value describes the direction in which the node is fixed
        (`0 = x`, `1 = y`, `2 = z`).
        The index has to be given in `(z, y, x)` format.

-   `/Image_Data/Fixed_Displacement_Values`
    -   Type: `H5T_IEEE_F32LE`
    -   Size: `k`
    -   Description: This dataset holds the displacements of the boundary
        condition. `Fixed_Displacement_Coordinates` indicates the position
        and the direction of the nodes and this dataset saves theirs
        displacements. This dataset must have the same number of rows as
        `Fixed_Displacement_Coordinates`.

Some load might also be applied on the nodes. Therefore following two datasets
can optionally defined:

-   `/Image_Data/Loaded_Nodes_Coordinates`
    -   Type: `H5T_STD_U16LE`
    -   Size: `(l, 4)`
    -   Description: This datasets is optional. It describes which nodes
        have loads on them.  The first 3 values are the indices of the
        boundary nodes. The fourth value describes the direction in which
        the load acts (`0 = x`, `1 = y`, `2 = z`).
        The index has to be given in `(z, y, x)` format.

-   `/Image_Data/Loaded_Nodes_Values`
    -   Type: `H5T_IEEE_F32LE`
    -   Size: `l`
    -   Description: This dataset is optional. It is strongly related
        to `Loaded_Nodes_Coordinates`. `Loaded_Nodes_Coordinates` stores
        the position and the direction of the loaded nodes and this dataset
        stores the amount of the load per node. This dataset must have the
        same number of rows as `Loaded_Nodes_Coordinates`.

## Output file format ##

ParOSol writes the solution and the mesh. The mesh is needed to
visualize the solution. The solution and the mesh are two different
groups of datasets. 
The number of elements is equal to the number of non-zero voxel in the input file.
For a completly filled image, it would be `m = x * y * z`.
The number of nodes in a completly filled image is equal to `n = (x+1) * (y+1) * (z+1)`.
The number might differ for images with zero-valued voxels.

The `Mesh` group contains following datasets:

-   `/Mesh/Coordinates`
    -   Type: `H5T_IEEE_F32LE`
    -   Size: `(n, 3)`
    -   Description: This dataset holds the coordinates of the vertices.
        These are not the indices anymore, but the real coordinates of the
        nodes, given in the unit of `/Image_Data/Voxelsize`.
        These coordinates are now in `(x, y, z)` format.

-   `/Mesh/Elements`
    -   Type: `H5T_STD_I64LE`
    -   Size: `(m, 8)`
    -   Description: This dataset stores the 8 indices of the adjacency nodes
        of an element.
        Beware: the node indices start with 1 and not with 0!

-   `/Mesh/Material IDs`
    -   Type: `H5T_IEEE_F32LE`
    -   Size: `(m, 1)`
    -   Description: This dataset stores for each element the Young’s modulus.

The results and some post processing values are stored in the group
`Solution`:

-   `/Solution/Nodal displacements`
    -   Type: `H5T_IEEE_F64LE`
    -   Size: `(n, 3)`
    -   Description: Displacement vector at each node.

-   `/Solution/Nodal forces`
    -   Type: `H5T_IEEE_F64LE`
    -   Size: `(n, 3)`
    -   Description: Force vector at each node.

-   `/Solution/SED`
    -   Type: `H5T_STD_F64LE`
    -   Size: `(m, 1)`
    -   Description: Strain Energy Density. It is computed at the center of
        the elements.

-   `/Solution/VonMises`
    -   Type: `H5T_STD_F64LE`
    -   Size: `(m, 1)`
    -   Description: Von Mises Stress. It is computed at the center of the
        elements.

-   `/Solution/EFF`
    -   Type: `H5T_STD_F64LE`
    -   Size: `(m, 1)`
    -   Description: Effective strain. It is computed at the center of the
        elements.

-   `/Solution/Element strain`
    -   Type: `HST_STD_F64LE`
    -   Size: `(n, 6)`
    -   Description: symmetric components of strain tensor for each element.

-   `/Solution/Element stress`
    -   Type: `HST_STD_F64LE`
    -   Size: `(n, 6)`
    -   Description: symmetric components stress tensor for each element.


