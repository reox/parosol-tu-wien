#!/usr/bin/python

"""
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
"""
# create the xmf file to the output of ParOsol. The xml file
# is needed by paraview

import h5py
import sys
import os
from xml.etree.ElementTree import Element, SubElement
from xml.etree import ElementTree
from xml.dom import minidom

def data_item(parent, dimensions, numbertype, precision, data):
    """
    Create a DataItem entry

    :param parent: ElementTree object to attach to
    :param dimensions: Dimensions of this Dataitem
    :param numbertype: type of numbers, e.g. Float, Int, ...
    :param precision: Precision of numbers
    :param data: path to HDF5 Resource
    """
    x = SubElement(parent, 'DataItem', {'Dimensions': ' '.join(map(str, dimensions)),
                                        'Format': 'HDF5',
                                        'NumberType': numbertype,
                                        'Precision': precision})
    x.text = data
    return x


def attribute(parent, name, atype, center, dimensions, numbertype, precision, data):
    """
    Create an Attribute entry in the XML

    :param name: Name of entry
    :param atype: type of the attribute, e.g. Scalar, Vector, Tensor, ...
    :param center: Center type, e.g. Node or Cell
    :param parent: ElementTree object to attach to
    :param dimensions: Dimensions of this Dataitem
    :param numbertype: type of numbers, e.g. Float, Int, Tensor, ...
    :param precision: Precision of numbers
    :param data: path to HDF5 Resource
    """
    a = SubElement(parent, 'Attribute', {'Name': name,
                                         'AttributeType': atype,
                                         'Center': center})
    data_item(a, dimensions, numbertype, precision, data)


def get_grid(parent, filepath):
    """
    Generate the Grid entry, which contains all the information of the H5 file.

    :param parent: parent ElementTree item
    :param nr_nodes: Number of nodes
    :param nr_elements: number of elements
    :param filename: h5 filename
    """
    with h5py.File(filepath,'r+') as f:
        nr_nodes, _ = f['Mesh']['Coordinates'].shape
        nr_elements, _ = f['Mesh']['Elements'].shape

        print("{}: nodes {}, elements {}".format(filepath, nr_nodes, nr_elements))

    # Split the name, so we have the correct path in the xmf file
    filename = os.path.basename(filepath)

    grid = SubElement(parent, 'Grid', {'Name': 'mesh', 'GridType': 'Uniform'})


    mesh_topo = SubElement(grid, 'Topology', {'TopologyType': 'Hexahedron',
                                              'NumberOfElements': str(nr_elements),
                                              'BaseOffset': "1"})
    data_item(mesh_topo, (nr_elements, 8), 'Int', "8", "{}:/Mesh/Elements".format(filename))

    geom = SubElement(grid, 'Geometry', {'GeometryType': 'XYZ'})
    data_item(geom, (nr_nodes, 3), 'Float', "4", "{}:/Mesh/Coordinates".format(filename))

    attribute(grid, 'Displacement', 'Vector', 'Node', (nr_nodes, 3), 'Float', "8", "{}:/Solution/Nodal displacements".format(filename))
    attribute(grid, 'Force', 'Vector', 'Node', (nr_nodes, 3), 'Float', "8", "{}:/Solution/Nodal forces".format(filename))
    attribute(grid, 'SED', 'Scalar', 'Cell', (nr_elements, 1), 'Float', "8", "{}:/Solution/SED".format(filename))
    attribute(grid, 'vonMises', 'Scalar', 'Cell', (nr_elements, 1), 'Float', "8", "{}:/Solution/VonMises".format(filename))
    attribute(grid, 'EFF', 'Scalar', 'Cell', (nr_elements, 1), 'Float', "8", "{}:/Solution/EFF".format(filename))
    attribute(grid, 'Strain', 'Tensor6', 'Cell', (nr_elements, 6), 'Float', "8", "{}:/Solution/Element strain".format(filename))
    attribute(grid, 'Stress', 'Tensor6', 'Cell', (nr_elements, 6), 'Float', "8", "{}:/Solution/Element stress".format(filename))
    attribute(grid, 'Material ID', 'Scalar', 'Cell', (nr_elements, 1), 'Float', "8", "{}:/Mesh/Material IDs".format(filename))

    return grid


if __name__ == "__main__":
    if len(sys.argv) < 2:
      print("usage: {} filename [filename ...]".format(sys.argv[0]))
      print("If you supply more than one filename, a time series will be "
            "created, if the filenames are like samename-xx.x.h5")
      sys.exit(0)

    xdmf = Element('Xdmf')
    xdmf.set('version', '2.2')
    domain = SubElement(xdmf, 'Domain')
    if len(sys.argv) == 2:
        filepath = sys.argv[1]
        get_grid(domain, filepath)
        outfilename = filepath.replace(".h5",".xmf")
    else:
        # Check if all filenames look the same...
        fnames = sys.argv[1:]

        # all files end to h5?
        if not all(map(lambda x: x.endswith(".h5"), fnames)):
            print("Not all files are h5 files?")
            sys.exit(1)

        # Test if all files have a "-" and a float:
        try:
            map(lambda x: float(x[:-3].rsplit("-", 1)[1]), fnames)
        except (ValueError, IndexError):
            print("the filenames does not match our schema!")
            sys.exit(1)

        # Test if all filenames are the same
        bname = list(map(lambda x: x.rsplit("-", 1)[0], fnames))
        if len(set(bname)) != 1:
            print("Not all filenames match each other!")
            sys.exit(1)

        coll = SubElement(domain, 'Grid', {'Name': 'Remodelling',
                                           'GridType': 'Collection',
                                           'CollectionType': 'Temporal'})
        for filepath in fnames:
            grid = get_grid(coll, filepath)
            SubElement(grid, 'Time', {'Value': filepath[:-3].rsplit("-", 1)[1]})

        outfilename = "{}_timeseries.xmf".format(bname[0])


    #write file
    with open(outfilename, 'w') as outfile:
        outfile.write(minidom.parseString(ElementTree.tostring(xdmf)).toprettyxml(indent="    "))
