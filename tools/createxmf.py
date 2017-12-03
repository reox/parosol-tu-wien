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

if len(sys.argv) != 2:
  print("usage: "+sys.argv[0]+" filename")
  sys.exit(0)

filepath = sys.argv[1]
try:
  f = h5py.File(filepath,'r+')
except:
  print("usage: "+sys.argv[0]+" filename")
  sys.exit(2)

# Split the name, so we have the correct path in the xmf file
filename = os.path.basename(filepath)

#get numbers
nr_nodes, _ = f['Mesh']['Coordinates'].shape
nr_elements, _ = f['Mesh']['Elements'].shape

print("nodes: "+str(nr_nodes)+" elements: "+str(nr_elements))

def data_item(parent, dimensions, numbertype, precision, data):
    x = SubElement(parent, 'DataItem', {'Dimensions': ' '.join(map(str, dimensions)),
                                        'Format': 'HDF5',
                                        'NumberType': numbertype,
                                        'Precision': precision})
    x.text = data
    return x

def attribute(parent, name, atype, center, dimensions, numbertype, precision, data):
    a = SubElement(parent, 'Attribute', {'Name': name,
                                         'AttributeType': atype,
                                         'Center': center})
    data_item(a, dimensions, numbertype, precision, data)

xdmf = Element('Xdmf')
xdmf.set('version', '2.2')
domain = SubElement(xdmf, 'Domain')
grid = SubElement(domain, 'Grid', {'Name': 'mesh', 'GridType': 'Uniform'})


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

#write file
outfilename = filepath.replace(".h5",".xmf")
with open(outfilename, 'w') as outfile:
    outfile.write(minidom.parseString(ElementTree.tostring(xdmf)).toprettyxml(indent="    "))
