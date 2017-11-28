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

# xml
sStart = "<?xml version=\"1.0\" ?>\n" \
         "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n" \
         "<Xdmf Version=\"2.0\">\n" \
         " <Domain>\n" \
         "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n" 

sMesh =  "     <Topology TopologyType=\"Hexahedron\" NumberOfElements=\""+repr(nr_elements)+"\" BaseOffset=\"1\" >\n" \
         "       <DataItem Dimensions=\""+repr(nr_elements)+" 8\" NumberType=\"Int\" Precision=\"8\" Format=\"HDF\">\n" \
         "         "+filename+":/Mesh/Elements\n" \
         "       </DataItem>\n" \
         "     </Topology>\n" \
         "     <Geometry GeometryType=\"XYZ\">\n" \
         "       <DataItem Dimensions=\""+repr(nr_nodes)+" 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n" \
         "         "+filename+":/Mesh/Coordinates\n" \
         "       </DataItem>\n" \
         "     </Geometry>\n"

sDisp =  "     <Attribute Name=\"Displacement\" AttributeType=\"Vector\" Center=\"Node\">\n" \
         "       <DataItem Dimensions=\""+repr(nr_nodes)+" 3\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" >\n" \
         "         "+filename+":/Solution/Nodal displacements\n" \
         "       </DataItem>\n" \
         "     </Attribute>\n"

sForce=  "     <Attribute Name=\"Force\" AttributeType=\"Vector\" Center=\"Node\">\n" \
         "       <DataItem Dimensions=\""+repr(nr_nodes)+" 3\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" >\n" \
         "         "+filename+":/Solution/Nodal forces\n" \
         "       </DataItem>\n" \
         "     </Attribute>\n"

sSED =   "     <Attribute Name=\"SED\" AttributeType=\"Scalar\" Center=\"Cell\">\n" \
         "       <DataItem Dimensions=\""+repr(nr_elements)+" 1\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" >\n" \
         "           "+filename+":/Solution/SED\n" \
         "       </DataItem>\n" \
         "     </Attribute>\n"

sVM =    "     <Attribute Name=\"S_vonMises\" AttributeType=\"Scalar\" Center=\"Cell\">\n" \
         "       <DataItem Dimensions=\""+repr(nr_elements)+" 1\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" >\n" \
         "           "+filename+":/Solution/VonMises\n" \
         "       </DataItem>\n" \
         "     </Attribute>\n"

eff =    "     <Attribute Name=\"EffectiveStrain\" AttributeType=\"Scalar\" Center=\"Cell\">\n" \
         "       <DataItem Dimensions=\""+repr(nr_elements)+" 1\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" >\n" \
         "           "+filename+":/Solution/EFF\n" \
         "       </DataItem>\n" \
         "     </Attribute>\n"

# FIXME size is nr_elements 8
eleStrain =    "     <Attribute Name=\"eleStrain\" AttributeType=\"Tensor6\" Center=\"Cell\">\n" \
         "       <DataItem Dimensions=\""+repr(nr_elements)+" 6\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" >\n" \
         "           "+filename+":/Solution/Element strain\n" \
         "       </DataItem>\n" \
         "     </Attribute>\n"

# FIXME size is nr_elements 7
eleStress =    "     <Attribute Name=\"eleStress\" AttributeType=\"Tensor6\" Center=\"Cell\">\n" \
         "       <DataItem Dimensions=\""+repr(nr_elements)+" 6\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" >\n" \
         "           "+filename+":/Solution/Element stress\n" \
         "       </DataItem>\n" \
         "     </Attribute>\n"

sMaterial ="     <Attribute Name=\"MaterialID\" AttributeType=\"Scalar\" Center=\"Cell\">\n" \
         "       <DataItem Dimensions=\""+repr(nr_elements)+" 1\" Format=\"HDF\" NumberType=\"Int\" Precision=\"4\" >\n" \
         "           "+filename+":/Mesh/Material IDs\n" \
         "       </DataItem>\n" \
         "     </Attribute>\n"

sEnd =   "   </Grid>\n" \
         " </Domain>\n" \
         "</Xdmf>\n" 
#write file
outfilename = filepath.replace(".h5",".xmf")
with open(outfilename, 'w') as outfile:
    outfile.write("\n".join([sStart,
                             sMesh,
                             sDisp,
                             sForce,
                             sSED,
                             sVM,
                             eff,
                             eleStrain,
                             eleStress,
                             sMaterial,
                             sEnd]))
