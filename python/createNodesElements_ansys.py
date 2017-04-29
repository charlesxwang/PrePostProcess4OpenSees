#!/usr/bin/env python
# This script is used to
# 1. comvert Gmsh file to paraview file
# 2. create nodes.tcl and elements.tcl for opensees

#vtk_writer = VTK_XML_Serial_Unstructured()
#vtk_writer.snapshot("filename.vtu", x, y, z, optional arguments...)
#vtk_writer.writePVD("filename.pvd")
from vtktools import VTK_XML_Serial_Unstructured
import numpy as np
import pandas as pd
from FEMtools import FEMtools

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



# ------------------------------------------------------------
# raw mesh file -> elements.txt & nodes.txt
# ------------------------------------------------------------

nodes = np.loadtxt('../data/nodes.txt')

print(nodes.shape, nodes.dtype)


elements = np.loadtxt('../data/elements.txt')
elements = elements + 1; # imported from Paraview, which strats from 0, so +1 to start from 1



print(elements.shape,np.array_str(elements.astype(int))[1:-1])

'''
# apply bc 1
with open('../data/bc1_raw.txt', 'r') as myfile:
    bc1=myfile.read().replace('\n', '')
bc1 = np.asarray(bc1.split(',')[0:-1]).astype(int)
bc1_df = pd.DataFrame(['fix']*len(bc1), columns=['cmd'])
bc1_df['node'] = np.asarray(bc1)
bc1_matrix = np.ones((len(bc1), 3))
bc1_matrix = np.append(bc1_matrix, np.zeros((len(bc1),1)), axis=1)
bc1_xyz_df = pd.DataFrame(bc1_matrix.astype(int), columns=['x','y','z','u'])
bc1_df[bc1_xyz_df.columns] = bc1_xyz_df
#bc1_df.drop(bc1_df.index[[len(bc1)-1]], inplace=True) # trim the last row
bc1_df.to_csv('../data/bc1.tcl',sep=' ',header=False,index=False)


# apply bc 2
with open('../data/bc2_raw.txt', 'r') as bc2file:
    bc2=bc2file.read().replace('\n', '')
bc2 = np.asarray(bc2.split(',')[0:-1]).astype(int)
bc2_df = pd.DataFrame(['fix']*len(bc2), columns=['cmd'])
bc2_df['node'] = np.asarray(bc2)
bc2_matrix_0 = np.zeros((len(bc2), 1))
bc2_matrix_1 = np.ones((len(bc2), 1))
bc2_matrix = np.append(bc2_matrix_0, bc2_matrix_1, axis=1)
bc2_matrix = np.append(bc2_matrix, bc2_matrix_0, axis=1)
bc2_matrix = np.append(bc2_matrix, bc2_matrix_1, axis=1)
bc2_xyz_df = pd.DataFrame(bc2_matrix.astype(int), columns=['x','y','z','u'])
bc2_df[bc2_xyz_df.columns] = bc2_xyz_df
#bc2_df.drop(bc2_df.index[[len(bc2)-1]], inplace=True) # trim the last row
bc2_df.to_csv('../data/bc2.tcl',sep=' ',header=False,index=False)
'''


# create nodes.tcl
nodes_df = pd.DataFrame(nodes[:,1:], columns=['x', 'y', 'z'])
nodes_df.insert(0, 'index', nodes[:,0].astype(int))
nodes_df.insert(0, 'cmd', ['node']*len(nodes[:,0]))
nodes_df.to_csv('../data/nodes.tcl',sep=' ',header=False,index=False)
with open('../data/nodes.tcl', 'a') as nodes_Tclfile:
    nodes_Tclfile.write('set numNodes '+str(len(nodes[:,0])))

if len(elements.shape) == 1:
    elements = elements.reshape([1,len(elements)])
elementsCount = len(elements[:,0])
# create elements.tcl
elements_df = pd.DataFrame(elements.astype(int))
elements_df.insert(0, 'index', range(1,elementsCount+1))
elements_df.insert(0, 'type', ['brickUP']*elementsCount)
elements_df.insert(0, 'cmd', ['element']*elementsCount)
for i in ['$matTag', '$Bfluid', '$rhoF', '$perm1', '$perm2', '$perm3', '$gravityX', '$gravityY', '$gravityZ']:
    elements_df[i] = [i]*elementsCount
elements_df.to_csv('../data/elements.tcl',sep=' ',header=False,index=False)
with open('../data/elements.tcl', 'a') as nodes_Tclfile:
    nodes_Tclfile.write('set numElems '+str(int(elementsCount)))


'''
# create equalDOF.tcl
equalDOF_df = pd.DataFrame(['3']*(len(bc2)-1), columns=['u'])
equalDOF_df.insert(0, 'y', ['1']*(len(bc2)-1))
equalDOF_df.insert(0, 'slave', bc2[1:])
equalDOF_df.insert(0, 'master', [str(bc2[0])]*(len(bc2)-1))
equalDOF_df.insert(0, 'cmd', ['equalDOF']*(len(bc2)-1))
equalDOF_df.to_csv('../data/equalDOF.tcl',sep=' ',header=False,index=False)
'''




x = nodes[:,1]
y = nodes[:,2]
z = nodes[:,3]
nodes_xyz = nodes[:,1:4]
#print(nodes_xyz.shape)

xmax = x.max()
xmin = x.min()
ymax = y.max()
ymin = y.min()
zmax = z.max()
zmin = z.min()
print(xmin,xmax,ymin,ymax,zmin,zmax)

radius = 8.0
height = 8.0


elements_fem = elements

#print(elements)
if len(elements.shape) == 2: # more than one element
    #elements = elements[:,[0, 1, 3, 2, 4, 5, 7, 6]]
    elements = np.concatenate(elements)
else: # only one element
    #elements = elements[[0, 1, 3, 2, 4, 5, 7, 6]]
    print(elements)

#elements = elements - 1



vtk_writer = VTK_XML_Serial_Unstructured()
#vtk_writer.snapshot("FEs.vtu", x, y, z, elements = elements, nodes = nodes_xyz,NumberOfComponentsResults="3",vtkType = "11",offsets="8")
disp_x = x*0
disp_y = y*0
disp_z = z*0
pwp = z*0

vtk_writer.snapshot("FEs.vtu", x+disp_x, y+disp_y, z+disp_z, elements = elements, nodes = nodes_xyz, disp=np.zeros(3*len(x)), NumberOfComponentsResults="3",vtkType = "12",offsets="8",pwp=pwp)







# FEM preprocess

FEM_tool = FEMtools(nodes,elements_fem)
'''
FEM_tool.apply_BC('bbbb')
print(FEM_tool.BCs)
selected = FEM_tool.getNbyLoc(x=(0,2)).getNbyLoc(z=(0,0.024))
print('kkk',FEM_tool.selectedNodes.shape)
'''

print('kkk',FEM_tool.selectedNodes.shape)

# fix the base
FEM_tool.selectAll()
print('kkk',FEM_tool.selectedNodes.shape)
FEM_tool.getNbyLocCyl(x = (0.00,radius), y=(0,360), z = (0,0), tol=0.001) # x is radius, y is angle
FEM_tool.setBCFix('baseBCfile',x=1,y=1,z=1,u=0)

selectedNodes = FEM_tool.selectedNodes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(selectedNodes[:,1], selectedNodes[:,2], selectedNodes[:,3])
plt.show()

# fix the top
FEM_tool.selectAll()
print('kkk',FEM_tool.selectedNodes.shape)
FEM_tool.getNbyLocCyl(x = (0.00,radius), y=(0,360), z = (height,height), tol=0.001) # x is radius, y is angle
FEM_tool.setBCFix('topBCfile',x=1,y=1,z=0,u=1)

selectedNodes = FEM_tool.selectedNodes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(selectedNodes[:,1], selectedNodes[:,2], selectedNodes[:,3])
plt.show()

# fix others 1 1 0 0
FEM_tool.selectAll()
print('kkk',FEM_tool.selectedNodes.shape)
FEM_tool.getNbyLocCyl(x = (0.00,radius), y=(0,360), z = (0.0011,height-0.0011), tol=0.001) # x is radius, y is angle
FEM_tool.setBCFix('othersBCfile',x=1,y=1,z=0,u=0)

selectedNodes = FEM_tool.selectedNodes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(selectedNodes[:,1], selectedNodes[:,2], selectedNodes[:,3])
plt.show()

# apply load at top
FEM_tool.selectAll()
print('kkk',FEM_tool.selectedNodes.shape)
FEM_tool.getNbyLocCyl(x = (0.00,radius), y=(0,360), z = (height,height), tol=0.001) # x is radius, y is angle
FEM_tool.applyCMD('LoadFile',cmd='sp',option='3 1.0')

selectedNodes = FEM_tool.selectedNodes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(selectedNodes[:,1], selectedNodes[:,2], selectedNodes[:,3])
plt.show()



'''
# fix the sides
FEM_tool.selectAll()
FEM_tool.getNbyLocCyl(x = (radius,radius), y=(0,360), z = (0.0+0.0011,height), tol=0.001) # x is radius, y is angle
FEM_tool.setBCFix('sideBCfile',x=1,y=1,z=0,u=0)

selectedNodes = FEM_tool.selectedNodes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(selectedNodes[:,1], selectedNodes[:,2], selectedNodes[:,3])
plt.show()

# fix the symmetric section at y=0
FEM_tool.selectAll()
FEM_tool.getNbyLoc(x = (0,radius-0.0011), y=(0,0), z = (0.0+0.0011,height), tol=0.001) # x is radius, y is angle
FEM_tool.setBCFix('sideY0BCfile',x=0,y=1,z=0,u=0)

selectedNodes = FEM_tool.selectedNodes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(selectedNodes[:,1], selectedNodes[:,2], selectedNodes[:,3])
plt.show()

# fix the symmetric section at x=0
FEM_tool.selectAll()
FEM_tool.getNbyLoc(y = (0,radius-0.0011), x=(0,0), z = (0.0+0.0011,height), tol=0.001) # x is radius, y is angle
FEM_tool.setBCFix('sideX0BCfile',x=1,y=0,z=0,u=0)

selectedNodes = FEM_tool.selectedNodes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(selectedNodes[:,1], selectedNodes[:,2], selectedNodes[:,3])
plt.show()
'''


'''
# fix the pwp for pile as zero
FEM_tool.selectAll()
FEM_tool.getNbyLocCyl(x = (0.0,0.0225), y=(0,360), z = (0.025,0.05), tol=0.001) # x is radius, y is angle
FEM_tool.setBCFix('pileBCfile',x=0,y=0,z=0,u=1)

selectedNodes = FEM_tool.selectedNodes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(selectedNodes[:,1], selectedNodes[:,2], selectedNodes[:,3])
plt.show()
'''

'''
# apply pore presure loades
FEM_tool.selectAll()
FEM_tool.getNbyLocCyl(x = (0.0225,0.1025), y=(0,360), z = (0.05,0.05), tol=0.001) # x is radius, y is angle
#FEM_tool.setBCFix('pileBCfile',x=0,y=0,z=0,u=1)
FEM_tool.writePwpLoad('pwpLoadNodes')

selectedNodes = FEM_tool.selectedNodes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(selectedNodes[:,1], selectedNodes[:,2], selectedNodes[:,3])
plt.show()

print(selectedNodes[:,0])
print(selectedNodes.shape)
'''
