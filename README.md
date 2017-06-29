# PrePostProcess4OpenSees


## Usage 

```python
from FEMtools import FEMtools

# Example: select nodes at the base (z = 0) of a FEM model, 
# and apply OpenSees command 'fix' to all selected nodes 


# select all nodes and elements
FEM_tool.selectAll()

# select nodes by location, getNbyLoc for Cartesian, getNbyLocCyl for Cylindrical
FEM_tool.getNbyLocCyl(x = (r1,r2), y=(angle1,angle2), z = (0,0), tol=0.001) # x is radius, y is angle

# fix nodeTage 1 1 1 0, stored in baseBCfile.tcl
FEM_tool.setBCFix('baseBCfile',x=1,y=1,z=1,u=0) 

# plot selected nodes
selectedNodes = FEM_tool.selectedNodes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(selectedNodes[:,1], selectedNodes[:,2], selectedNodes[:,3])
plt.show()
```


<img src="figures/1.png" height="300px">
