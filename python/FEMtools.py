#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class FEMtools:
    """
    USAGE:
    FEMtool = FEMtools()
    """
    def __init__(self,nodes=[],elements=[],testV=[]):
        self.fileNames = []
        self.BCs = []
        self.nodes = nodes
        self.elements = elements
        self.selectedNodes = np.copy(nodes)
    

    def coords_to_string(self, x,y,z):
        string = str()
        for i in range(len(x)):
            string = string + repr(x[i]) + ' ' + repr(y[i]) \
                    + ' ' + repr(z[i]) + ' '
        return string

    def array_to_string(self, a):
        string = str()
        for i in range(len(a)):
            string = string + repr(a[i]) + ' '
        return string

    def apply_BC(self, bc):
        self.BCs.append(bc)
        print(self.elements)
        return self
    
    def cart2pol(self,xy):
        x = xy[:,0]
        y = xy[:,1]
        rho = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x) * 180 / np.pi
        phi[np.where(phi < 0)] = phi[np.where(phi < 0)]+360.
        '''
        if phi < 0:
            phi = 360+phi
        '''
        return(rho, phi)

    def pol2cart(self,rp):
        rho = rp[:,0]
        phi = rp[:,1]
        phi = phi * np.pi / 180.
        x = rho * np.cos(phi)
        y = rho * np.sin(phi)
        return(x, y)

    def getNbyLoc(self,x=[],y=[],z=[],tol=0.00001):
        nodes = self.selectedNodes
        print('nodes: ',nodes.shape)
        try:
            print('xmin',x[0],'xmax',x[1])
            xloc = np.where((nodes[:,1]>=x[0] - tol) & (nodes[:,1]<=x[1] + tol))
            nodes = nodes[xloc]
        except:
            print("x location range not given.")
        try:
            print('ymin',y[0],'ymax',y[1])
            yloc = np.where((nodes[:,2]>=y[0] - tol) & (nodes[:,2]<=y[1] + tol))
            nodes = nodes[yloc]
        except:
            print("y location range not given.")
        try:
            print('zmin',z[0],'zmax',z[1])
            zloc = np.where((nodes[:,3]>=z[0] - tol) & (nodes[:,3]<=z[1] + tol))
            nodes = nodes[zloc]
        except:
            print("z location range not given.")
        print('nodes selected: ',nodes.shape)
        self.selectedNodes = nodes
        return self
        '''
        if len(nodes.shape) == 2: # more than one node
            selected_index = nodes[:,0]
        else: # only one element or no element
            selected_index = nodes[0]
        return selected_index
        '''

    def getNbyLocCyl(self,x=[],y=[],z=[],tol=0.00001): # x is radius, y is angle, z is height
        nodes = np.copy(self.selectedNodes)
        nodes[:,1], nodes[:,2] = self.cart2pol(nodes[:,1:3])
        print('nodes: ',nodes.shape)
        try:
            print('xmin',x[0],'xmax',x[1])
            xloc = np.where((nodes[:,1]>=x[0] - tol) & (nodes[:,1]<=x[1] + tol))
            nodes = nodes[xloc]
        except:
            print("x location range not given.")
        try:
            print('ymin',y[0],'ymax',y[1])
            yloc = np.where((nodes[:,2]>=y[0] - tol) & (nodes[:,2]<=y[1] + tol))
            nodes = nodes[yloc]
        except:
            print("y location range not given.")
        try:
            print('zmin',z[0],'zmax',z[1])
            zloc = np.where((nodes[:,3]>=z[0] - tol) & (nodes[:,3]<=z[1] + tol))
            nodes = nodes[zloc]
        except:
            print("z location range not given.")
        print('nodes selected: ',nodes.shape)
        nodes[:,1], nodes[:,2] = self.pol2cart(nodes[:,1:3])
        self.selectedNodes = np.copy(nodes)
        '''
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(nodes[:,1], nodes[:,2], nodes[:,3])
        plt.show()
        '''
        return self

    def selectAll(self):
        self.selectedNodes = self.nodes
        print('All nodes are selected')


    def setBCFix(self,BCname, nodes=[], x=[], y=[], z=[], u = []):
        if nodes==[]: nodes = self.selectedNodes
        if len(nodes.shape) == 2: # more than one element
            index = nodes[:,0]
        else: # only one node
            index = nodes[0]
        bc_df = pd.DataFrame(['fix']*len(index), columns=['cmd'])
        bc_df['node'] = np.asarray(index.astype(int))
        bc_matrix = np.ones((len(index), 1))*x
        bc_matrix = np.append(bc_matrix, np.ones((len(index),1))*y, axis=1)
        bc_matrix = np.append(bc_matrix, np.ones((len(index),1))*z, axis=1)
        bc_matrix = np.append(bc_matrix, np.ones((len(index),1))*u, axis=1)
        #bc_matrix = np.append(bc_matrix, np.zeros((len(index),1)), axis=1)
        bc_xyz_df = pd.DataFrame(bc_matrix.astype(int), columns=['x','y','z','u'])
        bc_df[bc_xyz_df.columns] = bc_xyz_df
        bc_df.to_csv('../data/'+BCname+'.tcl',sep=' ',header=False,index=False)


    def writePwpLoad(self,LoadName, nodes=[]):
        if nodes==[]: nodes = self.selectedNodes
        if len(nodes.shape) == 2: # more than one element
            index = nodes[:,0]
        else: # only one node
            index = nodes[0]
        bc_df = pd.DataFrame(["set"]*len(index), columns=['cmd'])
        
        bc_df['cmd2'] = ["layerThick("]*len(index)
        bc_df['IND'] = range(1,len(index)+1)
        bc_dfRP = pd.DataFrame([')']*len(index), columns=['RP'])
        bc_dfIndex = pd.DataFrame(index.astype(int), columns=['Index'])
        bc_df[bc_dfRP.columns] = bc_dfRP
        bc_df[bc_dfIndex.columns] = bc_dfIndex
        bc_df.to_csv('../data/'+LoadName+'.tcl',sep=' ',header=False,index=False)
        with open('../data/'+LoadName+'.tcl', 'a') as LoadName_Tclfile:
            LoadName_Tclfile.write('set numPwpLoad '+str(int(len(index))))

    def applyCMD(self,BCname, nodes=[], cmd=[], option=[]):
        if nodes==[]: nodes = self.selectedNodes
        if len(nodes.shape) == 2: # more than one element
            index = nodes[:,0]
        else: # only one node
            index = nodes[0]
        bc_df = pd.DataFrame([cmd]*len(index), columns=['cmd'])
        bc_df['node'] = np.asarray(index.astype(int))
        opt_matrix = pd.DataFrame([option]*len(index), columns=['opt'])
        bc_df[opt_matrix.columns] = opt_matrix
        bc_df.to_csv('../temp',sep=' ',header=False,index=False)
        with open('../temp', "rt") as fin:
            with open('../data/'+BCname+'.tcl', "wt") as fout:
                for line in fin:
                    fout.write(line.replace('"', ''))




#print(nodes_new.shape,self.elements.shape)





