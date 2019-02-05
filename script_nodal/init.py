#!/usr/bin/python

#Imports
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math as math

#class
from input import Input

class Init(Input):
    def __init__(self):
        """
        Class used to initialize the different arrays
   
        Variables:

        x           : variable containing the x position of the nodes
        y           : variable containing the y position of the nodes 
        
        nodes       : |id of node|x position of node|y position of node|
        neig        : |id of node|neighbour 1|neighbour 2| ...

        """
        
        Input.__init__(self)
        
        self.x = np.linspace(0,self.Lx/2,self.Nptsx)
        self.y = np.linspace(0,self.Ly/2,self.Nptsy)

        self.nodes = np.zeroes((self.Nptsx*self.Nptsy+self.Nptsx*self.ntheta,3))
        self.neig = np.zeroes((self.Nptsx*self.Nptsy+self.Nptsx*self.ntheta,5))

    def domain(self):
        """
        Creates the shape of the tank and fills it with nodes
        """

        #Create the boudary of the tank
        self.boundary=np.zeros((self.Nptsy,2))

        for j in range(0,self.Nptsy):
            self.boundary[j,0]=self.y[j]
            if (self.y[j]<self.Lx/2):
                self.boundary[j,1]=self.Lx/2-math.sqrt((self.Lx/2)**2-(self.y[j]-self.Lx/2)**2)
            elif (self.y[j]>self.Ly-self.Lx/2):
                self.boundary[j,1]=self.Lx/2-math.sqrt(abs((self.Lx/2)**2-(self.y[j]-self.Ly+self.Lx/2)**2))
            else :
                self.boundary[j,1]=0

        #Fill the tank with nodes
        self.xnodes=np.zeros((self.Nptsy,self.Nptsx))
        self.ynodes=np.zeros((self.Nptsy,self.Nptsx))

        for i in range(0,self.Nptsx):
            for j in range(0,self.Nptsy):
                if (self.y[j]<self.Lx/2):# and self.y[j]>math.sqrt(abs((self.Lx/2-self.x[i]/2)**2-(self.y[j]-self.Lx/2)**2))-self.Lx/2 ) :
                    self.ynodes[j,i]=self.y[j]
                    self.xnodes[j,i]=self.x[i]
                elif (self.y[j]<self.Ly-self.Lx/2 and self.y[j]>self.Lx/2):
                    self.ynodes[j,i]=self.y[j]
                    self.xnodes[j,i]=self.x[i]
                elif (self.y[j]>self.Ly+self.Lx/2): 
                    self.ynodes[j,i]=self.y[j]
                    self.xnodes[j,i]=self.x[i]#self.Lx/2-math.sqrt(abs((self.Lx/2-self.x[i]/2)**2-(self.y[j]-self.Ly+self.Lx/2)**2))


        
        plt.figure()
        plt.plot(self.boundary[:,1],self.boundary[:,0],'or')
        plt.xlim(-0.5,1.5)
        plt.ylim(-1,1)
        plt.plot(self.xnodes,self.ynodes,'o')
      #  plt.ylim(8,10.5)
      #  plt.xlim(0.5,1.4)
        plt.show()

    def init_domain(self):
        """
        Initialise the domain before the computation


        Variables :

        mesh        : contain the x and y position of the nodes (mesh[:,0] x positions,                      mesh[:,1] y positions)

        temp        : array containing the temperature field [K]
        pres        : array containing the pressure field [Pa]

        U           : array containing the x velocity field [m/s]
        V           : array containing the y velocity field [m/s]

        Rx          : array containing the conduction thermal resistance Rj,i+1/2 [K.m/W]
        Ry          : array containing the conduction thermal resistance Rj+1/2,i [K.m/W]

        """
        
        #Creation of the mesh
        self.mesh = np.array((self.y,self.x))
        
        #Creation of the temperature array
        self.temp = np.zeros((self.Nptsy,self.Nptsx)) 
        self.temp[1:self.Nptsy-1,1:self.Nptsx-1] = self.tfluid_init 
        self.temp[0,:] = self.T1 ; self.temp[self.Nptsy-1,:] = self.T1
        self.temp[:,0] = self.T1 ; self.temp[:,self.Nptsx-1] = self.T1
        
        #Creation of the pressure array
        self.pres = np.zeros((self.Nptsy,self.Nptsx)) 
        self.pres[:,:] = self.pfluid_init 
        
        #Creation of the velocity fields
        self.U = np.zeros((self.Nptsy,self.Nptsx))
        self.V = np.zeros((self.Nptsy,self.Nptsx))
        self.U[:,:]=self.ufluid_init
        self.V[:,:]=self.vfluid_init

        #Creation of thermal resistance arrays
        self.Rx = np.zeros((self.Nptsy,self.Nptsx+1))
        self.Ry = np.zeros((self.Nptsy+1,self.Nptsx))
    
    def resistance(self):
        if self.cond :
            for j in range(1,self.Nptsy-1):
                for i in range(1,self.Nptsx):      
                    self.Rx[j,i]=self.Rx[j,i] + (self.x[i]-self.x[i-1])/(self.k*(self.y[j+1]-self.y[j-1])/2)
            for i in range(1,self.Nptsx):
                self.Rx[0,i]=self.Rx[0,i] + (self.x[i]-self.x[i-1])/(self.k*(self.y[1]-self.y[0]))
                self.Rx[self.Nptsy-1,i]=self.Rx[self.Nptsy-1,i] + (self.x[i]-self.x[i-1])/(self.k*(self.y[self.Nptsy-1]-self.y[self.Nptsy-2]))
 
            for j in range(1,self.Nptsy):
                for i in range(1,self.Nptsx-1):      
                    self.Ry[j,i]=self.Ry[j,i] + (self.y[j]-self.y[j-1])/(self.k*(self.x[i+1]-self.x[i-1])/2)
            for j in range(1,self.Nptsy):
                self.Ry[j,0]=self.Ry[j,0] + (self.y[j]-self.y[j-1])/(self.k*(self.x[1]-self.x[0]))
                self.Ry[j,self.Nptsx-1]=self.Ry[j,self.Nptsx-1] + (self.y[j]-self.y[j-1])/(self.k*(self.x[self.Nptsx-1]-self.x[self.Nptsx-2]))
        
        if self.conv :
            print()
        #Creation on thermal capacity array
 

    def update_parameters(self):
        """
        Computes the thermodynamic parameters of the fluid using correlations
        """
