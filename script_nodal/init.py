#!/usr/bin/python3

#Imports
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math as math
import scipy.spatial as sp

#class
from input import Input

class Init(Input):
    def __init__(self):
        """
        Class used to initialize the different arrays
   
        Variables:

        x           : variable containing the x position of the nodes
        y           : variable containing the y position of the nodes 
        
        nodes       : |id of node|y position of node|x position of node|
        neig        : |id of node|neighbor 1|neighbor 2| ...  if neighbor = -1 : boundary

        """
        
        Input.__init__(self)
        
        self.x = np.linspace(0,self.Lx/2,self.Nptsx)
        self.y = np.linspace(0,self.Ly,self.Nptsy)

        self.nodes = np.zeros((self.Nptsx*self.Nptsy+self.Nptsx*self.ntheta,3))
        self.neig = np.zeros((self.Nptsx*self.Nptsy+self.Nptsx*self.ntheta,5))

    def domain_cart(self):
        """
        Rectangle domain to test
        """

        for j in range(0,self.Nptsy):
            for i in range(0,self.Nptsx):
                #nodes
                self.nodes[i+j*self.Nptsx,0]=i+j*self.Nptsx
                self.nodes[i+j*self.Nptsx,1]=self.y[j]
                self.nodes[i+j*self.Nptsx,2]=self.x[i]
                #neighbor
                self.neig[i+j*self.Nptsx,0]=i+j*self.Nptsx
                if (j>0 and j<self.Nptsy-1 and i>0 and i<self.Nptsx-1):
                    self.neig[i+j*self.Nptsx,1]=i-1+j*self.Nptsx 
                    self.neig[i+j*self.Nptsx,2]=i+1+j*self.Nptsx
                    self.neig[i+j*self.Nptsx,3]=i+(j-1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,4]=i+(j+1)*self.Nptsx
          
                elif (j==0 and i>0 and i<self.Nptsx-1):
                    self.neig[i+j*self.Nptsx,1]=i-1+j*self.Nptsx 
                    self.neig[i+j*self.Nptsx,2]=i+1+j*self.Nptsx
                    self.neig[i+j*self.Nptsx,3]=i+(j+1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,4]=-1

                elif (j==self.Nptsy-1 and i>0 and i<self.Nptsx-1):
                    self.neig[i+j*self.Nptsx,1]=i-1+j*self.Nptsx 
                    self.neig[i+j*self.Nptsx,2]=i+1+j*self.Nptsx
                    self.neig[i+j*self.Nptsx,3]=i+(j-1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,4]=-1

                elif (i==0 and j>0 and j<self.Nptsy-1):
                    self.neig[i+j*self.Nptsx,1]=i+1+j*self.Nptsx 
                    self.neig[i+j*self.Nptsx,2]=i+(j-1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,3]=i+(j+1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,4]=-1

                elif (i==self.Nptsx-1 and j>0 and j<self.Nptsy-1):
                    self.neig[i+j*self.Nptsx,1]=i-1+j*self.Nptsx 
                    self.neig[i+j*self.Nptsx,2]=i+(j+1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,3]=i+(j-1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,4]=-1

        self.neig[0,1]=1; self.neig[0,2]=self.Nptsx;self.neig[0,3]=-1                            
        self.neig[self.Nptsx-1,1]=self.Nptsx-2; self.neig[self.Nptsx-1,2]=2*self.Nptsx-1;self.neig[self.Nptsx-1,3]=-1
      
        self.neig[self.Nptsx*(self.Nptsy-1),1]=self.Nptsx*(self.Nptsy-1)+1; self.neig[self.Nptsx*(self.Nptsy-1),2]=self.Nptsx*(self.Nptsy-2);self.neig[self.Nptsx*(self.Nptsy-1),3]=-1
        self.neig[self.Nptsx*self.Nptsy-1,1]=self.Nptsx*self.Nptsy-2; self.neig[self.Nptsx*self.Nptsy-1,2]=self.Nptsx*(self.Nptsy-1)-1;self.neig[self.Nptsx*self.Nptsy-1,3]=-1
    
        
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

        temp        : array containing the temperature field [K]
        pres        : array containing the pressure field [Pa]

        U           : array containing the x velocity field [m/s]
        V           : array containing the y velocity field [m/s]

        R          : array containing the conduction thermal resistance [K.m/W]
        C          : array containing the conduction thermal capacity


        """
        
        #Creation of the temperature array
        self.temp = np.zeros((self.Nptsx*self.Nptsy+self.Nptsx*self.ntheta))
        
        #Creation of the pressure array
        self.pres = np.zeros((self.Nptsx*self.Nptsy+self.Nptsx*self.ntheta)) 
        self.pres[:] = self.pfluid_init 
        
        #Creation of the velocity fields
        self.U = np.zeros((self.Nptsx*self.Nptsy+self.Nptsx*self.ntheta))
        self.V = np.zeros((self.Nptsx*self.Nptsy+self.Nptsx*self.ntheta))
        self.U[:]=self.ufluid_init
        self.V[:]=self.vfluid_init

        #Creation of thermal resistance arrays
        self.Rx = np.zeros(self.Nptsy*(self.Nptsx-1))
        self.Ry = np.zeros((self.Nptsy-1)*self.Nptsx)
        self.R = np.zeros((self.Nptsy*self.Nptsx, 5))#.reshape((len(self.nodes),1))
        
        #Creation of thermal capacity array
        self.C = np.zeros((self.Nptsy*self.Nptsx))
    
    def initemp_cart(self):
        for k in range(0,self.Nptsx*self.Nptsy+self.Nptsx*self.ntheta):
            if int(self.nodes[k,1])<self.Ly/2 :
                self.temp[k]=self.T1
            else :
                self.temp[k]=self.tfluid_init
			
    def resistance_cart(self):
        dx=self.nodes[1,2] - self.nodes[0,2]
        dy=self.nodes[self.Nptsx,1] - self.nodes[0,1]
        for idnode in range(self.Nptsy*self.Nptsx) :
            j=1
            self.R[idnode,0] = idnode
            while (j<5 and (int(self.neig[idnode,j]) != -1)):
                ng=int(self.neig[idnode,j])
                dxx=abs(self.nodes[ng,2] - self.nodes[idnode,2])
                dyy=abs(self.nodes[ng,1] - self.nodes[idnode,1])
                if (dxx < 1e-6) :
                    res = dy /(self.k * dx)
                else :
                    res= dx / (self.k * dy)
                self.R[idnode,j]= res
                j+=1




