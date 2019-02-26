#!/usr/bin/python2.7

#Imports
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math as math
# ~ import scipy.spatial as sp

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
        neig        : |id of node|neighbor 1|neighbor 2| ...  
                      in neig : -1 = no more neighbors, -2 = wall, -3 = symetry
        """
        
        Input.__init__(self)
        
        self.x = np.linspace(0,self.Lx/2,self.Nptsx)
        self.y = np.linspace(self.Lx/2,self.Ly-self.Lx/2,self.Nptsy)

        if (self.mesh_type=='cart'):
            self.nodes = np.zeros((self.Nptsx*self.Nptsy,3))
            self.neig = np.zeros((self.Nptsx*self.Nptsy,6))
        else:
            self.nodes = np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta,3))
            self.neig = np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta,5+self.ntheta))


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
                    self.neig[i+j*self.Nptsx,5]=-1 
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

        self.neig[0,1]=1
        self.neig[0,2]=self.Nptsx
        self.neig[0,3]=-1                            
        
        self.neig[self.Nptsx-1,1]=self.Nptsx-2
        self.neig[self.Nptsx-1,2]=2*self.Nptsx-1
        self.neig[self.Nptsx-1,3]=-1
      
        self.neig[self.Nptsx*(self.Nptsy-1),1]=self.Nptsx*(self.Nptsy-1)+1
        self.neig[self.Nptsx*(self.Nptsy-1),2]=self.Nptsx*(self.Nptsy-2)
        self.neig[self.Nptsx*(self.Nptsy-1),3]=-1
        
        self.neig[self.Nptsx*self.Nptsy-1,1]=self.Nptsx*self.Nptsy-2
        self.neig[self.Nptsx*self.Nptsy-1,2]=self.Nptsx*(self.Nptsy-1)-1
        self.neig[self.Nptsx*self.Nptsy-1,3]=-1
    
        
    def domain_tank(self):
        """
        Creates the shape of the tank and fills it with nodes
        
        creation in 3 steps :
            -creation on the 'rectangle part' of the tank and storage of the neighbors (see domain_cart)
            -creation of the bottom 'circle part' centered on (Lx,Lx) using polar coordinates
            -creation of the upper 'circle part' centered on (Ly+Lx,Lx) using polar coordinates


        angle   : angle between two consecutive lines of points for polar coordinats

        local variables :
       
            id_node : id of the node (used to create the 'cicrle parts') 
            r       : distance to the center of the 'circle' for polar coordinates
            theta   : angle relative to the border of the 'rectangle part' for polar coordinates
        """
       
        assert (self.ntheta > 0), "ntheta must be >0 to use the tank shape" 

        #Create the 'rectangle part' of the tank

        for j in range(0,self.Nptsy):
            for i in range(0,self.Nptsx):
                #nodes
                self.nodes[i+j*self.Nptsx,0]=i+j*self.Nptsx
                self.nodes[i+j*self.Nptsx,1]=self.y[j]
                self.nodes[i+j*self.Nptsx,2]=self.x[i]
                #neighbor
                self.neig[i+j*self.Nptsx,0]=i+j*self.Nptsx
                #middle
                if (j>0 and j<self.Nptsy-1 and i>0 and i<self.Nptsx-1):
                    self.neig[i+j*self.Nptsx,1]=i-1+j*self.Nptsx 
                    self.neig[i+j*self.Nptsx,2]=i+1+j*self.Nptsx
                    self.neig[i+j*self.Nptsx,3]=i+(j-1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,4]=i+(j+1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,5]=-1 
                #lower boudary
                elif (j==0 and i>0 and i<self.Nptsx-1):
                    self.neig[i+j*self.Nptsx,1]=i-1+j*self.Nptsx 
                    self.neig[i+j*self.Nptsx,2]=i+1+j*self.Nptsx
                    self.neig[i+j*self.Nptsx,3]=-1 #temporary -1
                    self.neig[i+j*self.Nptsx,4]=i+(j+1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,5]=-1
                #upper boundary
                elif (j==self.Nptsy-1 and i>0 and i<self.Nptsx-1):
                    self.neig[i+j*self.Nptsx,1]=i-1+j*self.Nptsx 
                    self.neig[i+j*self.Nptsx,2]=i+1+j*self.Nptsx
                    self.neig[i+j*self.Nptsx,3]=i+(j-1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,4]=-1 #temporary -1
                    self.neig[i+j*self.Nptsx,5]=-1
                #left boundary
                elif (i==0 and j>0 and j<self.Nptsy-1):
                    self.neig[i+j*self.Nptsx,1]=-2
                    self.neig[i+j*self.Nptsx,2]=i+1+j*self.Nptsx 
                    self.neig[i+j*self.Nptsx,3]=i+(j-1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,4]=i+(j+1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,5]=-1
                #right boundary
                elif (i==self.Nptsx-1 and j>0 and j<self.Nptsy-1):
                    self.neig[i+j*self.Nptsx,1]=i-1+j*self.Nptsx 
                    self.neig[i+j*self.Nptsx,2]=-3
                    self.neig[i+j*self.Nptsx,3]=i+(j-1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,4]=i+(j+1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,5]=-1
        #corners
        self.neig[0,1]=-2
        self.neig[0,2]=1
        self.neig[0,3]=-1 #temporary -1
        self.neig[0,4]=self.Nptsx
        self.neig[0,5]=-1
                            
        self.neig[self.Nptsx-1,1]=self.Nptsx-2
        self.neig[self.Nptsx-1,2]=-3
        self.neig[self.Nptsx-1,3+self.ntheta]=2*self.Nptsx-1
        self.neig[self.Nptsx-1,4+self.ntheta]=-1

        self.neig[self.Nptsx*(self.Nptsy-1),1]=-2     
        self.neig[self.Nptsx*(self.Nptsy-1),2]=self.Nptsx*(self.Nptsy-1)+1
        self.neig[self.Nptsx*(self.Nptsy-1),3]=-1 #temporary -1
        self.neig[self.Nptsx*(self.Nptsy-1),4]=self.Nptsx*(self.Nptsy-2)
        self.neig[self.Nptsx*(self.Nptsy-1),5]=-1
        
        self.neig[self.Nptsx*self.Nptsy-1,1]=self.Nptsx*self.Nptsy-2
        self.neig[self.Nptsx*self.Nptsy-1,2]=-3
        self.neig[self.Nptsx*self.Nptsy-1,3]=self.Nptsx*(self.Nptsy-1)-1
        self.neig[self.Nptsx*self.Nptsy-1,4]=-1
   
        id_node=self.Nptsy*self.Nptsx-1 
        #Create lower 'cricle part'
       
        self.angle=math.pi/(self.ntheta*2)

        for theta in range(1,self.ntheta+1):
            for r in range(1,self.Nptsx):
                #update the id of the node
                id_node+=1
                #nodes
                self.nodes[id_node,0]=id_node
                self.nodes[id_node,1]=self.Lx/2-self.x[r]*math.sin(theta*self.angle)
                self.nodes[id_node,2]=self.Lx/2-self.x[r]*math.cos(theta*self.angle)
                #neighbor
                self.neig[id_node,0]=id_node
                if (r>1 and r<self.Nptsx-1 and theta>1 and theta<self.ntheta):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=id_node-1
                    self.neig[id_node,3]=id_node+self.Nptsx-1
                    self.neig[id_node,4]=id_node-self.Nptsx+1
                    self.neig[id_node,5]=-1
                #border of the circle
                elif (r==self.Nptsx-1 and theta>1 and theta<self.ntheta):
                    self.neig[id_node,1]=-2
                    self.neig[id_node,2]=id_node-1
                    self.neig[id_node,3]=id_node+self.Nptsx-1
                    self.neig[id_node,4]=id_node-self.Nptsx+1
                    self.neig[id_node,5]=-1
                #right boundary
                elif (r>1 and r<self.Nptsx-1 and theta==self.ntheta):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=id_node-1
                    self.neig[id_node,3]=-3
                    self.neig[id_node,4]=id_node-self.Nptsx+1
                    self.neig[id_node,5]=-1
                #upper boundary
                elif (r>1 and r<self.Nptsx-1 and theta==1 ):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=id_node-1
                    self.neig[id_node,3]=id_node+self.Nptsx-1
                    self.neig[id_node,4]=self.Nptsx-1-r
                    self.neig[id_node,5]=-1
                    self.neig[self.Nptsx-1-r,3]=id_node
                #center of the circle
                elif (r==1 and theta>1 and theta<self.ntheta):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=self.Nptsx-1
                    self.neig[id_node,3]=id_node+self.Nptsx-1
                    self.neig[id_node,4]=id_node-self.Nptsx+1
                    self.neig[id_node,5]=-1
                    self.neig[self.Nptsx-1,2+theta]=id_node
                #corners
                elif (r==self.Nptsx-1 and theta==1):
                    self.neig[id_node,1]=-2
                    self.neig[id_node,2]=id_node-1
                    self.neig[id_node,3]=id_node+self.Nptsx-1
                    self.neig[id_node,4]=0
                    self.neig[id_node,5]=-1
                    self.neig[0,3]=id_node
                elif (r==self.Nptsx-1 and theta==self.ntheta):
                    self.neig[id_node,1]=-2
                    self.neig[id_node,2]=id_node-1
                    self.neig[id_node,3]=-3
                    self.neig[id_node,4]=id_node-self.Nptsx+1
                    self.neig[id_node,5]=-1
                elif (r==1 and theta==1):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=self.Nptsx-1
                    self.neig[id_node,3]=id_node+self.Nptsx-1
                    self.neig[id_node,4]=self.Nptsx-2
                    self.neig[id_node,5]=-1
                    self.neig[self.Nptsx-2,3]=id_node
                    self.neig[self.Nptsx-1,3]=id_node
                elif (r==1 and theta==self.ntheta):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=self.Nptsx-1
                    self.neig[id_node,3]=-3
                    self.neig[id_node,4]=id_node-self.Nptsx+1
                    self.neig[id_node,5]=-1
                    self.neig[self.Nptsx-1,2+theta]=id_node
        
        # Create upper circle part
        for theta in range(1,self.ntheta+1):
            for r in range(1,self.Nptsx):
                #update the id of the node
                id_node+=1
                #nodes
                self.nodes[id_node,0]=id_node
                self.nodes[id_node,1]=self.Ly-self.Lx/2+self.x[r]*math.sin(theta*self.angle)
                self.nodes[id_node,2]=self.Lx/2-self.x[r]*math.cos(theta*self.angle)
                #neighbor
                self.neig[id_node,0]=id_node
                if (r>1 and r<self.Nptsx-1 and theta>1 and theta<self.ntheta):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=id_node-1
                    self.neig[id_node,3]=id_node-self.Nptsx+1
                    self.neig[id_node,4]=id_node+self.Nptsx-1
                    self.neig[id_node,5]=-1
                #border of the circle
                elif (r==self.Nptsx-1 and theta>1 and theta<self.ntheta):
                    self.neig[id_node,1]=-2
                    self.neig[id_node,2]=id_node-1
                    self.neig[id_node,3]=id_node-self.Nptsx+1
                    self.neig[id_node,4]=id_node+self.Nptsx-1
                    self.neig[id_node,5]=-1
                #right boundary
                elif (r>1 and r<self.Nptsx-1 and theta==self.ntheta):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=id_node-1
                    self.neig[id_node,3]=id_node-self.Nptsx+1
                    self.neig[id_node,4]=-3
                    self.neig[id_node,5]=-1
                #lower boundary
                elif (r>1 and r<self.Nptsx-1 and theta==1 ):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=id_node-1
                    self.neig[id_node,3]=self.Nptsy*self.Nptsx-1-r
                    self.neig[id_node,4]=id_node+self.Nptsx-1
                    self.neig[id_node,5]=-1
                    self.neig[self.Nptsy*self.Nptsx-1-r,3]=id_node
                #center of the circle
                elif (r==1 and theta>1 and theta<self.ntheta):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=self.Nptsx*self.Nptsy-1
                    self.neig[id_node,3]=id_node-self.Nptsx+1
                    self.neig[id_node,4]=id_node+self.Nptsx-1
                    self.neig[id_node,5]=-1
                    self.neig[self.Nptsy*self.Nptsx-1,3+theta]=id_node
                #corners
                elif (r==self.Nptsx-1 and theta==1):
                    self.neig[id_node,1]=-2
                    self.neig[id_node,2]=id_node-1
                    self.neig[id_node,3]=self.Nptsx*(self.Nptsy-1)
                    self.neig[id_node,4]=id_node+self.Nptsx-1
                    self.neig[id_node,5]=-1
                    self.neig[self.Nptsx*(self.Nptsy-1),3]=id_node
                elif (r==self.Nptsx-1 and theta==self.ntheta):
                    self.neig[id_node,1]=-2
                    self.neig[id_node,2]=id_node-1
                    self.neig[id_node,3]=id_node-self.Nptsx+1
                    self.neig[id_node,4]=-3
                    self.neig[id_node,5]=-1
                elif (r==1 and theta==1):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=id_node+self.Nptsx-1
                    self.neig[id_node,3]=self.Nptsy*self.Nptsx-2
                    self.neig[id_node,4]=self.Nptsy*self.Nptsx-1
                    self.neig[id_node,5]=-1
                    self.neig[self.Nptsy*self.Nptsx-2,4]=id_node
                    self.neig[self.Nptsy*self.Nptsx-1,4]=id_node
                elif (r==1 and theta==self.ntheta):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=self.Nptsy*self.Nptsx-1
                    self.neig[id_node,3]=id_node-self.Nptsx+1
                    self.neig[id_node,4]=-3
                    self.neig[id_node,5]=-1
                    self.neig[self.Nptsy*self.Nptsx-1,3+theta]=id_node
                    self.neig[self.Nptsy*self.Nptsx-1,4+theta]=-1
                    
    def init_phase_cart(self):
        """
        phi         : array describing whether the fluid is liquid or gas 
        """
        #Creation of the phase array
        self.phi = np.zeros((self.Nptsx*self.Nptsy))

        #Initialisation of the phase array
        self.phi[:]=0.

    def init_phase_tank(self):
        """
        phi         : array describing whether the fluid is liquid or gas 
        """
        #Creation of the phase array
        self.phi = np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta))

        #Initialisation of the phase array (phi = 1 => gas)
        self.phi[:]=0.

    def init_domain_cart(self):
        """
        Initialise the domain before the computation


        Variables :

        temp        : array containing the temperature field [K]
        pres        : array containing the pressure field [Pa]

        U           : array containing the x velocity field [m/s]
        V           : array containing the y velocity field [m/s]

        R           : array containing the conduction thermal resistance [K.m/W]
        C           : array containing the conduction thermal capacity

        """
        #Creation of the temperature array
        self.temp = np.zeros((self.Nptsx*self.Nptsy))
        
        #Creation of the pressure array
        self.pres = np.zeros((self.Nptsx*self.Nptsy)) 
        self.pres[:] = (self.phi[:]*self.pgas_init+(1-self.phi[:])*self.pliq_init)
        
        #Creation of the velocity fields
        self.U = np.zeros((self.Nptsx*self.Nptsy))
        self.V = np.zeros((self.Nptsx*self.Nptsy))
        self.U[:]=(self.phi[:]*self.ugas_init+(1-self.phi[:])*self.uliq_init)
        self.V[:]=(self.phi[:]*self.vgas_init+(1-self.phi[:])*self.vliq_init)

        #Creation of thermal resistance arrays
        self.R = np.zeros((self.Nptsy*self.Nptsx, 5))#.reshape((len(self.nodes),1))
        #~ np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta))
        
        #Creation of thermal capacity array
        self.C = np.zeros((self.Nptsy*self.Nptsx))
   

    def init_domain_tank(self):
        """
        Initialise the domain before the computation


        Variables :

        temp        : array containing the temperature field [K]
        pres        : array containing the pressure field [Pa]

        U           : array containing the x velocity field [m/s]
        V           : array containing the y velocity field [m/s]

        R           : array containing the conduction thermal resistance [K.m/W]
        C           : array containing the conduction thermal capacity

        """
       
        #Creation of the temperature array
        self.temp = np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta))
        
        #Creation of the pressure array
        self.pres = np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta)) 
        self.pres[:] = (self.phi[:]*self.pgas_init+(1-self.phi[:])*self.pliq_init) 
        
        #Creation of the velocity fields
        self.U = np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta))
        self.V = np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta))
        self.U[:]=(self.phi[:]*self.ugas_init+(1-self.phi[:])*self.uliq_init)
        self.V[:]=(self.phi[:]*self.vgas_init+(1-self.phi[:])*self.vliq_init)

        #Creation of thermal resistance arrays
        self.R = np.zeros((self.Nptsy*self.Nptsx, 5))#.reshape((len(self.nodes),1))
        #~ np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta))
        
        #Creation of thermal capacity array
        self.C = np.zeros((self.Nptsy*self.Nptsx+2*(self.Nptsx-1)*self.ntheta))
   
    def initemp_cart_y(self):
        for k in range(0,self.Nptsx*self.Nptsy):
            if self.nodes[k,1]<self.Ly/2 :
                self.temp[k]=self.T1
            else :
                self.temp[k]=(self.phi[k]*self.tgas_init+(1-self.phi[k])*self.tliq_init)
              
    def initemp_cart_x(self):
        for k in range(0,self.Nptsx*self.Nptsy):
            if self.nodes[k,2]<self.Lx/4:
                self.temp[k]=self.T1
            else :
                self.temp[k]=(self.phi[k]*self.tgas_init+(1-self.phi[k])*self.tliq_init)
			
    def initemp_tank_y(self):
        for k in range(0,self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta):
            if self.nodes[k,1]<self.Ly/2 :
                self.temp[k]=self.T1
            else :
                self.temp[k]=(self.phi[k]*self.tgas_init+(1-self.phi[k])*self.tliq_init)
                

    def initemp_tank_x(self):
        for k in range(0,self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta):
            if self.nodes[k,2]<self.Lx/4:
                self.temp[k]=self.T1
            else :
                self.temp[k]=(self.phi[k]*self.tgas_init+(1-self.phi[k])*self.tliq_init)
 
    def initemp_tank(self):
        for k in range(0,self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta):
            self.temp[k]=(self.phi[k]*self.tgas_init+(1-self.phi[k])*self.tliq_init)
    
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
                    res = dy /((self.phi[idnode]*self.k_gas+(1-phi[idnode])*self.k_liq)* dx)
                else :
                    res= dx / ((self.phi[idnode]*self.k_gas+(1-phi[idnode])*self.k_liq) * dy)
                self.R[idnode,j]= res
                j+=1

			
    def resistance_tank(self):
        print('Computing Thermal Resistance for the domain...')
        dx=self.nodes[1,2] - self.nodes[0,2]
        dy=self.nodes[self.Nptsx,1] - self.nodes[0,1]
        self.R=np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta, 4+self.ntheta))
        for idnode in range(self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta) :
            k_diph=(self.phi[idnode]*self.k_gas+(1-self.phi[idnode])*self.k_liq)
            #POINT CENTRE DU BAS
            if (idnode == self.Nptsx - 1) :
                self.R[idnode,0] = idnode
                #~ voisin gauche
                j=1
                ng=int(self.neig[idnode,j])
                self.R[idnode,j] = dx / (k_diph * dy)
                #~ voisins bas
                for j in range(3,self.ntheta+3) :
                    ng=int(self.neig[idnode,j])
                    rng=np.sqrt( (self.Lx-self.nodes[ng,1])**2+ (self.Lx-self.nodes[ng,2])**2 )
                    self.R[idnode,j] = rng / (k_diph * rng*np.sin(self.angle))
                j+=1
                #~ voisin haut
                ng=int(self.neig[idnode,j])
                self.R[idnode,j] = dx / (k_diph * dy)


            #POINT CENTRE DU HAUT
            elif (idnode == self.Nptsy*self.Nptsx - 1) :
                self.R[idnode,0] = idnode
                #~ voisin gauche
                j=1
                ng=int(self.neig[idnode,j])
                self.R[idnode,j] = dx / (k_diph * dy)
                #~ voisin bas
                j+=2
                ng=int(self.neig[idnode,j])
                self.R[idnode,j] = dx / (k_diph * dy)
                #~ voisins haut
                for j in range(4,self.ntheta+4) :
                    ng=int(self.neig[idnode,j])
                    rng=np.sqrt( (self.Ly-self.Lx-self.nodes[ng,1])**2+ (self.Lx-self.nodes[ng,2])**2 )
                    self.R[idnode,j] = rng / (k_diph * rng*np.sin(self.angle))


                
			#~ PARTIE RECTANGLE
            elif (idnode>=0 and idnode<self.Nptsx*self.Nptsy) :
                j=1
                self.R[idnode,0] = idnode
                while ((int(self.neig[idnode,j]) != -1)):
                    ng=int(self.neig[idnode,j])
                    if (ng != -3 and ng !=-2) :
                        if j<3:
                            res = dx /(k_diph * dy)
                        else :
                            dxx=abs(self.nodes[ng,2] - self.nodes[idnode,2])
                            dyy=abs(self.nodes[ng,1] - self.nodes[idnode,1])
                            l=np.sqrt(dxx**2 + dyy**2)
                            res= dy / (k_diph * dx)
                        self.R[idnode,j]= res
                    j+=1
                    
			#~ PARTIE POLAIRE
            else :
                rnode=min(math.sqrt( (self.Ly-self.Lx-self.nodes[idnode,1])**2+ (self.Lx-self.nodes[idnode,2])**2 ), math.sqrt( (self.Lx-self.nodes[idnode,1])**2+ (self.Lx-self.nodes[idnode,2])**2 ))
                j=1
                self.R[idnode,0] = idnode
                while ((int(self.neig[idnode,j]) != -1)):
                    ng=int(self.neig[idnode,j])
                    if (ng != -3 and ng !=-2) :
                        if j<3:
                            rng=min(math.sqrt( (self.Ly-self.Lx-self.nodes[ng,1])**2+ (self.Lx-self.nodes[ng,2])**2 ), math.sqrt( (self.Lx-self.nodes[ng,1])**2+ (self.Lx-self.nodes[ng,2])**2 )) 
                             
                            r2=max(rnode,rng)
                            r1=min(rnode,rng)
                            res = np.log(r2/r1)/(2*np.pi*k_diph)
                        else :
                            dxx=abs(self.nodes[ng,2] - self.nodes[idnode,2])
                            dyy=abs(self.nodes[ng,1] - self.nodes[idnode,1])		
                            l=np.sqrt(dxx**2 + dyy**2)
                            res= dxx / (self.k_liq * l)
                        self.R[idnode,j]= res
                    j+=1




    def capacite_cart(self):
        dx=self.nodes[1,2] - self.nodes[0,2]
        dy=self.nodes[self.Nptsx,1] - self.nodes[0,1]
        for idnode in range(self.Nptsy*self.Nptsx) :
            rho_diph=self.phi[idnode]*self.rho_gas+(1-self.phi[idnode])*self.rho_liq
            cp_diph=self.phi[idnode]*self.cp_gas+(1-self.phi[idnode])*self.cp_liq
            self.C[idnode] = rho_diph*cp_diph*dx*dy

    def capacite_tank(self):
        dx_cart=self.nodes[1,2] - self.nodes[0,2]
        dy_cart=self.nodes[self.Nptsx,1] - self.nodes[0,1]
        for idnode in range(0,self.Nptsy*self.Nptsx+2*(self.Nptsx-1)*self.ntheta):
            rho_diph=self.phi[idnode]*self.rho_gas+(1-self.phi[idnode])*self.rho_liq
            cp_diph=self.phi[idnode]*self.cp_gas+(1-self.phi[idnode])*self.cp_liq
            #middle of the 'rectangle part'
            if (self.nodes[idnode,1]>self.Lx and self.nodes[idnode,1]<self.Ly-self.Lx):
                self.C[idnode] = rho_diph*cp_diph*dx_cart*dy_cart
            #lower boundary of the 'rectangle part'
            elif (idnode<self.Nptsx-1):
                self.C[idnode] = rho_diph*cp_diph*(dx_cart*dy_cart/2+self.angle/2*(2*dx_cart*math.sqrt( (self.Lx-self.nodes[idnode,1])**2+ (self.Lx-self.nodes[idnode,2])**2 ))/2)
            #lower right corner of the 'rectangle part'
            elif (idnode==self.Nptsx-1):
                self.C[idnode]=rho_diph*cp_diph*(dx_cart*dy_cart/2+math.pi/4*(dx_cart/2)**2)
            #lower 'circle part'
            elif (self.nodes[idnode,1]>self.Ly-self.Lx):
                self.C[idnode] = rho_diph*cp_diph*(self.angle/2*(2*dx_cart*math.sqrt( (self.Lx-self.nodes[idnode,1])**2+ (self.Lx-self.nodes[idnode,2])**2 )))
            #upper boundary of the 'rectangle part' 
            elif (idnode>self.Nptsx*(self.Nptsy-1)-1 and idnode<self.Nptsx*self.Nptsy-1):
                self.C[idnode] = rho_diph*cp_diph*(dx_cart*dy_cart/2+self.angle/2*(2*dx_cart*math.sqrt( (self.Ly-self.Lx-self.nodes[idnode,1])**2+ (self.Lx-self.nodes[idnode,2])**2 ))/2)
            #upper right corner of the 'rectangle part'
            elif (idnode==self.Nptsx*self.Nptsy-1):
                self.C[idnode]=rho_diph*cp_diph*(dx_cart*dy_cart/2+math.pi/4*(dx_cart/2)**2)
            #upper 'circle part'
            elif (self.nodes[idnode,1]<self.Lx):
                self.C[idnode] = rho_diph*cp_diph*(self.angle/2*(2*dx_cart*math.sqrt( (self.Ly-self.Lx-self.nodes[idnode,1])**2+ (self.Lx-self.nodes[idnode,2])**2 )))





#test=Init()
#test.domain_cart()
#test.init_domain()
#test.capacite_tank()
#test.resistance_tank()
#print(self.R)
