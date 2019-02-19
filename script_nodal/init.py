#!/usr/bin/python3

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
        neig        : |id of node|neighbor 1|neighbor 2| ...  if neighbor = -1 : boundary

        """
        
        Input.__init__(self)
        
        self.x = np.linspace(0,self.Lx/2,self.Nptsx)
        self.y = np.linspace(self.Lx/2,self.Ly-self.Lx/2,self.Nptsy)

        if (self.mesh_type=='cart'):
            self.nodes = np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta,3))
            self.neig = np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta,6+self.ntheta))
        else:
            self.nodes = np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta,3))
            self.neig = np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta,4+self.ntheta))


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

        self.neig[0,1]=1; self.neig[0,2]=self.Nptsx;self.neig[0,3]=-1                            
        self.neig[self.Nptsx-1,1]=self.Nptsx-2; self.neig[self.Nptsx-1,2]=2*self.Nptsx-1;self.neig[self.Nptsx-1,3]=-1
      
        self.neig[self.Nptsx*(self.Nptsy-1),1]=self.Nptsx*(self.Nptsy-1)+1; self.neig[self.Nptsx*(self.Nptsy-1),2]=self.Nptsx*(self.Nptsy-2);self.neig[self.Nptsx*(self.Nptsy-1),3]=-1
        self.neig[self.Nptsx*self.Nptsy-1,1]=self.Nptsx*self.Nptsy-2; self.neig[self.Nptsx*self.Nptsy-1,2]=self.Nptsx*(self.Nptsy-1)-1;self.neig[self.Nptsx*self.Nptsy-1,3]=-1
    
        
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
                    self.neig[i+j*self.Nptsx,3]=i+(j+1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,4]=-1
                #upper boundary
                elif (j==self.Nptsy-1 and i>0 and i<self.Nptsx-1):
                    self.neig[i+j*self.Nptsx,1]=i-1+j*self.Nptsx 
                    self.neig[i+j*self.Nptsx,2]=i+1+j*self.Nptsx
                    self.neig[i+j*self.Nptsx,3]=i+(j-1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,4]=-1
                #left boundary
                elif (i==0 and j>0 and j<self.Nptsy-1):
                    self.neig[i+j*self.Nptsx,1]=i+1+j*self.Nptsx 
                    self.neig[i+j*self.Nptsx,2]=i+(j-1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,3]=i+(j+1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,4]=-1
                #right boundary
                elif (i==self.Nptsx-1 and j>0 and j<self.Nptsy-1):
                    self.neig[i+j*self.Nptsx,1]=i-1+j*self.Nptsx 
                    self.neig[i+j*self.Nptsx,2]=i+(j+1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,3]=i+(j-1)*self.Nptsx
                    self.neig[i+j*self.Nptsx,4]=-1
        #corners
        self.neig[0,1]=1; self.neig[0,2]=self.Nptsx;self.neig[0,3]=-1                            
        self.neig[self.Nptsx-1,1]=self.Nptsx-2; self.neig[self.Nptsx-1,2]=2*self.Nptsx-1;self.neig[self.Nptsx-1,3]=-1
      
        self.neig[self.Nptsx*(self.Nptsy-1),1]=self.Nptsx*(self.Nptsy-1)+1; self.neig[self.Nptsx*(self.Nptsy-1),2]=self.Nptsx*(self.Nptsy-2);self.neig[self.Nptsx*(self.Nptsy-1),3]=-1
        self.neig[self.Nptsx*self.Nptsy-1,1]=self.Nptsx*self.Nptsy-2; self.neig[self.Nptsx*self.Nptsy-1,2]=self.Nptsx*(self.Nptsy-1)-1;self.neig[self.Nptsx*self.Nptsy-1,3]=-1
   
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
                    self.neig[id_node,1]=id_node-1
                    self.neig[id_node,2]=id_node+1
                    self.neig[id_node,3]=id_node+self.Nptsx-1
                    self.neig[id_node,4]=id_node-self.Nptsx+1
                    self.neig[id_node,5]=-1
                #border of the circle
                elif (r==self.Nptsx-1 and theta>1 and theta<self.ntheta):
                    self.neig[id_node,1]=id_node-1
                    self.neig[id_node,2]=id_node+self.Nptsx-1
                    self.neig[id_node,3]=id_node-self.Nptsx+1
                    self.neig[id_node,4]=-1
                #right boundary
                elif (r>1 and r<self.Nptsx-1 and theta==self.ntheta):
                    self.neig[id_node,1]=id_node-1
                    self.neig[id_node,2]=id_node+1
                    self.neig[id_node,3]=id_node+self.Nptsx-1
                    self.neig[id_node,4]=-1
                #upper boundary
                elif (r>1 and r<self.Nptsx-1 and theta==1 ):
                    self.neig[id_node,1]=id_node-1
                    self.neig[id_node,2]=id_node+1
                    self.neig[id_node,3]=self.Nptsx-1-r
                    self.neig[id_node,4]=id_node+self.Nptsx-1
                    self.neig[self.Nptsx-1-r,4]=id_node
                    self.neig[self.Nptsx-1-r,5]=-1
                #center of the circle
                elif (r==1 and theta>1 and theta<self.ntheta):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=id_node+self.Nptsx-1
                    self.neig[id_node,3]=id_node-self.Nptsx+1
                    self.neig[id_node,4]=self.Nptsx-1
                    self.neig[id_node,5]=-1
                    self.neig[self.Nptsx-1,2+theta]=id_node
                #corners
                elif (r==self.Nptsx-1 and theta==1):
                    self.neig[id_node,1]=id_node-1
                    self.neig[id_node,2]=id_node+self.Nptsx-1
                    self.neig[id_node,3]=0
                    self.neig[id_node,4]=-1
                    self.neig[0,3]=id_node
                    self.neig[0,4]=-1
                elif (r==self.Nptsx-1 and theta==self.ntheta):
                    self.neig[id_node,1]=id_node-1
                    self.neig[id_node,2]=id_node-self.Nptsx+1
                    self.neig[id_node,3]=-1
                elif (r==1 and theta==1):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=id_node+self.Nptsx-1
                    self.neig[id_node,3]=self.Nptsx-2
                    self.neig[id_node,4]=self.Nptsx-1
                    self.neig[id_node,5]=-1
                    self.neig[self.Nptsx-2,4]=id_node
                    self.neig[self.Nptsx-2,5]=-1
                    self.neig[self.Nptsx-1,3]=id_node
                elif (r==1 and theta==self.ntheta):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=id_node-self.Nptsx+1
                    self.neig[id_node,3]=self.Nptsx-1
                    self.neig[id_node,4]=-1
                    self.neig[self.Nptsx-1,2+theta]=id_node
                    self.neig[self.Nptsx-1,3+theta]=-1
        
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
                    self.neig[id_node,1]=id_node-1
                    self.neig[id_node,2]=id_node+1
                    self.neig[id_node,3]=id_node+self.Nptsx-1
                    self.neig[id_node,4]=id_node-self.Nptsx+1
                    self.neig[id_node,5]=-1
                #border of the circle
                elif (r==self.Nptsx-1 and theta>1 and theta<self.ntheta):
                    self.neig[id_node,1]=id_node-1
                    self.neig[id_node,2]=id_node+self.Nptsx-1
                    self.neig[id_node,3]=id_node-self.Nptsx+1
                    self.neig[id_node,4]=-1
                #right boundary
                elif (r>1 and r<self.Nptsx-1 and theta==self.ntheta):
                    self.neig[id_node,1]=id_node-1
                    self.neig[id_node,2]=id_node+1
                    self.neig[id_node,3]=id_node-self.Nptsx+1
                    self.neig[id_node,4]=-1
                #lower boundary
                elif (r>1 and r<self.Nptsx-1 and theta==1 ):
                    self.neig[id_node,1]=id_node-1
                    self.neig[id_node,2]=id_node+1
                    self.neig[id_node,3]=self.Nptsy*self.Nptsx-1-r
                    self.neig[id_node,4]=id_node+self.Nptsx-1
                    self.neig[id_node,5]=-1
                    self.neig[self.Nptsy*self.Nptsx-1-r,4]=id_node
                    self.neig[self.Nptsy*self.Nptsx-1-r,5]=-1
                #center of the circle
                elif (r==1 and theta>1 and theta<self.ntheta):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=id_node+self.Nptsx-1
                    self.neig[id_node,3]=id_node-self.Nptsx+1
                    self.neig[id_node,4]=self.Nptsx*self.Nptsy-1
                    self.neig[id_node,5]=-1
                    self.neig[self.Nptsy*self.Nptsx-1,2+theta]=id_node
                #corners
                elif (r==self.Nptsx-1 and theta==1):
                    self.neig[id_node,1]=id_node-1
                    self.neig[id_node,2]=id_node+self.Nptsx-1
                    self.neig[id_node,3]=self.Nptsx*(self.Nptsy-1)
                    self.neig[id_node,4]=-1
                    self.neig[self.Nptsx*(self.Nptsy-1),3]=id_node
                    self.neig[self.Nptsx*(self.Nptsy-1),4]=-1
                elif (r==self.Nptsx-1 and theta==self.ntheta):
                    self.neig[id_node,1]=id_node-1
                    self.neig[id_node,2]=id_node-self.Nptsx+1
                    self.neig[id_node,3]=-1
                elif (r==1 and theta==1):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=id_node+self.Nptsx-1
                    self.neig[id_node,3]=self.Nptsy*self.Nptsx-2
                    self.neig[id_node,4]=self.Nptsy*self.Nptsx-1
                    self.neig[id_node,5]=-1
                    self.neig[self.Nptsy*self.Nptsx-2,4]=id_node
                    self.neig[self.Nptsy*self.Nptsx-2,5]=-1
                    self.neig[self.Nptsy*self.Nptsx-1,3]=id_node
                elif (r==1 and theta==self.ntheta):
                    self.neig[id_node,1]=id_node+1
                    self.neig[id_node,2]=id_node-self.Nptsx+1
                    self.neig[id_node,3]=self.Nptsy*self.Nptsx-1
                    self.neig[id_node,4]=-1
                    self.neig[self.Nptsy*self.Nptsx-1,2+theta]=id_node
                    self.neig[self.Nptsy*self.Nptsx-1,3+theta]=-1

  #      print(self.nodes)
  #      print(self.neig)
 
        fig=plt.figure()
        ax = fig.add_subplot(111)
        plt.plot(self.nodes[:,2],self.nodes[:,1],'or')
        for i, txt in enumerate(self.nodes[:,0]):
            ax.annotate(txt, (self.nodes[i,2], self.nodes[i,1]))
        plt.xlim(0,0.5)
        plt.ylim(0,0.5)
      #  plt.ylim(8,10.5)
      #  plt.xlim(0.5,1.4)
#        plt.show()

    def init_domain(self):
        """
        Initialise the domain before the computation


        Variables :

        temp        : array containing the temperature field [K]
        pres        : array containing the pressure field [Pa]

        U           : array containing the x velocity field [m/s]
        V           : array containing the y velocity field [m/s]

        R           : array containing the conduction thermal resistance [K.m/W]
        C           : array containing the conduction thermal capacity

        phi         : array describing whether the fluid is liquid or gas 
        """
        
        #Creation of the temperature array
        self.temp = np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta))
        
        #Creation of the pressure array
        self.pres = np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta)) 
        self.pres[:] = self.pfluid_init 
        
        #Creation of the velocity fields
        self.U = np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta))
        self.V = np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta))
        self.U[:]=self.ufluid_init
        self.V[:]=self.vfluid_init

        #Creation of thermal resistance arrays
        self.Rx = np.zeros(self.Nptsy*(self.Nptsx-1))
        self.Ry = np.zeros((self.Nptsy-1)*self.Nptsx)
        self.R = np.zeros((self.Nptsy*self.Nptsx, 5))#.reshape((len(self.nodes),1))
        
        #Creation of thermal capacity array
        self.C = np.zeros((self.Nptsy*self.Nptsx+2*(self.Nptsx-1)*self.ntheta))
   
        #Creation of the phase array
        self.phi = np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta))
 
    def initemp_cart_y(self):
        for k in range(0,self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta):
            if self.nodes[k,1]<self.Ly/2 :
                self.temp[k]=self.T1
            else :
                self.temp[k]=self.tfluid_init
    
    def initemp_cart_x(self):
        for k in range(0,self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta):
            if self.nodes[k,2]<self.Lx/4:
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

    def capacite_cart(self):
        dx=self.nodes[1,2] - self.nodes[0,2]
        dy=self.nodes[self.Nptsx,1] - self.nodes[0,1]
        for idnode in range(self.Nptsy*self.Nptsx) :
            self.C[idnode] = self.rho*self.cp*dx*dy

    def capacite_tank(self):
        dx_cart=self.nodes[1,2] - self.nodes[0,2]
        dy_cart=self.nodes[self.Nptsx,1] - self.nodes[0,1]
        for idnode in range(0,self.Nptsy*self.Nptsx+2*(self.Nptsx-1)*self.ntheta):
            #middle of the 'rectangle part'
            if (self.nodes[idnode,1]>self.Lx and self.nodes[idnode,1]<self.Ly-self.Lx):
                self.C[idnode] = self.rho*self.cp*dx_cart*dy_cart
            #lower boundary of the 'rectangle part'
            elif (idnode<self.Nptsx-1):
                self.C[idnode] = self.rho*self.cp*(dx_cart*dy_cart/2+self.angle/2*(2*dx_cart*math.sqrt( (self.Lx-self.nodes[idnode,1])**2+ (self.Lx-self.nodes[idnode,2])**2 ))/2)
            #lower right corner of the 'rectangle part'
            elif (idnode==self.Nptsx-1):
                self.C[idnode]=self.rho*self.cp*(dx_cart*dy_cart/2+math.pi/4*(dx_cart/2)**2)
            #lower 'circle part'
            elif (self.nodes[idnode,1]>self.Ly-self.Lx):
                self.C[idnode] = self.rho*self.cp*(self.angle/2*(2*dx_cart*math.sqrt( (self.Lx-self.nodes[idnode,1])**2+ (self.Lx-self.nodes[idnode,2])**2 )))
            #upper boundary of the 'rectangle part' WIP
            elif (idnode>self.Nptsx*(self.Nptsy-1)-1 and idnode<self.Nptsx*self.Nptsy-1):
                self.C[idnode] = self.rho*self.cp*(dx_cart*dy_cart/2+self.angle/2*(2*dx_cart*math.sqrt( (self.Ly-self.Lx-self.nodes[idnode,1])**2+ (self.Lx-self.nodes[idnode,2])**2 ))/2)
            #upper right corner of the 'rectangle part'
            elif (idnode==self.Nptsx*self.Nptsy-1):
                self.C[idnode]=self.rho*self.cp*(dx_cart*dy_cart/2+math.pi/4*(dx_cart/2)**2)
            #upper 'circle part'
            elif (self.nodes[idnode,1]<self.Lx):
                self.C[idnode] = self.rho*self.cp*(self.angle/2*(2*dx_cart*math.sqrt( (self.Ly-self.Lx-self.nodes[idnode,1])**2+ (self.Lx-self.nodes[idnode,2])**2 )))




#test=Init()
#test.domain_tank()
#test.init_domain()
#test.capacite_tank()
