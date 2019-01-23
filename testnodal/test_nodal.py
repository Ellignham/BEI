#!/usr/bin/python


"""Script to test the implementation of the nodal approach and compare it to finite diference method """

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from fonctions import *

##Input variables##

#Definition of the domain
Lx=1.
Ly=1.

#Resolution
Nptsx=20
Nptsy=2

##Definition of the constants needed for the calculation##
T0=293
T1=313
#Water at 20C
alpha=0.143*10**(-5)
conductance=598.03*10**(-3)
resistance=1/conductance
#Computation constants
dt=0.4*((Lx**2)/(Nptsx**2))/(2*alpha)
residu=10**(-5)



##Initialisation##

#create arrays
mesh_df=np.zeros((Nptsy,Nptsx))
mesh_md=np.copy(mesh_df)
x=np.linspace(0,Lx,Nptsx)
y=np.linspace(0,Ly,Nptsy)
dx=x[1]-x[0]

#Initial values of arrays
mesh_df[:,:]=T0
mesh_df[:,0]=T1

mesh_md[:,:]=T0
mesh_md[:,0]=T1

##Solver##

#reference
temp_df=solveur1D_df(mesh_df,x,dt,alpha,Nptsx,Lx,residu,T1,T0)
temp_md=solveur1D_md(mesh_md,x,dx,dt,alpha,Nptsx,Lx,residu,T1,T0)



##Plot##

#2Dplot of the finite element method
plt.figure()

plt.imshow(temp_df.reshape((Nptsy,Nptsx)),extent=(x.min(), x.max(), y.max(), y.min()), 
    interpolation='nearest',cmap=cm.gray)

plt.colorbar()
plt.grid()
plt.show(block=False)

#2Dplot of the nodal method
plt.figure()

plt.imshow(temp_md.reshape((Nptsy,Nptsx)),extent=(x.min(), x.max(), y.max(), y.min()), 
    interpolation='nearest',cmap=cm.gray)

plt.colorbar()
plt.grid()
plt.show(block=False)


#1D plot to compare the 2 methods
plt.figure()
plt.plot(x,temp_df[0,:])
plt.plot(x,temp_md[0,:])
plt.xlabel('x')
plt.ylabel('Temp')
plt.show()



