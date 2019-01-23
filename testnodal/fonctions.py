#!/usr/bin/python

import numpy as np
import math as math

def solveur1D_df(mesh,x,dt,alpha,Nptsx,Lx,residu,T1,T0):
    """Solves the 1D temperture equation using finite difference method """
    t=dt
    temp=np.copy(mesh)
    temp_tmp=np.copy(mesh)
    residu_i=1
    
    #convergence loop
    while (residu < residu_i):
        #print(residu_i)
        for i in range(1,Nptsx-1):
            #print(i)
            temp_tmp[:,i]=temp[:,i]-alpha*dt*(2*temp[:,i]-temp[:,i-1]-temp[:,i+1])/((Lx*Lx)/(Nptsx*Nptsx))
        
        temp_tmp[:,Nptsx-1]=T0-(T0-T1)*math.erf(0.5*x[Nptsx-1]/math.sqrt(alpha*t))
    
        residu_i=max(abs(temp_tmp[0,:]-temp[0,:]))
        temp=np.copy(temp_tmp)
        t=t+dt
        #print(residu_i)
    return temp


def solveur1D_md(mesh,x,dx,dt,alpha,Nptsx,Lx,residu,T1,T0):
    """Solves the 1D temperture equation using the nodal approach """
    t=dt
    temp=np.copy(mesh)
    temp_tmp=np.copy(mesh)
    residu_i=1

    #convergence loop
    while (t<500*dt):#(residu < residu_i):
        
        for i in range(1,Nptsx-1): #find the range
            
            temp_tmp[:,i]=temp[:,i]-alpha*(temp[:,i+1]-temp[:,i-1])*dt/(dx**2) #find the right model
        
        temp_tmp[:,Nptsx-1]=T0-(T0-T1)*math.erf(0.5*x[Nptsx-1]/math.sqrt(alpha*t))
 
        residu_i=max(abs(temp_tmp[0,:]-temp[0,:]))
        temp=np.copy(temp_tmp)
        t=t+dt
        
    return temp
