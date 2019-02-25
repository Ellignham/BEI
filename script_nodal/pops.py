#!/usr/bin/python2.7

#Imports
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math as math
import sys
LibPath='../../Sujet_29/'
sys.path.append(LibPath)
import NumericalProblemClass as PbDef

#Class
from input import Input
from init import Init
# ~ from visu import lecture_champs
# ~ from visu import *



class Reservoir(Init) :

    def __init__(self):
        Init.__init__(self)
        self.meshx, self.meshy = np.meshgrid(self.x, self.y)
        self.domain=np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta))
		
		
	
    def systeme_init(self, temp) :
        if (self.mesh_type=='cart'):
            self.init_phase_cart()
            self.domain_cart()
            self.init_domain_cart()
            #~ self.domain[:,0]=self.nodes[:,2]
            #~ self.domain[:,1]=self.nodes[:,1]
            self.resistance_cart()
            #self.initemp_cart_x()
            self.initemp_cart_y()
            self.capacite_cart()
            temp[:]=np.copy(self.temp)
        elif (self.mesh_type=='tank'):
            self.init_phase_tank()
            self.domain_tank()
            self.init_domain_tank()
            #~ self.domain[:,0]=self.nodes[:,2]
            #~ self.domain[:,1]=self.nodes[:,1]
            self.resistance_tank()
            self.initemp_tank_x()
            #~ self.initemp_tank_y()
            self.capacite_tank()
            temp[:]=np.copy(self.temp)

    def systeme_cond(self, T, dT_dt, time=0.0):
		taille=len(T)
		for idnode in range(taille) :
			j=1
			dT_dt[idnode]=0
			#~ C=1./self.C[idnode]
			C=1.
			#~ print(int(self.neig[idnode,j]))
			while ((int(self.neig[idnode,j]) != -1)):
				#~ print('toto')
				if self.R[idnode,j] !=0. :
					G=1./self.R[idnode,j]
					ng=int(self.neig[idnode,j])
					deltaT=T[ng] - T[idnode]
					dT_dt[idnode]+= G*deltaT				
				j+=1
			dT_dt[idnode]=dT_dt[idnode] * C


    def interface(self,time):
        interface=np.loadtxt("LOX_Height_vs_Time_BEI.txt")
        #print(interface)

        loop=True
        i=0

        while loop :
            if interface[i,0]>time:
                loop=False
            else :
                i+=1            
        
        if interface[i-1,0]==time:
            self.height=interface[i-1,1]
        else :
            t1=interface[i-1,0]
            t2=interface[i,0]

            self.height=(t2-time)/(t2-t1)*interface[i-1,1]+(time-t1)/(t2-t1)*interface[i,1]
    
        print(self.height)
   
    def update_constant(self,T):
        pass

test=Reservoir()
test.interface(900)
