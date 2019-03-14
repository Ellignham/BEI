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
		'''
		Calls all functions required to initialise either the cartesian domain or the tank
		'''
        if (self.mesh_type=='cart'):
            #~ Generation du maillage
            self.domain_cart()
            #~ Initialisation de la fonction de phase (uniforme a 0 => liq)
            self.init_phase()
            #~ Trouve la hauteur initiale de l'interface
            self.hauteur_interface(0.0)
            #~ Calcul de phi initial avec interface
            self.update_phi()
            #~ Cree les vecteurs de temp, pres, vitesse, resistance et capa
            #~ Initialise pres, vitesse
            self.init_domain()
            #~ initialise le champs de temp avec l'interface
            self.initemp()
			#~ Initialise les resistance
            self.resistance_cart()
            #~ Initialise les capacites des noeuds dans le tank
            self.capacite_cart()
            temp[:]=np.copy(self.temp)
                
        elif (self.mesh_type=='tank'):
            #~ Generation du maillage
            self.domain_tank()
            #~ Initialisation de la fonction de phase (uniforme a 0 => liq)
            self.init_phase()
            #~ Trouve la hauteur initiale de l'interface
            self.hauteur_interface(0.0)
            #~ Calcul de phi initial avec interface
            self.update_phi()
            #~ Cree les vecteurs de temp, pres, vitesse, resistance et capa
            #~ Initialise pres, vitesse
            self.init_domain()
            #~ initialise le champs de temp avec l'interface
            self.initemp()
            #~ Initialise les resistance
            self.resistance_tank()
            #~ Initialise les capacites des noeuds dans le tank
            self.capacite_tank()
            temp[:]=np.copy(self.temp)
                
            
            
            

    def systeme_cond(self, T, dT_dt, time=0.0):
        taille=len(T)
        for idnode in range(taille) :
            j=1
            dT_dt[idnode]=0
            C=1./self.C[idnode]
            #~ C=1.
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
            
            
            
            
    def systeme_diph(self, T, dT_dt, time=0.0):
        
        
        self.hauteur_interface(time)
        self.width_interface()
        phi_old = np.copy(self.phi)
        #~ update tableau des phi
        self.update_phi()
        # ~ Update des resistances et capa
        self.resistance_tank()
        self.capacite_tank()
        #~ calcul des flux 
        flux_pc = np.zeros(self.dom_size)
        #~ print('time', time)
        #~ print(self.height)

        taille=len(T)
        for idnode in range(taille) :
            j=1
            dT_dt[idnode]=0
            C=1./self.C[idnode]
            flux_pc[idnode] = -self.rho_diph[idnode] * self.Hlv * (self.phi[idnode] - phi_old[idnode])/self.dt
            while ((int(self.neig[idnode,j]) != -1)):
                if self.R[idnode,j] !=0. :
                    #~ flux provenant des autres noeuds par conduction 
                    G=1./self.R[idnode,j]
                    ng=int(self.neig[idnode,j])
                    deltaT=T[ng] - T[idnode]
                    dT_dt[idnode]+= G*deltaT				
                j+=1
            #~ bilan des flux
            dT_dt[idnode]=C * (dT_dt[idnode] + flux_pc[idnode])
          

    def hauteur_interface(self,time):
        '''
        Computes the position of the interface
        '''
        time=time+self.time_init 
        interface=np.loadtxt("LOX_Height_vs_Time_BEI.txt")
   
        loop=True
        i=0

        while loop :
            if interface[i,0]>time:
                loop=False
            else :
                i+=1            
        
        if interface[i-1,0]==time:
            self.height=interface[i-1,1]*self.Ly
        else :
            t1=interface[i-1,0]
            t2=interface[i,0]

            self.height=((t2-time)/(t2-t1)*interface[i-1,1]+(time-t1)/(t2-t1)*interface[i,1])*self.Ly
   
    def width_interface(self):
		'''
		Computes the width of the interface
		'''
        #Points around the interface
        pts=[]
        epaisseur=self.dy
        for idnode in range(0,self.dom_size):
            if abs(self.nodes[idnode,1]-self.height)<=epaisseur :
                pts.append(idnode)
        
        #Width of the interface
        gradTL=[]
        gradTG=[]
 
        for k in range(0,len(pts)):
            #below the interface
            if self.nodes[pts[k],1]<self.height :
                grad=(self.Tint-self.temp[pts[k]])/(self.height-self.nodes[pts[k],1])
                gradTL.append(grad)
            #above the interface
            else :
                grad=(self.Tint-self.temp[pts[k]])/(self.height-self.nodes[pts[k],1])
                gradTG.append(grad)
        moyL=0
        if len(gradTL)!=0:
            for k in range(0,len(gradTL)): 
                moyL+=gradTL[k]
            moyL=moyL/len(gradTL)

        moyG=0
        if len(gradTG)!=0 :
            for k in range(0,len(gradTG)):
                moyG+=gradTG[k]
            moyG=moyG/len(gradTG)

        sect=self.Lx
        	
        self.dz=abs((1/self.Hlv*sect*self.rho_liq)*(self.k_liq*moyL-self.k_gas*moyG))
        
        
    def update_phi(self):
		'''
		Updates the void fraction function
		'''
        
        #~ Creation du domaine de changement de phase avec interface epaisse
        yimin=self.height-self.dz
        yimax=self.height+self.dz
        a= 1. / (yimax - yimin)
        
        #~ Calcul dans domaine pc
        imin=self.nodes[:,1]>yimin
        imax=self.nodes[:,1]<yimax
        it=np.where(imin*imax)
        self.phi[it] = a * ( self.nodes[it,1] - yimin)

        #~ Domaine liquide (bas)
        imin=self.nodes[:,1]<yimin
        it=np.where(imin)
        self.phi[it] = 0.

        #~ Domaine gaz (haut)
        imax=self.nodes[:,1]>yimax
        it=np.where(imax)
        self.phi[it] = 1.
