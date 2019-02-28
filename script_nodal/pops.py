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
            self.init_phase()
            self.domain_cart()
            self.init_domain()
            #~ self.domain[:,0]=self.nodes[:,2]
            #~ self.domain[:,1]=self.nodes[:,1]
            self.resistance_cart()
            #self.initemp_cart_x()
            self.initemp_y()
            self.capacite_cart()
            temp[:]=np.copy(self.temp)
        elif (self.mesh_type=='tank'):
            #~ Generation du maillage
            self.domain_tank()
            #~ Initialisation de la fonction de phase (uniforme a 0 => liq)
            self.init_phase()
            #~ Trouve la hauteur initiale de l'interface
            self.hauteur_interface(self.time_init)
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
                
            
            
            

    def systeme_cond(self, T, dT_dt):
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
            
            
            
            
            
            
            
            
    def systeme_diph(self, T, dT_dt,time=0.0):
        self.hauteur_interface(time+self.time_init)
        self.width_interface()
        phi_old = np.copy(self.phi)
        #~ update tableau des phi
        self.update_phi()
        #~ calcul des flux  (isentropique)
        flux_pc = -self.Hlv * (self.phi - phi_old)
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
                    #~ flux provenant des autres noeuds par conduction 
                    G=1./self.R[idnode,j]
                    ng=int(self.neig[idnode,j])
                    deltaT=T[ng] - T[idnode]
                    dT_dt[idnode]+= G*deltaT				
                j+=1
            #~ bilan des flux
            dT_dt[idnode]=dT_dt[idnode] * C + flux_pc[idnode]
            

    def hauteur_interface(self,time):
        '''
        Computes the position of the interface
        '''
        print(time)
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
            self.height=interface[i-1,1]*self.Ly
        else :
            t1=interface[i-1,0]
            t2=interface[i,0]

            self.height=(t2-time)/(t2-t1)*interface[i-1,1]+(time-t1)/(t2-t1)*interface[i,1]*self.Ly
   
    def width_interface(self):
        #Points around the interface
        pts=[]
        for idnode in range(0,len(self.temp)):
            if abs(self.nodes[idnode,1]-self.height)<=self.dz :
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
        for k in range(0,len(gradTL)): 
            moyL+=gradTL[k]
        moyL=moyL/len(gradTL)

        moyG=0
        for k in range(0,len(gradTG)):
            moyG+=gradTG[k]
        moyG=moyG/len(gradTG)
	
        self.dz=abs((1/self.Hlv)*(self.k_liq*moyL-self.k_gas*moyG))
        
        
    def update_phi(self):
        
        #~ Creation du domaine de changement de phase 
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

test=Reservoir()

if (test.mesh_type=='tank'):		
    ProblemSize=test.Nptsx*test.Nptsy+2*(test.Nptsx-1)*test.ntheta
elif (test.mesh_type=='cart'):
    ProblemSize=test.Nptsx*test.Nptsy

prout = np.zeros(ProblemSize)
test.systeme_init(prout)
#~ test.update_phi()

#~ print(test.temp)
#~ 

exec(open("./visu.py").read())
temps, temperature, x, y = reconstruct_champs(test, 'ResultArray.dat')
plot_temp_int(test, x, y, test.temp, temps[0])
