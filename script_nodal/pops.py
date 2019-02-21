#!/usr/bin/python3

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
		#~ self.domain_cart()
		self.domain_tank()
		self.init_domain()
		#~ self.domain[:,0]=self.nodes[:,2]
		#~ self.domain[:,1]=self.nodes[:,1]
		#~ self.resistance_cart()
		self.resistance_tank()
	#	self.initemp_cart_y()
		self.initemp_cart_x()
		self.initemp_cart_y()
		self.capacite_tank()
		temp[:]=np.copy(self.temp)


		

		
	def systeme_cond(self, T, dT_dt, time=0.0):
		taille=self.Nptsy*self.Nptsx
		for idnode in range(taille) :
			j=1
			dT_dt[idnode]=0
			C=1./self.C[idnode]
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



