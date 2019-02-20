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
		self.temp2d=np.zeros((self.Nptsx,self.Nptsy))
		self.meshx, self.meshy = np.meshgrid(self.x, self.y)
		self.domain=np.zeros((self.Nptsx*self.Nptsy+2*(self.Nptsx-1)*self.ntheta))
	def systeme_init_cart(self, temp) :
		self.domain_cart()
		self.init_domain()
		self.domain[:,0]=self.nodes[:,2]
		self.domain[:,1]=self.nodes[:,1]
		self.resistance_cart()
	#	self.initemp_cart_y()
		self.initemp_cart_x()
		self.initemp_cart_y()
		self.capacite_cart()
		temp[:]=np.copy(self.temp)
		#~ print()
		#~ Rajout condition initiale

		

		
	def systeme_cond(self, T, dT_dt, time=0.0):
		taille=self.Nptsy*self.Nptsx
		for idnode in range(taille) :
			j=1
			dT_dt[idnode]=0
			C=1./self.C[idnode]
			while (j<5 and (int(self.neig[idnode,j]) != -1)):
				G=1./self.R[idnode,j]
				ng=int(self.neig[idnode,j])
				deltaT=T[ng] - T[idnode]
				dT_dt[idnode]+= G*deltaT				
				j+=1
			dT_dt[idnode]=dT_dt[idnode] * C

	def champs_reconstruct(self):
		pass


