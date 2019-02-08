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




class Reservoir(Init) :

	def __init__(self):
		Init.__init__(self)
		self.Tinit=0.

	def systeme_init(self) :
		T0=10.
		self.Tinit = T0*np.ones(self.Nptsy*self.Nptsx)
		print('Tinit dasn init', self.Tinit)
		
		
		
	def systeme(self,T, dT_dt, time=0.0):
		pass
		
Reservoir1=Reservoir()
print(Reservoir1.Nptsx)
print(Reservoir1.Nptsy)
Reservoir1.systeme_init()
print('Tinit=', Reservoir1.Tinit)
