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

from init import Init

class Reservoir(Init) :

	def __init__(self):
		Init.__init__(self)
		
	def systeme_init(self, Tinit) :
		T0=0.
		Tinit = np.zeros(self.Nptsy*Nptsx)
		Tinit=T0
		print
	
	def systeme(self,T, dT_dt ):
		pass
		
