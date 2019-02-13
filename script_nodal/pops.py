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

	def systeme_init_cart(self, temp) :
		self.domain_cart()
		self.init_domain()
		self.resistance_cart()
		self.initemp_cart()
		temp[:]=np.copy(self.temp)
		#~ print()
		#~ Rajout condition initiale

		

		
	def systeme_cond(self, T, dT_dt, time=0.0):
		taille=self.Nptsy*self.Nptsx
		for idnode in range(taille) :
			j=1
			dT_dt[idnode]=0
			while (j<5 and (int(self.neig[idnode,j]) != -1)):
				C=1. #inverse capacite a implementer
				G=1./self.R[idnode,j]
				ng=int(self.neig[idnode,j])
				deltaT=T[ng] - T[idnode]
				dT_dt[idnode]+= G*deltaT *C
				j+=1


		#~ INITIALISATION 

Reservoir1=Reservoir()
#~ Reservoir1.systeme_init_cart()


		#~ INTEGRATION TEMPORELLE ET ECRITURE
		
ProblemSize=Reservoir1.Nptsx*Reservoir1.Nptsy
listVar=[i for i in range(ProblemSize)]
# ~ listVar=[0,1]
AdapTimeStepBool=False
Duration=1.
MAXNTIMESTEP=1001 
TIMESTEP=0.002
METHODE='Euler'
ListOfVariablesToSave=listVar
SavedIteration=100
Error=1e-8

Problem=PbDef.NumericalProblem(Reservoir1.systeme_init_cart,Reservoir1.systeme_cond,METHODE,\
MAXNTIMESTEP,ProblemSize,TimeStep=TIMESTEP,\
Duration=Duration,AdaptativeTimeStep=AdapTimeStepBool,\
NbIterationSaved=SavedIteration,\
ListOfIdVar=ListOfVariablesToSave,AdatativeTimeStep_Error=Error)

Problem.SolveDifferentialProblem()

