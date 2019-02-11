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

	def systeme_init_cart(self, temp) :
		self.domain_cart()
		self.init_domain()
		self.resistance_cart()
		self.initemp_cart()
		temp=self.temp
		#~ Rajout condition initiale

		

		
	def systeme_cond(self, T, dT_dt, time=0.0):
		taille=self.Nptsy*self.Nptsx
		dT_dt=np.zeros(taille)
		for idnode in range(taille) :
			j=1
			dT_dt[idnode]=0
			while (j<5 and (int(self.neig[idnode,j]) != -1)):
				G=1./self.R[idnode,j]
				ng=int(self.neig[idnode,j])
				deltaT=self.temp[ng] - self.temp[idnode]
				dT_dt[idnode]+= G*deltaT
				j+=1
		print(time)
		print(dT_dt)


		#~ INITIALISATION 

Reservoir1=Reservoir()
#~ Reservoir1.systeme_init_cart()


		#~ INTEGRATION TEMPORELLE ET ECRITURE
		
listVar=[i for i in range(Reservoir1.Nptsx*Reservoir1.Nptsy)]
		
ProblemSize=Reservoir1.Nptsx*Reservoir1.Nptsy
AdapTimeStepBool=False
Duration=6.
MAXNTIMESTEP=101 
TIMESTEP=0.0002
METHODE='RK451'
ListOfVariablesToSave=listVar
SavedIteration=10
Error=1e-8

Problem=PbDef.NumericalProblem(Reservoir1.systeme_init_cart,Reservoir1.systeme_cond,METHODE,\
MAXNTIMESTEP,ProblemSize,TimeStep=TIMESTEP,\
Duration=Duration,AdaptativeTimeStep=AdapTimeStepBool,\
NbIterationSaved=SavedIteration,\
ListOfIdVar=ListOfVariablesToSave,AdatativeTimeStep_Error=Error)

Problem.SolveDifferentialProblem()

