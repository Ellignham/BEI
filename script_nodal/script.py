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
from pops import Reservoir


		#~ INITIALISATION 

Reservoir1=Reservoir()

		#~ TIME INTEGRATION AND WRITING
if (Reservoir1.mesh_type=='tank'):		
    ProblemSize=Reservoir1.Nptsx*Reservoir1.Nptsy+2*(Reservoir1.Nptsx-1)*Reservoir1.ntheta
elif (Reservoir1.mesh_type=='cart'):
    ProblemSize=Reservoir1.Nptsx*Reservoir1.Nptsy


listVar=[i for i in range(ProblemSize)]
AdapTimeStepBool=False
Duration=50.
MAXNTIMESTEP=1000000
TIMESTEP=0.002
METHODE='Euler'
ListOfVariablesToSave=listVar
SavedIteration=1000
Error=1e-8

Reservoir1.dt = TIMESTEP

Problem=PbDef.NumericalProblem(Reservoir1.systeme_init,Reservoir1.systeme_cond,METHODE,\
MAXNTIMESTEP,ProblemSize,TimeStep=TIMESTEP,\
Duration=Duration,AdaptativeTimeStep=AdapTimeStepBool,\
NbIterationSaved=SavedIteration,\
ListOfIdVar=ListOfVariablesToSave,\
AdatativeTimeStep_Error=Error)


Problem.SolveDifferentialProblem()


		# ~ TRACE DES DONNEES
				
exec(open("./visu.py").read())
TypeFrequency=1 # =0 : frequence d'iteration , =1 : frequence de temps
TimeFrequency=1.0
IterationFrequency=200
temps, temp = lecture_champs('ResultArray.dat',TIMESTEP,SavedIteration,Duration,TimeFrequency,IterationFrequency,TypeFrequency)


temps, temperature, x, y = reconstruct_champs(Reservoir1, 'ResultArray.dat')
#~ Reservoir1.temp2d = temp[i].reshape((50, 5))
plot_temp_int(Reservoir1, x, y, Reservoir1.phi, temps[-1]+Reservoir1.time_init)




#~ ecriture_csv(ProblemSize,temps,Reservoir1)	


#for i in range(len(temps)):


#for i in range(len(temps)) :
#	
#	Reservoir1.temp2d = temp[i].reshape((50, 5))
#	plot_champs(Reservoir1, Reservoir1.temp2d, temps[i])
