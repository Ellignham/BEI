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
from pops import Reservoir


		#~ INITIALISATION 

Reservoir1=Reservoir()

		#~ INTEGRATION TEMPORELLE ET ECRITURE
		
ProblemSize=Reservoir1.Nptsx*Reservoir1.Nptsy
listVar=[i for i in range(ProblemSize)]
# ~ listVar=[0,1]
AdapTimeStepBool=False
Duration=60.
MAXNTIMESTEP=100001 
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
