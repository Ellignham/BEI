#!/usr/bin/python2.7
import os
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

Reservoir1=Reservoir()

if (Reservoir1.mesh_type=='cart'):
    Reservoir1.domain_cart()
    Reservoir1.init_domain_cart()
    Reservoir1.initemp_cart_x()
    Reservoir1.resistance_tank()
elif (Reservoir1.mesh_type=='tank'):
    Reservoir1.domain_tank()
    Reservoir1.init_domain_tank()
    Reservoir1.initemp_tank_x()
    Reservoir1.resistance_tank()



exec(open("./visu.py").read())
if os.path.isfile('./ResultArray.dat') : 
	if Reservoir1.mesh_type=='cart' :
		#~ Version cartesienne
		temps, temp = lecture_champs('ResultArray.dat')
		Reservoir1.temp2d = temp[-1].reshape((Reservoir1.Nptsy, Reservoir1.Nptsx))
		plot_champs_cart(Reservoir1, Reservoir1.temp2d, temps[-1])
	else :
		#~Version totale 
		temps, temperature, x, y = reconstruct_champs(Reservoir1, 'ResultArray.dat')
		#~ Reservoir1.temp2d=Reservoir1.temp
		plot_champs_res(x, y, Reservoir1.temp, temps[-1])
	

else :
	print('Run python pops.py before visualisation')
