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

plt.close('all')

Reservoir1=Reservoir()
#~ 
#~ if (Reservoir1.mesh_type=='cart'):
    #~ Reservoir1.domain_cart()
    #~ Reservoir1.init_domain_cart()
    #~ Reservoir1.initemp_cart_x()
    #~ Reservoir1.resistance_cart()
#~ elif (Reservoir1.mesh_type=='tank'):
    #~ Reservoir1.domain_tank()
    #~ Reservoir1.init_domain_tank()
    #~ Reservoir1.initemp_tank_x()
    #~ Reservoir1.resistance_tank()
	#~ 


exec(open("./visu.py").read())
if os.path.isfile('./ResultArray.dat') : 
	if Reservoir1.mesh_type=='cart' :
		#~ Version cartesienne
		temps, temp = lecture_champs('ResultArray.dat')
		Reservoir1.temp2d = temp[-1].reshape((Reservoir1.Nptsy, Reservoir1.Nptsx))
		plot_champs_cart(Reservoir1, Reservoir1.temp2d, temps[-1])
	else :
		
		#~ Reservoir1=Reservoir()
		Reservoir1.domain_tank()
		Reservoir1.init_phase()
		Reservoir1.init_domain()
		
		#~Version totale 
		temps, temperature, x, y = reconstruct_champs(Reservoir1, 'ResultArray.dat')
		print(temps[-1])
		Reservoir1.hauteur_interface(temps[-1])
		save_png(Reservoir1, x,y,temperature,temps)
		# ~ plot_temp_int(Reservoir1, x, y, temperature[-1], temps[-1]+Reservoir1.time_init)

else :
	print('Run python script.py before visualisation')
