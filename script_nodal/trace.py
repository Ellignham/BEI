
import os
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math as math
import sys
LibPath='../../Sujet_29/'
sys.path.append(LibPath)
import NumericalProblemClass as PbDef


print 'toto'
#Class
from input import Input
from init import Init
from pops import Reservoir


Reservoir1=Reservoir()
typ='cart'

exec(open("./visu.py").read())
if os.path.isfile('./ResultArray.dat') : 
	if typ=='cart' :
		#~ Version cartesienne
		temps, temp = lecture_champs('ResultArray.dat')
		Reservoir1.temp2d = temp[-1].reshape((Reservoir1.Nptsy, Reservoir1.Nptsx))
		plot_champs_cart(Reservoir1, Reservoir1.temp2d, temps[-1])
	else :
		#~Version totale 
		Reservoir1.temp2d = temp[-1]
		plot_champs_res(Reservoir1, Reservoir1.temp2d, temps[-1])
	

else :
	print('Run python pops.py before visualisation')
