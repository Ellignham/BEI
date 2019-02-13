
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

exec(open("./visu.py").read())
if os.path.isfile('./ResultArray.dat') : 
	temps, temp = lecture_champs('ResultArray.dat')
	Reservoir1.temp2d = temp[-1].reshape((50, 5))
	plot_champs(Reservoir1, Reservoir1.temp2d, temps[-1])

else :
	print('Run python pops.py before visualisation')
