#!/usr/bin/python3

#Imports
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math as math

#Class
from input import Input
from init import Init
from pops import Reservoir
from debug import Debug
from visu import Visu

simu=Visu()
simu.domain_tank()
#simu.resistance()
#simu.write_temp()
#simu.write_pres()
#simu.plot_temp()
#simu.plot_pres()

#simu.domain_cart()
#simu.detect_neig()

if simu.debug :
    test=0
#    simu.write_debug()
plt.show()
