#!/usr/bin/python3

#Imports
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math as math

#Class
# ~ from input import Input
# ~ from init import Init
# ~ from pops import Reservoir
# ~ from debug import Debug


# ~ def __init__(self):
	# ~ """
	# ~ Class used to plot and write to files the different arrays
	# ~ """
	
	

def lecture_champs(fname):
	"""
	Reads the temporal evolution of the field contained in ResultArray.dat
	"""
	data = np.loadtxt(fname, dtype='float', comments='#', delimiter=' ')
	print data
	# ~ return data
	


def plot_champs(champs):
	""" 
	Plots a surface view of a field. The field must be of the form defines in Reservoir
	"""
	plt.figure() 
	plt.imshow(self.temp.reshape((self.Nptsy,self.Nptsx)),extent=(self.x.min(), self.x.max(), self.y.max(), self.y.min()), interpolation='nearest',cmap=cm.gist_rainbow)
	plt.colorbar()
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Temperature')



def plot_pres(self):
	""" 
	Plots a surface view of the pressure
	"""
	plt.figure()
	plt.imshow(self.pres.reshape((self.Nptsy,self.Nptsx)),extent=(self.x.min(), self.x.max(), self.y.max(), self.y.min()), interpolation='nearest',cmap=cm.gist_rainbow)
	plt.colorbar()
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Pressure')








