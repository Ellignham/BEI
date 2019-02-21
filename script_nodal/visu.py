#!/usr/bin/python2.7

#Imports
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math as math
from matplotlib.mlab import griddata
import csv

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
	return temps, and temperature table
	"""
	data = np.loadtxt(fname, dtype='float', comments='#', delimiter=' ')
	temps=data[:,0]
	temperature=data[:,1:]
	# ~ print temperature
	return temps, temperature
	
	
	

def reconstruct_champs(dom,fname):
	"""
	Reads the temporal evolution of the field contained in ResultArray.dat
	return temps, and temperature table
	"""
	data = np.loadtxt(fname, dtype='float', comments='#', delimiter=' ')
	temps=data[:,0]
	temperature=data[:,1:]
	for idnode in range(len(temperature[0])):
		x=np.copy(dom.nodes[:,2])
		y=np.copy(dom.nodes[:,1])
	return temps, temperature, x, y
	


def plot_champs_cart(dom, champs, time):
	""" 
	Plots a surface view of a cartesian field. The field must be of the form defined in Reservoir
	"""
	plt.figure() 
	plt.contourf(dom.meshx, dom.meshy, champs)
	plt.colorbar()
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Temperature au bout de %(g) secondes'%{'g' : time})
	plt.show()
	#~ 
	#~ 
#~ def plot_coupe_lx2(x, y, champs, dom, time) :
	#~ plt.figure('Profil de temperature en coupe')
	#~ milieu = (


def plot_champs_res(x, y, champs, time):
	""" 
	Plots a surface view of a reservoir-like geometry field. The field must be of the form defined in Reservoir
	"""

	
	fig1=plt.figure()	

	# define grid.
	xi = np.linspace(-2.5, 2.5, 100)
	yi = np.linspace(-.1, 10.1, 200)
	# grid the data.
	zi = griddata(x,y, champs, xi,yi, interp='linear')
	# contour the gridded data, plotting dots at the nonuniform data points.
	#~ CS = plt.contour(xi, yi, zi, 15, linewidths=0.5, colors='k')
	#~ CS = plt.contourf(xi, yi, zi, 100,
					  #~ vmax=abs(zi).max(), vmin=-abs(zi).max(), cmap=plt.cm.bone, origin='upper')
	CS = plt.contourf(xi, yi, zi, 100,
					   cmap=plt.cm.rainbow)
	cbar = fig1.colorbar(CS)
	cbar.ax.set_ylabel('champs')
	#~ plt.colorbar()  # draw colorbar
	# plot data points.
	plt.scatter(x, y, marker='o', s=5, zorder=10)
	#~ plt.xlim(-2, 2)
	#~ plt.ylim(-2, 2)

	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Temperature au bout de %(g) secondes'%{'g' : time})
	plt.show()





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

def ecriture_csv(ProblemSize,temps,Reservoire):
	"""
	Reads the array for each time
	return an array of the field, and temperature for each coordinate
	"""
	ArrayTemp=np.zeros((ProblemSize,3),dtype=float)
	liste=["x","y","Temperature"]
	for j in range(0,len(temps)):
		f=open("ArrayTemp_{}.csv".format(j),"wb")
		ArrayTemp[:,0]=Reservoire.domain[:,0]
		ArrayTemp[:,1]=Reservoire.domain[:,1]
		ArrayTemp[:,2]=temp[j]
		writer=csv.writer(f,delimiter=',')
		writer.writerow(liste)
		writer.writerows(ArrayTemp)



#~ vmax=abs(zi).max(), vmin=-abs(zi).max(),


