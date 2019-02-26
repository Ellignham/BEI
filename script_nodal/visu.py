#!/usr/bin/python2.7

#Imports
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math as math
from matplotlib.mlab import griddata 
import matplotlib.colors as mplc

import csv
import os



#Class
# ~ from input import Input
# ~ from init import Init
# ~ from pops import Reservoir
# ~ from debug import Debug


# ~ def __init__(self):
	# ~ """
	# ~ Class used to plot and write to files the different arrays
	# ~ """
	
	

def lecture_champs(fname,TIMESTEP,SavedIteration,Duration,TimeFrequency,IterationFrequency,TypeFrequency):
	"""
	Reads the temporal evolution of the field contained in ResultArray.dat
	return temps, and temperature table
	"""
	data = np.loadtxt(fname, dtype='float', comments='#', delimiter=' ')
	if(TypeFrequency==1):
		temps=np.zeros((int(Duration/TimeFrequency)+1))
		temperature=np.zeros((int(Duration/TimeFrequency)+1))
		for i in range(0,int(Duration/TimeFrequency)):
			temps[i]=data[i*int(TimeFrequency/(TIMESTEP*SavedIteration)),0]
			temperature[i]=data[i*int(TimeFrequency/(TIMESTEP*SavedIteration)),1]
	else :
		temps=np.zeros((int(Duration/(TIMESTEP*IterationFrequency))+1))
		temperature=np.zeros((int(Duration/(TIMESTEP*IterationFrequency))+1))
		for i in range(0,int(Duration/(TIMESTEP*IterationFrequency))):
			temps[i]=data[i*IterationFrequency/SavedIteration,0]
			temperature[i]=data[i*IterationFrequency/SavedIteration,1]
	#temps=data[:,0]
	#temperature=data[:,1:]
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
	Plots a surface view of a cartesian field. The field must be written in the non-structured form defined in Reservoir
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
    Plots a surface view of a reservoir-like geometry field. The field must be written in the non-structured form defined in Reservoir
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

def save_champs_res(x, y, champs, time, cmax, cmin, j):
    """ 
    Plots a surface view of a reservoir-like geometry field. The field must be written in the non-structured form defined in Reservoir
    """
    fig1=plt.figure()
    # define grid.
    xi = np.linspace(-0.1, 0.6, 800)
    yi = np.linspace(-.1, 1.1, 1000)
    norm = mplc.Normalize(cmin, cmax)
    #~ v = np.linspace(cmin, cmax, 100, endpoint=True)
    bounds=np.linspace(cmin,cmax,100)
    # grid the data.
    zi = griddata(x,y, champs, xi,yi, interp='linear')
    # contour the gridded data, plotting dots at the nonuniform data points.
    CS = plt.contourf(xi, yi, zi,100, vmax=cmax, vmin=cmin, norm=norm, levels=bounds, extend='max')
    #~ plt.autumn()
    cbar = fig1.colorbar(CS)
    plt.clim(vmin=cmin, vmax=cmax)
    cbar.ax.set_ylabel('Temperature')
    # plot data points.
    #~ plt.scatter(x, y, marker='o', s=5, zorder=10)
    #~ plt.xlim(-0.1, 0.6)
    #~ plt.ylim(-0.1, 1.1)

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Temperature au bout de %(g)s secondes'%{'g' : time})
    plt.savefig('res_%(g)s.png'%{'g' : '%06d' %  j})
    plt.close()

def save_png(x,y,champs_tab,time_tab):
	nframes=30
	mult=len(time_tab) // nframes
	j=0
	if os.path.isdir('data'):
		print('Writing figures in data directory')
	else :
		print('Creating data directory')
		os.mkdir('data')
	os.chdir('data')
	while j <  len(time_tab) :
		print('Saving figure %(g)s out of %(h)s in data directory'%{'g' : '%08d' % j, 'h' :  str(len(time_tab))})
		save_champs_res(x, y, champs_tab[j], time_tab[j], champs_tab.max(), champs_tab.min(), j)
		j+=mult

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

def ecriture_csv(ProblemSize,temps,Reservoir):
	"""
	Reads and write the array in csv format for Paraview
	"""



	ArrayTemp=np.zeros((ProblemSize,4),dtype=float)
	liste=["x","y","z","Temperature"]
	
	if(Reservoir.mesh_type=='cart'):
		for j in range(0,len(temps)):
			f=open("ArrayTemp_cart/ArrayTemp_cart.csv.{}".format(j),"wb")
			ArrayTemp[:,0]=Reservoir.nodes[:,2]
			ArrayTemp[:,1]=Reservoir.nodes[:,1]
			ArrayTemp[:,2]=0
			ArrayTemp[:,3]=temp[j]
			writer=csv.writer(f,delimiter=',')
			writer.writerow(liste)
			writer.writerows(ArrayTemp)
	else:
		for j in range(0,len(temps)):
			f=open("ArrayTemp_tank/ArrayTemp_tank.csv.{}".format(j),"wb")
			ArrayTemp[:,0]=Reservoir.nodes[:,2]
			ArrayTemp[:,1]=Reservoir.nodes[:,1]
			ArrayTemp[:,2]=0
			ArrayTemp[:,3]=temp[j]
			writer=csv.writer(f,delimiter=',')
			writer.writerow(liste)
			writer.writerows(ArrayTemp)


#~ vmax=abs(zi).max(), vmin=-abs(zi).max(),


