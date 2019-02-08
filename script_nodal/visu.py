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

class Visu(Debug):
    def __init__(self):
        """
        Class used to plot and write to files the different arrays
        """
        
        Debug.__init__(self)

    def plot_temp(self):
        """ 
        Plots a surface view of the temperature
        """
        plt.figure() 
        plt.imshow(self.temp.reshape((self.Nptsy,self.Nptsx)),extent=(self.x.min(), self.x.max(), self.y.max(), self.y.min()), interpolation='nearest',cmap=cm.gist_rainbow)
        plt.colorbar()
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Temperature')

    def write_temp(self):
        """
        Writes the temperature array to a file
        """
        file=open("temp.txt",'w+') 

        for j in range(0,self.Nptsy):
            for i in range(0,self.Nptsx):
                file.write(str(self.temp[j,i]) + ' ')
            file.write('\n')
        file.close()

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

    def write_pres(self):
        """
        Writes the pressure array to a file
        """
        file=open("press.txt",'w+') 

        for j in range(0,self.Nptsy):
            for i in range(0,self.Nptsx):
                file.write(str(self.pres[j,i]) + ' ')
            file.write('\n')
        file.close()







