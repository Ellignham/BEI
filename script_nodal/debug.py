#!/usr/bin/python3

#Imports
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math as math  

#Class
from input import Input
from init import Init


class Debug(Init):
    def __init__(self):
        """
        Class used to write to files the different arrays to hel debug the code
        """
        
        Init.__init__(self)
   
    def write_debug(self):
        """ 
        Write a file containing most varibles and constants to help debuging
        """

        file=open("debug.txt",'w+')
        file.write("---------------------------------------------------------" + "\n")
        file.write("| File containing most variables to help debuging" + "\n")
        file.write("---------------------------------------------------------" + "\n")
        file.write("\n")
        file.write("||  USER SET VARIBLES ||" + "\n")
        file.write("computed_parameters : " + str(self.computed_parameters) + "\n")
        file.write("fluid : " + str(self.fluid) + "\n")
        file.write("tfluid : " + str(self.tfluid_init) + "\n")
        file.write("pfluid : " + str(self.pfluid_init) + "\n")
        if not self.computed_parameters :
            file.write("rho : " + str(self.rho) + "\n")
        else :
            file.write("rho : IGNORED" + "\n")
        file.write("Lx : " + str(self.Lx) + "\n")
        file.write("Ly : " + str(self.Ly) + "\n")
        file.write("Nptsx : " + str(self.Nptsx) + "\n")
        file.write("Nptsy : " + str(self.Nptsy) + "\n")
        file.write("T1 : " + str(self.T1) + "\n")
        file.write("\n")
        
        file.write("||  COMPUTED VARIBLES ||" + "\n")
        file.write("\n")

        file.write("||  COMPUTED ARRAYS  ||" + "\n")
        
        file.write("x array" + "\n")
        for i in range(0,self.Nptsx):
            file.write(str(self.x[i]) + " ")
        file.write("\n")
        file.write("\n")
        
        file.write("y array" + "\n")
        for j in range(0,self.Nptsy):
            file.write(str(self.y[j]) + " ")
        file.write("\n")
        file.write("\n") 

        file.write("mesh array" + "\n")
        for j in range(0,self.Nptsy):
            file.write(str(self.mesh[0][j]) + " " )
        file.write("\n")
        for i in range(0,self.Nptsx):
            file.write(str(self.mesh[1][i]) + " " )
        file.write("\n")
        file.write("\n")
 
        file.write("temperature array" + "\n")
        for j in range(0,self.Nptsy):
            for i in range(0,self.Nptsx):
                file.write(str(self.temp[j,i]) + ' ')
            file.write('\n')
        file.write("\n")
 
        file.write("Pressure array" + "\n")
        for j in range(0,self.Nptsy):
            for i in range(0,self.Nptsx):
                file.write(str(self.pres[j,i]) + ' ')
            file.write('\n')
        file.write("\n")
        
        file.write("Velocity x array" + "\n")
        for j in range(0,self.Nptsy):
            for i in range(0,self.Nptsx):
                file.write(str(self.U[j,i]) + ' ')
            file.write('\n')
        file.write("\n")
        
        file.write("Velocity y array" + "\n")
        for j in range(0,self.Nptsy):
            for i in range(0,self.Nptsx):
                file.write(str(self.V[j,i]) + ' ')
            file.write('\n')
        file.write("\n")

        file.write("Thermal resistance x" + "\n")
        for j in range(0,self.Nptsy):
            for i in range(0,self.Nptsx+1):
                file.write(str(self.Rx[j,i]) + ' ')
            file.write('\n')
        file.write("\n")
        
        file.write("Thermal resistance y" + "\n")
        for j in range(0,self.Nptsy+1):
            for i in range(0,self.Nptsx):
                file.write(str(self.Ry[j,i]) + ' ')
            file.write('\n')
        file.write("\n")
        file.write("\n")

        file.close()

