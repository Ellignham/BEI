#!/usr/bin/python
"""
TO DO
-changer ordre boucles dans le remplissage du maillage, raffiner les extremites (avec une fonction ?) 
- essayer les coord polaires pour les extremites du maillage
- OXYGENE LIQUIDE
-ajouter fonction de phase
- gaz a psat
- on neglige la convection en 1e approximation

-Implementer les CL sur la resistance thermique + ajouter la resistance par convection
-Creer le bon maillage + adapter x,y,mesh
-Calculer h coeff de convection avec nusselt (ou Gr Pr cf. handbook)
-Calculer nusselt avec corelation adaptee au probleme (Gnielinski ? Dittus-Boelter ?)
-Lire le handbook

-Implementer des noeuds mobiles pour le suivi d'interface
-Implementer la biblioteque d integration (pas encore fournie)
-Implementer le calcul des param thermo en fonction des conditions P,T (params fixes pour le moment)
-Fichier lisibles par paraview

-sortie environ 10cm
"""

#Imports
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math as math


class Input:
    def __init__(self):
        """
        mesh : cell centered for temperature, pressure    

        Variables :

        debug       : if true the program will run on debug mode

        fluid       : type of fluid in the tank (used to compute rho, Cp, ...)
        tfluid_init : initial temperature of the fluid [K]
        pfluid_init : initial pressure of the fluid [Pa]
        ufluid_init : initial x velocity of the fluid [m/s]
        vfluid_init : initial y velocity of the fluid [m/s]

        cond        : If true thermal conduction will be computed
        conv        : If true thermal convection will be computed
    
        rho         : density [kg/m^3] (ignored if computed_parameters=True )
        k           : thermal condutivity [W/K/m] (ignored if computed_parameters=True )

        Lx          : size of the domain in the x direction [m]
        Ly          : size of the domain in the y direction [m]
        
        Nptsx       : number of nodes in the x direction
        Nptsy       : number of nodes in the y direction

        T1          : temperature of the walls [K]
        """
        #Debug Mode
        self.debug=True

        #Parameters of the fluid
        self.computed_parameters=False
        self.fluid="liquid_water"
        self.tfluid_init=293
        self.pfluid_init=10**5
        self.ufluid_init=1
        self.vfluid_init=1
        
        #Parameters of the heat transfer
        self.cond=True
        self.conv=False

        #Parameters of the fluid (with thermodynamic parameters locked)
        if not self.computed_parameters :
            self.rho=10**3
            self.k=0.5918

        #Parameters of the mesh
        self.Lx=1.
        self.Ly=10.
        self.Nptsx=10
        self.Nptsy=100

        #Parameters of the problem
        self.T1=313
