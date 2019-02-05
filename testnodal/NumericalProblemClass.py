import os, sys
import numpy as np
import copy
import time

class NumericalProblem:
    """
    Class NumericalProblem
    Built to allow integrating differential equations written in an 
    input file by users    
    To display the complete description of the class     
    use method DisplayHelp()    
    To use it, some requirements are mandatory
    1 - Define two functions to specify initial conditions and 
    differential systems 
    to put it as arguments into __init__ function 
    2 - Define every variable required to define time integration 
    method (scheme, timestep, duration, variables to save,...)
    3 - Solve the problem by calling SolveDifferentialProblem function
    
    Major functions :
    - __init__ : allow to define every variables and arrays to set the problem
    - SolveDifferentialProblem : Main function of the class to build 
    time integration
    
    Functions related to Data storage :
    - RecordData
    - PickData
    - FormateData
    
    Functions related to time integration step by step
    - NextStepComputation
    - UpdateTimeIntegration
    
    Functions for every explicit scheme 
    - ExplicitEulerStep
    - ExplicitRK2
    - ExplicitRK4
    - ExplicitRK41
    - ExplicitRK45
    - ExplicitRK451
    
    Function to display help :
    - DisplayHelp   
    
    """
    
    def __init__(self,InitFunc,DiffFunc,SchemeName,NbIteration,NbVariables,\
                 TimeStep=1e-3,\
                 Duration=0.0,AdaptativeTimeStep=False,NbIterationSaved=1000,\
                 ListOfIdVar=[],AdatativeTimeStep_Error=1e-8):
        """     
        
        Function to initialize the numerical problem by specifying :
        Inputs       : InitFunc ; function to define initial conditions written 
                       by users into InputFile.py
        Inputs       : DiffFunc ; function to define differential 
        equations written by users into InputFile.py
        Inputs       : SchemeName ; STR Name of the scheme into 
        explicit and implicit schemes available 
        Inputs       : NbIteration ; INT Max iteration number
        Inputs (OPT) : TimeStep ; FLT Timestep (1e-3 by default)
        Inputs (OPT) : NbVariables ; INT Size of the problem
        Inputs (OPT) : Duration ; FLT Duration of the simulation 
        Inputs (OPT) : AdaptativeTimeStep ; BOOL Boolean to make active adapatative time step
        Inputs (OPT) : NbIterationSaved : INT Step of iteration to save data
        Inputs (OPT) : ListOfIdVar ; List of Id of varables to save
        
        Outputs : nothing () 
        
        $$$$$$$$$$$$$$$$$$$$$
        """    
        
        self.NbIterationSavedMax=NbIteration/NbIterationSaved+1
    
        self.AlphatimeStep=1.0
        self.AlphaRK2=0.5

        self.TimeStep=TimeStep

        self.NbIterationMax=100000


        self.CurrentIteration=0
        self.Clock=0.0
        self.NbOfSavedVar=0

        self.zero=1e-30
        
        self.AdaptativeError=AdatativeTimeStep_Error
        self.JacobianError=1e-15
        self.MaxIte=50
        self.hf=0.5**20
        self.ha=0.5**40
        self.LastCallRecord=0
        self.LastCallRecordResult=0
      
        self.InitializationFunction=InitFunc
        self.DifferentialFunction  =DiffFunc
        
        if (SchemeName[:2]=='RK' or SchemeName=='Euler'):
            self.ExplicitBool=True
            self.ImplicitBool=False
        else:
            self.ExplicitBool=False
            self.ImplicitBool=True            
        
        self.SchemeName=SchemeName
        self.NbIteration=NbIteration
        self.Duration=Duration

        self.NbIterationSaved=NbIterationSaved
        self.AdaptativeTimeStep=AdaptativeTimeStep
        self.NbVariables=NbVariables
        self.dY=np.zeros(self.NbVariables,dtype='float64',order='F')
        


        
        
        if (self.SchemeName=='ImplicitRK3'):
            self.K1=np.zeros(self.NbVariables,dtype='float64',order='F')
            self.K2=np.zeros(self.NbVariables,dtype='float64',order='F')
            self.K1[:]=1.0
            self.K2[:]=1.0
            
            self.K_1Flux=np.zeros(self.NbVariables,dtype='float64',order='F')
            self.K_2Flux=np.zeros(self.NbVariables,dtype='float64',order='F')
            self.K_Flux=np.zeros(2*self.NbVariables,dtype='float64',order='F')
                        
        self.dy1=np.zeros(self.NbVariables,dtype='float64',order='F')
        self.dy2=np.zeros(self.NbVariables,dtype='float64',order='F')
        self.dy3=np.zeros(self.NbVariables,dtype='float64',order='F')
        self.dy4=np.zeros(self.NbVariables,dtype='float64',order='F')
        self.dy5=np.zeros(self.NbVariables,dtype='float64',order='F')
        self.dy6=np.zeros(self.NbVariables,dtype='float64',order='F')
        self.dy5=np.zeros(self.NbVariables,dtype='float64',order='F')
        self.dy_tmp_45=np.zeros(self.NbVariables,dtype='float64',order='F')        
        self.dy_tmp_451=np.zeros(self.NbVariables,dtype='float64',order='F')     
        
        self.Fn=np.zeros(self.NbVariables,dtype='float64',order='F')
        self.Fnp1=np.zeros(self.NbVariables,dtype='float64',order='F')
        self.Fnm1=np.zeros(self.NbVariables,dtype='float64',order='F')
        self.Fnm2=np.zeros(self.NbVariables,dtype='float64',order='F')

        self.Flux_m1=np.zeros(self.NbVariables,dtype='float64',order='F')
        self.Flux_p1=np.zeros(self.NbVariables,dtype='float64',order='F')
        self.Flux_TMP=np.zeros(self.NbVariables,dtype='float64',order='F')
        
        self.dYtmp=np.zeros(self.NbVariables,dtype='float64',order='F')   
        self.Ytmp=np.zeros(self.NbVariables,dtype='float64',order='F')        
        self.Ynm1=np.zeros(self.NbVariables,dtype='float64',order='F')
        self.Ynm2=np.zeros(self.NbVariables,dtype='float64',order='F')
        self.Yn=np.zeros(self.NbVariables,dtype='float64',order='F')

        self.dy_TestOrder=np.zeros(self.NbVariables,dtype='float64',order='F')

        self.NbIterationSaved=NbIterationSaved      
        
        self.Ynp1=np.zeros(self.NbVariables)
        if (self.ImplicitBool):
            self.JacobianMatrix=np.zeros((self.NbVariables,self.NbVariables),\
                                         dtype='float64',order='F')
                
        self.NbOfSavedVar=len(ListOfIdVar)
        self.ListOfIdVar=ListOfIdVar
        self.TimeVector=np.zeros(self.NbIterationSavedMax,dtype='float64')
        self.Residuals=np.zeros(self.NbIterationSavedMax,dtype='float64')
        
        self.ResultArray=np.zeros((self.NbIterationSavedMax,\
                                   self.NbOfSavedVar+1),dtype='float64')
        self.HastTable=np.zeros((len(ListOfIdVar),2)).astype(int)
        i=0
        for id in ListOfIdVar:
            self.HastTable[i,1]=id
            self.HastTable[i,0]=i
            i=i+1
            
        self.DefineDict() 
        self.DefineSchemeParams()

    def DefineDict(self):
        """
        
        Function to define dictionnaries used to select integration 
        schemes according to user inputs
        
        $$$$$$$$$$$$$$$$$$$$$
        """
        self.SchemeFunctionsDic={}
        
        self.SchemeFunctionsDic['Euler']=self.ExplicitEulerStep
        self.SchemeFunctionsDic['RK2']=self.ExplicitRK2
        self.SchemeFunctionsDic['RK3']=self.ExplicitRK3
        self.SchemeFunctionsDic['RK4']=self.ExplicitRK4
        self.SchemeFunctionsDic['RK41']=self.ExplicitRK41
        self.SchemeFunctionsDic['RK45']=self.ExplicitRK45
        self.SchemeFunctionsDic['RK451']=self.ExplicitRK451
        
        self.SchemeFunctionsDic['ImplicitEuler']=self.ImplicitIntegration
        self.SchemeFunctionsDic['ImplicitRK2']=self.ImplicitIntegration
        self.SchemeFunctionsDic['ImplicitAdam3o']=self.ImplicitIntegration
        self.SchemeFunctionsDic['ImplicitAdam4o']=self.ImplicitIntegration
        self.SchemeFunctionsDic['ImplicitBDF2']=self.ImplicitIntegration
        self.SchemeFunctionsDic['ImplicitBDF3']=self.ImplicitIntegration
        self.FluxFunctionsForImplicit={}
        
        self.FluxFunctionsForImplicit['Euler']=None
        self.FluxFunctionsForImplicit['RK2']=None
        self.FluxFunctionsForImplicit['RK3']=None
        self.FluxFunctionsForImplicit['RK4']=None
        self.FluxFunctionsForImplicit['RK41']=None
        self.FluxFunctionsForImplicit['RK45']=None
        self.FluxFunctionsForImplicit['RK451']=None
        
        self.FluxFunctionsForImplicit['ImplicitEuler']=self.ImplicitEulerFlux
        self.FluxFunctionsForImplicit['ImplicitRK2']=self.ImplicitRK2Flux        
        self.FluxFunctionsForImplicit['ImplicitAdam3o']=self.ImplicitAdam3o   
        self.FluxFunctionsForImplicit['ImplicitAdam4o']=self.ImplicitAdam4o
        self.FluxFunctionsForImplicit['ImplicitBDF2']=self.ImplicitBDF2
        self.FluxFunctionsForImplicit['ImplicitBDF3']=self.ImplicitBDF3
        
        self.SchemeFunctionsDicForAdaptativeTimeStep={}
        
        self.SchemeFunctionsDicForAdaptativeTimeStep['Euler']=self.ExplicitRK2
        self.SchemeFunctionsDicForAdaptativeTimeStep['RK2']=self.ExplicitRK3
        self.SchemeFunctionsDicForAdaptativeTimeStep['RK3']=self.ExplicitRK45
        self.SchemeFunctionsDicForAdaptativeTimeStep['RK4']=self.ExplicitRK451
        self.SchemeFunctionsDicForAdaptativeTimeStep['RK41']=self.ExplicitRK451
        self.SchemeFunctionsDicForAdaptativeTimeStep['RK45']=self.ExplicitRK451

        self.OrderSupForAdaptativeTimeStep={}

        self.OrderSupForAdaptativeTimeStep['Euler']=1.0/2.0
        self.OrderSupForAdaptativeTimeStep['RK2']=1.0/3.0
        self.OrderSupForAdaptativeTimeStep['RK3']=1.0/4.0
        self.OrderSupForAdaptativeTimeStep['RK4']=1.0/5.0
        self.OrderSupForAdaptativeTimeStep['RK41']=1.0/5.0
        self.OrderSupForAdaptativeTimeStep['RK45']=1.0/5.0

        self.OrderInfForAdaptativeTimeStep={}

        self.OrderInfForAdaptativeTimeStep['Euler']=1.0
        self.OrderInfForAdaptativeTimeStep['RK2']=1.0/2.0
        self.OrderInfForAdaptativeTimeStep['RK3']=1.0/3.0
        self.OrderInfForAdaptativeTimeStep['RK4']=1.0/4.0
        self.OrderInfForAdaptativeTimeStep['RK41']=1.0/4.0
        self.OrderInfForAdaptativeTimeStep['RK45']=1.0/4.0
        
        if (self.AdaptativeTimeStep):
            self.AdaptativeOrderSup=\
            self.OrderSupForAdaptativeTimeStep[self.SchemeName]
            self.AdaptativeOrderInf=\
            self.OrderInfForAdaptativeTimeStep[self.SchemeName]
            self.ExplicitFunctionForAdaptTimeStep=\
            self.SchemeFunctionsDicForAdaptativeTimeStep[self.SchemeName]
        
        self.TimeIntegrationFunction=self.SchemeFunctionsDic[self.SchemeName]
        self.BuildFluxFunction=self.FluxFunctionsForImplicit[self.SchemeName]

    def DefineSchemeParams(self):
        """
        
        Function to define dictionnaries used to set scheme parameters 
        especially for divisions
        
        $$$$$$$$$$$$$$$$$$$$$
        """    
        # ParamScheme
        
        self.inv2=0.5
        self.inv4=0.25
        self.inv6=1.0/6.0
        self.inv3=1.0/3.0
        self.inv8=1.0/8.0
        self.inv13=1.0/13.0
        self.inv32=1.0/32.0
        self.inv2197=1.0/2197.0
        self.inv216=1.0/216.0
        self.inv513=1.0/513.0
        self.inv4104=1.0/4104.0
        self.inv27=1.0/27.0
        self.inv2565=1.0/2565.0
        self.inv40=1.0/40.0
        
        self.inv11=1.0/11.0
        
        self.inv135=1.0/135.0
        self.inv12825=1.0/12825.0
        self.inv56430=1.0/56430.0
        self.inv50=1.0/50.0
        self.inv55=1.0/55.0
        self.inv12=1.0/12.0
        self.inv24=1.0/24.0

    def RecordData(self):
        """
        
        Function to record data selected by the user every 
        iteration (also selected by the user)
        
        Store data into self.ResultArray (first column : time)
        
        Inputs : none
        Outputs : none
        
        $$$$$$$$$$$$$$$$$$$$$
        """
        if (self.LastCallRecord<self.NbIterationSavedMax-1):
            self.TimeVector[self.LastCallRecord]=self.Clock      
            self.Residuals[self.LastCallRecord]=np.sum(np.abs(self.dY))
            if (self.LastCallRecord%1==0):
                i=0
                self.ResultArray[self.LastCallRecordResult,0]=self.Clock
                for Id in self.ListOfIdVar:
                    self.ResultArray[self.LastCallRecordResult,i+1]=self.Yn[self.HastTable[i,1]]
                    i=i+1
                self.LastCallRecordResult=self.LastCallRecordResult+1    
            self.LastCallRecord=self.LastCallRecord+1
                 
    def FormateData(self):
        """
        
        Function to properly cut and close result arrays once integration performed.
        Write data into ResultArray.dat file and the residuals (sum of differential vector)
        into Residual.dat file.
        
        Store data into self.ResultArray (first column : time)
        
        Inputs : none
        Outputs : none
        
        $$$$$$$$$$$$$$$$$$$$$
        """
        Flag=1

        self.TimeVector=copy.deepcopy(self.TimeVector[:self.LastCallRecord+1-Flag])
        self.ResultArray=copy.deepcopy(self.ResultArray[:self.LastCallRecordResult+1-Flag,:])

        self.Residuals=copy.deepcopy(self.Residuals[:self.LastCallRecord+1-Flag])
        f=open('Residual.dat','w+')
        np.savetxt(f,np.vstack((self.TimeVector,self.Residuals)).T)
        f.close()
        f=open('ResultArray.dat','w+')
        np.savetxt(f,self.ResultArray)
        f.close()  
 
    def UpdateTimeIntegration(self,dt):
        """
        
        Function to properly cut and close result arrays once integration performed.
        Write data into ResultArray.dat file and the residuals (sum of differential vector)
        into Residual.dat file.
        
        Store data into self.ResultArray (first column : time)
        
        Inputs : none
        Outputs : none
        
        $$$$$$$$$$$$$$$$$$$$$
        """    
        self.TimeStep = dt	
        self.Ynm2[:]  = self.Ynm1[:]	
        self.Ynm1[:]  = self.Yn[:]	
        self.Yn[:]    = self.Ynp1[:]

        

        
        self.Fnm2[:]=self.Fnm1[:]
        
        self.Fnm1[:]=self.Fn[:]
        self.Fn[:]=self.Fnp1[:]
        
        self.Clock    = self.Clock + self.TimeStep
        self.CurrentIteration+=1	
 
        self.DifferentialFunction(self.Yn[:],self.dY[:],time=self.Clock)
   
    def DefineInitialCondition(self):
        """
        
        Function to initialize the solution vector with the function 
        specified by the user
        
        Inputs : none
        Outputs : none
        
        $$$$$$$$$$$$$$$$$$$$$
        """ 
        self.InitializationFunction(self.Yn)
                
    def NextStepComputation(self,y_vector,dy_vector,dh):
        """
        
        Function to compute a next vector according to a step, 
        the previous vector and the differential one
        
        y(n+1) = y(n) + h dy(n)
        
        Function called by explicit schemes only
        
        Inputs : y(:) array of FLT which is the vector at the current position
        Inputs : dy(:) array of FLT which is the differential vector 
        Inputs : h the integration step
        Outputs : y(:) + h dy(:) array of float which is the vector at the predicted position
        
        $$$$$$$$$$$$$$$$$$$$$
        """ 
        return y_vector[:]+dy_vector[:]*dh
        
    def ExplicitEulerStep(self):
        """
        
        Function to build explicit euler scheme
        1 step
        1th order

        Inputs : none
        Outputs : none
        
        $$$$$$$$$$$$$$$$$$$$$
        """
        self.DifferentialFunction(self.Yn,self.dy1,time=self.Clock)
        self.Ynp1=self.NextStepComputation(self.Yn,self.dy1,self.TimeStep)
    
    def ExplicitRK2(self):
        """
        
        Function to build explicit Runge Kutta 2 scheme
        2 step
        2th order

        Inputs : none
        Outputs : none
        
        $$$$$$$$$$$$$$$$$$$$$
        """    
        self.DifferentialFunction(self.Yn,self.dy1,time=self.Clock)
        self.dYtmp[:]=self.NextStepComputation(self.Yn,self.dy1[:]*self.inv2,self.TimeStep)
        self.DifferentialFunction(self.dYtmp,self.dy2,time=self.Clock+self.TimeStep*self.inv2)
        self.Ynp1[:]=self.NextStepComputation(self.Yn,self.dy2,self.TimeStep)
  
    def ExplicitRK3(self):
        """
        
        Function to build explicit Runge Kutta 3 scheme
        3 step
        3th order
        
        Inputs : none
        Outputs : none
        
        $$$$$$$$$$$$$$$$$$$$$
        """     
        self.DifferentialFunction(self.Yn,self.dy1,time=self.Clock)
        self.dYtmp[:]=self.dy1[:]
        self.Ytmp[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep*self.inv2)
        
        self.DifferentialFunction(self.Ytmp,self.dy2,time=self.Clock+self.TimeStep*self.inv2)
        self.dYtmp[:]=-self.dy1[:]+2.0*self.dy2[:]
        
        self.Ytmp[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep)
        self.DifferentialFunction(self.Ytmp,self.dy3,time=self.Clock+self.TimeStep)
        self.dYtmp[:]=(self.dy1[:]+4.0*self.dy2[:]+self.dy3)*self.inv6

         
        self.Ynp1[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep) 
  
    def ExplicitRK4(self):
        """
        
        Function to build explicit Runge Kutta 4 scheme
        4 step
        4th order

        Inputs : none
        Outputs : none
        
        $$$$$$$$$$$$$$$$$$$$$
        """      
        self.DifferentialFunction(self.Yn,self.dy1,time=self.Clock)
        self.dYtmp[:]=self.dy1[:]*self.inv2
        self.Ytmp[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep)
        self.DifferentialFunction(self.Ytmp,self.dy2,time=self.Clock+self.TimeStep*self.inv2)
        self.dYtmp[:]=self.dy2[:]*self.inv2
        self.Ytmp[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep)
        self.DifferentialFunction(self.Ytmp,self.dy3,time=self.Clock+self.TimeStep*self.inv2)
        self.Ytmp[:]=self.NextStepComputation(self.Yn,self.dy3,self.TimeStep)
        self.DifferentialFunction(self.Ytmp,self.dy4,time=self.Clock+self.TimeStep)
        self.dYtmp[:]=self.inv6*(self.dy1[:]+2.0*self.dy2[:]+2.0*self.dy3[:]+self.dy4[:])
        self.Ynp1[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep)	 
               
    def ExplicitRK41(self):
        """
        
        Function to build explicit Runge Kutta 4 scheme
        (formulation slightly different)
        4 step
        4th order

        Inputs : none
        Outputs : none
        
        $$$$$$$$$$$$$$$$$$$$$
        """      
        self.DifferentialFunction(self.Yn,self.dy1,time=self.Clock)
        self.dYtmp[:]=self.dy1[:]*self.inv3
        self.Ytmp[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep)
        
        self.DifferentialFunction(self.Ytmp,self.dy2,time=self.Clock+self.TimeStep*self.inv3)
        self.dYtmp[:]=-self.dy1[:]*self.inv3+self.dy2[:]
        
        self.Ytmp[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep)
        self.DifferentialFunction(self.Ytmp,self.dy3,time=self.Clock+self.TimeStep*2.0*self.inv3)
        self.dYtmp[:]=(self.dy1[:]-self.dy2[:]+self.dy3)
        self.Ytmp[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep)
        self.DifferentialFunction(self.Ytmp,self.dy4,time=self.Clock+self.TimeStep)
        self.dYtmp[:]=self.inv8*(self.dy1[:]+3.0*self.dy2[:]+3.0*self.dy3[:]+self.dy4[:])
         
        self.Ynp1[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep) 
             
    def ExplicitRK45(self):
        """
        
        Function to build explicit Runge Kutta 45 scheme
        5 step
        4-5th order

        Inputs : none
        Outputs : none
        
        $$$$$$$$$$$$$$$$$$$$$
        """  
            
        self.DifferentialFunction(self.Yn,self.dy1,time=self.Clock)
        self.dYtmp[:]=self.inv4*self.dy1[:]
        self.Ytmp[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep)
        
        self.DifferentialFunction(self.Ytmp,self.dy2,time=self.Clock+self.inv4*self.TimeStep)
        self.dYtmp[:]=3.0*self.inv32*self.dy1[:]+9.0*self.inv32*self.dy2[:]
        self.Ytmp[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep)
        
        self.DifferentialFunction(self.Ytmp,self.dy3,time=self.Clock+3.0*self.inv8*self.TimeStep)
        self.dYtmp[:]=1932.0*self.inv2197*self.dy1[:]-7200.0*self.inv2197*self.dy2[:]+7296.0*self.inv2197*self.dy3[:]
        self.Ytmp[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep)
        
        self.DifferentialFunction(self.Ytmp,self.dy4,time=self.Clock+12.0*self.inv13*self.TimeStep)
        self.dYtmp[:]=439.0*self.inv216*self.dy1[:]-8.0*self.dy2[:]+3680.0*self.inv513*self.dy3[:]-845.0*self.inv4104*self.dy4[:]
        self.Ytmp[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep)    
        
        self.DifferentialFunction(self.Ytmp,self.dy5,time=self.Clock+self.TimeStep)

        self.dYtmp[:]=-8.0*self.inv27*self.dy1[:]+2.0*self.dy2[:]-3544.0*self.inv2565*self.dy3[:]+1859*self.inv4104*self.dy4[:]-11.0*self.inv40*self.dy5[:]
        self.Ytmp[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep)    
        self.DifferentialFunction(self.Ytmp,self.dy6,time=self.Clock+self.TimeStep*self.inv2)
        
        self.dYtmp[:]=(25.0*self.inv216*self.dy1[:]+1408.0*self.inv2565*self.dy3[:]+2197.0*self.inv4104*self.dy4[:]-0.20*self.dy5[:])
        self.dy_tmp_45[:]=self.dYtmp[:]
        self.Ynp1[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep)
        
    def ExplicitRK451(self):
        """
        
        Function to build explicit Runge Kutta 451 scheme
        (this function calls RK 45 scheme)
        6 step
        5th order

        Inputs : none
        Outputs : none
        
        $$$$$$$$$$$$$$$$$$$$$
        """    
        self.ExplicitRK45()
        self.dYtmp[:]=(16.0*self.inv135*self.dy1[:]+6656.0*self.inv12825*self.dy3[:]+28561.0*self.inv56430*self.dy4[:]-9.0*self.inv50*self.dy5[:]+2.0*self.inv55*self.dy6[:])    
        self.dy_tmp_451[:]=self.dYtmp[:]
        self.Ynp1[:]=self.NextStepComputation(self.Yn,self.dYtmp,self.TimeStep)

    def AdaptTimeStep(self):
        """
        
        Function to compute the best time step if Adatapative time step is active
        based on two explicit schemes, the function computes the best time step to limit 
        the troncature error to an user's defined value
        
        The use of this option is costly because it duplicates the number of calls
        The efficiency depends on the scheme and the error (the error has to be decreased with the order
        of the scheme). 

        Inputs : none
        Outputs : none
        
        $$$$$$$$$$$$$$$$$$$$$
        """    
        
        estimator_LowOrder=self.dy_TestOrder[:]*self.TimeStep
        
        estimator_HighOrder=self.dYtmp[:]*self.TimeStep   
        
        AbsError=np.abs(estimator_HighOrder-estimator_LowOrder)
        MaxError=np.max(AbsError) 
        scale=1.0

        
        if (MaxError>self.zero):
            
            if(self.AdaptativeError>=MaxError):
                
                scale=((self.AdaptativeError/(MaxError))**(self.AdaptativeOrderSup)*(self.AlphatimeStep))
            else:    
                scale=((self.AdaptativeError/(MaxError))**(self.AdaptativeOrderInf)*(self.AlphatimeStep))    
        self.TimeStep=scale*self.TimeStep

    def ImplicitEulerFlux(self,Ynp1):
        """
        
        Function to build implicit Euler flux
        Write a flux vector to be equal to zero
        1th order

        Inputs : Ynp1[:] ; array of FLT : Variables vector at the computed position x+h
        Outputs : none
        
        $$$$$$$$$$$$$$$$$$$$$
        """   
        self.DifferentialFunction(Ynp1,self.Fnp1,time=self.Clock+self.TimeStep)
        self.DifferentialFunction(self.Yn,self.Fn,time=self.Clock)      
        self.Fnp1[:]=self.Fnp1[:]*self.TimeStep
        self.Flux_TMP[:]=Ynp1[:]-self.Yn[:]-self.Fnp1[:]

    def ImplicitRK2Flux(self,Ynp1):
        """
        
        Function to build implicit RK2 flux
        Write a flux vector to be equal to zero
        2th order
 
        Inputs : Ynp1[:] ; array of FLT : Variables vector at the computed position x+h
        Outputs : none
         
        $$$$$$$$$$$$$$$$$$$$$
        """      
        self.DifferentialFunction(Ynp1,self.Fnp1,time=self.Clock+self.TimeStep)
        self.DifferentialFunction(self.Yn,self.Fn,time=self.Clock)            
        self.Flux_TMP[:]=Ynp1-self.Yn[:]-self.TimeStep*(self.Fnp1[:]+self.Fn[:])*self.inv2

    def ImplicitAdam3o(self,Ynp1):
        """
        
        Function to build implicit Third order Adams-Moulton flux
        Write a flux vector to be equal to zero
        3th order
 
        Inputs : Ynp1[:] ; array of FLT : Variables vector at the computed position x+h
        Outputs : none
         
        $$$$$$$$$$$$$$$$$$$$$
        """  
        self.DifferentialFunction(Ynp1,self.Fnp1,time=self.Clock+self.TimeStep)
        self.DifferentialFunction(self.Yn,self.Fn,time=self.Clock)       
   

        self.Flux_TMP[:]=Ynp1[:]-self.Yn[:]-self.TimeStep*(self.Fnp1[:]*5.0*self.inv12+8.0*self.inv12*self.Fn[:]-1.0*self.inv12*self.Fnm1[:])
 
    def ImplicitAdam4o(self,Ynp1):
        """
        
        Function to build implicit Fourth order Adams-Moulton flux
        Write a flux vector to be equal to zero
        4th order

        Inputs : Ynp1[:] ; array of FLT : Variables vector at the computed position x+h
        Outputs : none
                
        $$$$$$$$$$$$$$$$$$$$$
        """      
        self.DifferentialFunction(Ynp1,self.Fnp1,time=self.Clock+self.TimeStep)
        
        self.DifferentialFunction(self.Yn,self.Fn,time=self.Clock)  
        
  
        self.Flux_TMP[:]=Ynp1[:]-self.Yn[:]-self.TimeStep*self.inv24*(self.Fnp1[:]*9.0+19.0*self.Fn[:]-5.0*self.Fnm1[:]+self.Fnm2[:])

    def ImplicitBDF2(self,Ynp1):
        """
        
        Function to build implicit Second order BDF (backward differentiation formula) flux
        Write a flux vector to be equal to zero
        2th order
 
        Inputs : Ynp1[:] ; array of FLT : Variables vector at the computed position x+h
        Outputs : none
         
        $$$$$$$$$$$$$$$$$$$$$
        """      
         
        self.DifferentialFunction(Ynp1,self.Fnp1,time=self.Clock+self.TimeStep)
        self.DifferentialFunction(self.Yn,self.Fn,time=self.Clock)  

        self.Fnp1[:]=self.Fnp1[:]*2.0*self.inv3*self.TimeStep
        self.Flux_TMP[:]=Ynp1[:]-4.0*self.inv3*self.Yn[:]+1.0*self.inv3*self.Ynm1[:]-self.Fnp1[:]
    
    def ImplicitBDF3(self,Ynp1):
        """
        
        Function to build implicit Third order BDF (backward differentiation formula) flux
        Write a flux vector to be equal to zero
        3th order
 
        Inputs : Ynp1[:] ; array of FLT : Variables vector at the computed position x+h
        Outputs : none
         
        $$$$$$$$$$$$$$$$$$$$$
        """ 
        self.DifferentialFunction(Ynp1,self.Fnp1,time=self.Clock+self.TimeStep)
        self.DifferentialFunction(self.Yn,self.Fn,time=self.Clock)  
        
        self.Fnp1[:]=self.Fnp1[:]*6.0*self.inv11*self.TimeStep
        self.Flux_TMP[:]=Ynp1[:]-18.0*self.inv11*self.Yn[:]+9.0*self.inv11*self.Ynm1[:]-2.0*self.inv11*self.Ynm2[:]-self.Fnp1[:]
        
    def BuildJacobianMatrix(self):
        """
        
        Function to build Jacobian matrix according to a flux function
        Every derivative is computed using two evaluation around the targetted position
        For a x position, fluxes are computed at
        (1.0+hf).x+ha and (1.0-hf).x-ha
        dh = 2.0.hf.x+2.0.ha
        hf and ha are two constants very small (1e-6 and 1e-5) specified in the initialization 
        of the class
 
        Inputs : none
        Outputs : none
         
        $$$$$$$$$$$$$$$$$$$$$
        """ 
        hf=self.hf
        ha=self.ha
        
        TMP_NumProb=copy.deepcopy(self)
        
        
        self.Ytmp[:]=self.Ynp1[:]
        for i in range(self.NbVariables):
            # Construction du dY
            dh=(2.0*hf)*self.Ytmp[i]+2.0*ha
            
            self.Ytmp[i]=((1.0-hf)*self.Ytmp[i]-ha)
            self.BuildFluxFunction(self.Ytmp)        
            self.Flux_m1[:]=self.Flux_TMP[:]
            
            self.Ytmp[i]=self.Ytmp[i]+dh
            self.BuildFluxFunction(self.Ytmp)
            self.Flux_p1[:]=self.Flux_TMP[:]
            inv_dY=1.0/dh
            self.JacobianMatrix[:,i]=(self.Flux_p1[:]-self.Flux_m1[:])*inv_dY
            self.Ytmp[i]=self.Ynp1[i]        

    def ImplicitIntegration(self):
        """
        
        Function to perform time integration using implicit approach
        For every time step, the flux vector is evaluated at a new position equal to 
        x_new=x_old - inv(J)*Flux*limitor
        where J is the Jacobian matrix of the flux vector and limitor is a constant 
        equal to 0.95 specified in the initialization function
 
        Inputs : none
        Outputs : none
        
        $$$$$$$$$$$$$$$$$$$$$
        """ 
    
    
        # Posons le vecteur Ynp1 etat a l'instant n+1
        # Posons le vecteur Yn etat a l'instant n
        # Posons Fn et Fnp1 les differentes contributions aux etats n et n+1
        # Le vecteur Flux est donc egal a Ynp1 - Yn - h(Fnp1) (pour Euler Implicit)
        
        converged = False
        TotError=1e30
        StopIterativeMathod=0
        k=0
        
        while( StopIterativeMathod == 0):
            
            
            self.BuildJacobianMatrix()
            
            Inv_JacobianMatrix=np.linalg.inv(self.JacobianMatrix)
            
            self.BuildFluxFunction(self.Ynp1)
            
            Delta=np.linalg.solve(Inv_JacobianMatrix,-self.Flux_TMP[:])
            
            self.Ynp1[:]=self.Ynp1[:]+Delta
            
            Tmp_Error=np.sum(abs(Delta))
            
            if (k == (self.MaxIte-1) or Tmp_Error<self.JacobianError):
                StopIterativeMathod=1  
    
            
            TotError=Tmp_Error    
            k=k+1

        restart = False
    
    def SolveDifferentialProblem(self):
        """
        
        Main funciton of the class which make the integration over the number of iteration by calling
        the selected scheme. 
        This function allows computing ODE with both implicit and explicit schemes
        
        This function is called into the InputFile.py once every parameter of the problem is well-defined.
        
        $$$$$$$$$$$$$$$$$$$$$
        """    
        StartTime=time.time()
        self.DefineInitialCondition()
        
        EndIfDuration=False
        EnfIfIteration=True
        EndIfStopDat=False
        for iter in range(self.NbIteration):
            
            if (iter%self.NbIterationSaved==0):
                self.RecordData()
                print('Iteration/TimeStep/Duration :'+str(iter)+', '+str(self.TimeStep)+', '+str(self.Clock))
            if (self.ExplicitBool):    
                if (self.AdaptativeTimeStep):
                    self.TimeIntegrationFunction()
                    self.dy_TestOrder[:]=self.dYtmp[:]
                    self.ExplicitFunctionForAdaptTimeStep()
                    self.AdaptTimeStep()
                    self.TimeIntegrationFunction()
                else:
                    self.TimeIntegrationFunction()
            if (self.ImplicitBool):        
                if (iter==0):
                    self.DefineDict()
                    if (self.SchemeName.find('ImplicitAdam')!=-1 or self.SchemeName.find('ImplicitBDF')!=-1):
                        self.TimeIntegrationFunction=self.SchemeFunctionsDic['ImplicitEuler']
                        self.BuildFluxFunction=self.FluxFunctionsForImplicit['ImplicitEuler'] 
                elif (iter==1):
                    self.DefineDict()
                    if (self.SchemeName.find('ImplicitAdam4o')!=-1 or self.SchemeName.find('ImplicitBDF3')!=-1):
                        self.TimeIntegrationFunction=self.SchemeFunctionsDic['ImplicitAdam3o']
                        self.BuildFluxFunction=self.FluxFunctionsForImplicit['ImplicitAdam3o'] 
                else:
                    self.DefineDict()
                        
                 
                self.TimeIntegrationFunction()
                    
                
                    
            
            self.UpdateTimeIntegration(self.TimeStep)
            
            if (self.Clock>self.Duration):
                self.RecordData()
                EndIfDuration=True
                EnfIfIteration=False
                EndIfStopDat=False
                break
            if (os.path.isfile('stop.dat')):
                self.RecordData()
                EndIfDuration=False
                EnfIfIteration=False
                EndIfStopDat=True
                break  
        EndTime=time.time()-StartTime
        TimePerIter=EndTime/iter
        print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        print('Problem solved')
        print('Iterations number = '+str(iter))
        print('CPU Time          = '+str(EndTime)+' seconds')
        print('CPU Time Per Iter = '+str(TimePerIter)+' seconds/iter')
        if (EndIfDuration):
            print('Computation stopped because of duration')
        if (EnfIfIteration):
            print('Computation stopped because of iteration number')
        if (EndIfStopDat):  
            print('Computation stopped because of stop.dat file')
        print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        self.FormateData()         
        
    def DisplayHelp(self):
        """     
        
        Function containing every call of help() methods for the funcitons
        Inputs  : nothing (class itself)
        Outputs : nothing (terminal output) 
        
        $$$$$$$$$$$$$$$$$$$$$
        """
        
        help(self.__init__)
        help(self.RecordData)
        help(self.FormateData)
        help(self.UpdateTimeIntegration)
        help(self.DefineInitialCondition)
        help(self.NextStepComputation)
        help(self.ExplicitEulerStep)
        help(self.ExplicitRK2)
        help(self.ExplicitRK4)
        help(self.ExplicitRK41)
        help(self.ExplicitRK45)
        help(self.ExplicitRK451)
        help(self.AdaptTimeStep)
        help(self.ImplicitEulerFlux)
        help(self.ImplicitRK2Flux)
        help(self.ImplicitAdam3o)
        help(self.ImplicitAdam4o)
        help(self.ImplicitBDF2)
        help(self.ImplicitBDF3)
        help(self.ImplicitIntegration)
        help(self.BuildJacobianMatrix)
        help()
        help(self.SolveDifferentialProblem)
        help(self.DisplayHelp)
