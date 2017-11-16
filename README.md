# pyRK4
This module implements the 4th order Runge-Kutta method in Python.

To install simply clone the repository or download the source file into its own dirctory and run the setup script with "python3 setup.py install"

The module contains two functions, "pyRK4_step" and "pyRK4_integrate"
The step function implements a single iteration of the RK4 method.
The integrate function does the full integration when provided the bounds of integration.

This module is intended for use with Python 3 but should work in Python 2 with no changes other than installing with python2 instead of python3.

This software is written by Christopher Phillips and implies no gurantees of functionality.
The software may be freely distributed provided this Readme remains unmodified and is distributed with the software.

Requirements:
Python
Numpy

Examples for using the code are below.

#-------------Example 1 pyRK4_step --------------

 from pyRK4 import pyRK4
 import numpy
 
 #Defining the bounds of integration
 t0 = 0
 tf = 10
 dt = 1

 #Initial position
 x = 0

 #Function of the DE to solve.
 def dxdt(x,t): #Even if x is not used, pyRK4 expects two arguments
   return t**2
   

 #Integrating position
 for t in numpy.arange(t0,tf,dt) #The upper bound is supposed to be open.
   x = pyRK4(dxdt,x,t,dt)

 print("Final position is: {}".format(x))

#------------------------------------

#-------------Example 2 pyRK4_integrate--------------

 from pyRK4 import pyRK4_integrate
 
 #Defining the bounds of integration
 x0 = 0
 xf = 100
 dx = 0.01

 #initial condition
 y0 = 0

 #The DE to solve
 def dydx(y,x): #Even if y is not used, pyRK4 expects two arguments
   return x**2

 #Integrate
 [ys,xs] = pyRK4_integrate(dydx,y0,x0,xf,dx)

 print("Final value is: {}".format(ys[-]))

#------------------------------------
