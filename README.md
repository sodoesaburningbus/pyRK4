# pyRK4
This module implements the 4th order Runge-Kutta method in Python.

To install simply clone the repository or download the source file into its own dirctory and run the setup script with "python3 setup.py install" or

This module is intended for use with Python 3 but should work in Python 2 with no changes.
An example of using the code is below.

This software is written by Christopher Phillips and implies no gurantees of functionality.
The software may be freely distributed provided this Readme remains unmodified and is distributed with the software.


-------------Example--------------

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

#----------------------------------
