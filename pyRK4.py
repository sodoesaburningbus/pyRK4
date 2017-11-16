#This module implements the Runge-Kutta 4th order numerical
#integration scheme in Python.
#The sole function in this module is pyRK4.
#To use this module, the user will need to loop over their variable
#of integration. An example is provided below using position and velocity:
#
#-------------Example--------------
#
# from pyRK4 import pyRK4
# import numpy
# 
# #Defining the bounds of integration
# t0 = 0
# tf = 10
# dt = 1
#
# #Initial position
# x = 0
#
# #Function of the DE to solve.
# def dxdt(x,t): #Even if x is not used, pyRK4 expects two arguments
#   return t**2
#   
#
# #Integrating position
# for t in numpy.arange(t0,tf,dt) #The upper bound is supposed to be open.
#   x = pyRK4(dxdt,x,t,dt)
#
# print("Final position is: {}".format(x))
#
#----------------------------------
#
#Written by Chris Phillips


#Defining the 4th Order Runge-Kutta function
#This function returns the value of y(i+1) when given y(i)
#where i is the current iteration of the integration loop.
#Inputs:
#dydx, this is the function handle that corresponds to dy/dx
#y, this is the current value of the variable being integrated
#x, this is the current value of the variable that y is being integrated over
#dx, this is the integration step size
#
#Outputs:
#y(i+1), the next value of y in the integration.
def pyRK4(dydx,y,x,dx):
    k1 = dydx(y,x)*dx
    k2 = dydx(y+0.5*k1,x+0.5*dx)*dx
    k3 = dydx(y+0.5*k2,x+0.5*dx)*dx
    k4 = dydx(y+k3,x+dx)*dx
    return y+(k1+2.0*k2+2.0*k3+k4)/6.0
