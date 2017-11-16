#This module implements the Runge-Kutta 4th order numerical
#integration scheme in Python.
#The pyRK4_step function in this module implements the basic RK4 method.
#It is deisgned to be incorporated into a loop defined by the user and
#step through the integration.
#Example 1 shows how this function might be used.
#
#The pyRK4_integrate function is a stand-alone function that does the full
#integration in a single function call.
#Example 2 illustrates how to use this function.
#
#Module requirements:
#numpy
#
#-------------Example 1--------------
#
# from pyRK4 import pyRK4_step
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
#   x = pyRK4_step(dxdt,x,t,dt)
#
# print("Final position is: {}".format(x))
#
#------------------------------------
#
#-------------Example 2--------------
#
# from pyRK4 import pyRK4_integrate
# 
# #Defining the bounds of integration
# x0 = 0
# xf = 100
# dx = 0.01
#
# #initial condition
# y0 = 0
#
# #The DE to solve
# def dydx(y,x): #Even if y is not used, pyRK4 expects two arguments
#   return x**2
#
# #Integrate
# [ys,xs] = pyRK4_integrate(dydx,y0,x0,xf,dx)
#
# print("Final value is: {}".format(ys[-]))
#
#------------------------------------
#
#This software may be freely shared provided that it is accompanied by
#its README and this header.
#
#Written by Christopher Phillips


#-------------BEGIN MODULE-------------

#Importing required modules
import numpy

#pyRK4_step
#Defining the basic 4th Order Runge-Kutta function
#This function returns the value of y(i+1) when given y(i)
#where i is the current iteration of the integration loop.
#Inputs:
#dydx, this is the function handle of dy/dx (the function to integrate.)
#y, this is the current value of the variable being solved for.
#x, this is the current value of the variable that y is being integrated over
#dx, this is the integration step size
#
#Outputs:
#y(i+1), the next value of y in the integration.
def pyRK4_step(dydx,y,x,dx):
    #The RK4 method
    k1 = dydx(y,x)*dx
    k2 = dydx(y+0.5*k1,x+0.5*dx)*dx
    k3 = dydx(y+0.5*k2,x+0.5*dx)*dx
    k4 = dydx(y+k3,x+dx)*dx
    return y+(k1+2.0*k2+2.0*k3+k4)/6.0

#-------------------------------------------------------------

#pyRK4_integrate
#Defining the full integration function using RK4
#This function returns a ist of lists containing the values of y and x at
#each step of the integration.
#Inputs:
#dydx, this is the function handle of dy/dx (the function to integrate.)
#y0, this is the initial value of y, the variable being solved for.
#x0, this is the initial value of x, the lower bound of integration.
#xf, this is the final value of x, the upper bound of integration.
#dx, this is the integration step size.
#Outputs:
#[y,x_steps]
#   y is a list containing the y values at each step of the integration
#   x_steps is a list of the x points that match each y
def pyRK4_integrate(dydx,y0,x0,xf,dx):
    #Computing the value of x at each step of the integration
    x_steps = numpy.arange(x0,xf,dx)  
    
    #Integrating
    y = [] #List o hold y values during integration
    y.append(y0) #Starting with the initial y.
    for x in x_steps:
        k1 = dydx(y[-1],x)*dx
        k2 = dydx(y[-1]+0.5*k1,x+0.5*dx)*dx
        k3 = dydx(y[-1]+0.5*k2,x+0.5*dx)*dx
        k4 = dydx(y[-1]+k3,x+dx)*dx
        y.append(y[-1]+(k1+2.0*k2+2.0*k3+k4)/6.0)
    
    #Appending final value to x_steps, so that each x matches the proper y.
    x_steps = numpy.append(x_steps,x_steps[-1]+dx)

    #Returning values
    return [y,x_steps]
