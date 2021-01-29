# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 09:24:57 2015

@author: 140001917
"""

import sympy
from numpy import *
from pylab import *

print 'Q1: on paper'
print 'Q2 pt. 1: '

# Predefined variables
a1 = a2 = b1 = b2= 1
c1 = d1 = 0.1
c2 = 0.2 
d2 = 0.3

# Define the functions for which we require the roots.
def g(x,y):
    g = a1-b1*x-c1*exp(d1*y)
    return g
    
def h(x,y):
    h = a2-b2*y-c2*exp(d2*x)
    return h
    
# Define the derivatives of g(x,y) and h(x,y)
def dgdx(x,y):
    dgdx = -b1
    return dgdx
    
def dgdy(x,y):
    dgdy = -c1*d1*y*exp(d1*y)
    return dgdy
    
def dhdx(x,y): 
    dhdx = -c2*d2*x*exp(d2*x)
    return dhdx
    
def dhdy(x,y): 
    dhdy = -b2
    return dhdy

def root(x,y,err,nmax):
    errorx = 1.0
    errory = 1.0
    n = 0
    while (((abs(errorx) > abs(err)) or (abs(errory) > abs(err))) and (n < nmax)):
        n = n+1
        errorx = (dgdy(x,y)*h(x,y)-g(x,y)*dhdy(x,y))/(dhdy(x,y)*dgdx(x,y)-dgdy(x,y)*dhdx(x,y)) 
        errory = (dgdx(x,y)*h(x,y)-g(x,y)*dhdx(x,y))/(dgdy(x,y)*dhdx(x,y)-dgdx(x,y)*dhdy(x,y))
        x = x + errorx
        y = y + errory
        return round(x, 6), round(y,6)
      
# Initial estimates
x=0.95
y=0.75
# Accuracy of the estimate
err = 1.0e-6
# Maximum number of iterations in case it does not converge.
nmax = 30
root = root(x,y,err,nmax)

print 'The root of g(x,y) = h(x,y) is', root

print 'Q2 pt. 2: on paper'
print 'Q2 pt. 3: '

def f(x,y):
    f = 1.0
    return f

n = 1000
step = 0.892/float(n)
x = zeros(n+1,dtype=float)
Integrand = zeros(n+1,dtype=float)
x[0] = 0.0
x[n] = 0.892
for i in range(0, n):
    x[i] = 0.0 + float(i)*step
    ylim1 = (1-0.2*exp(0.3*x[i]))
    ylim2 = (10*log(10*(1-x[i])))
    y = zeros(n+1,dtype=float)
    fvalue = zeros(n+1,dtype=float)
    y[0] = ylim1 
    y[n] = ylim2
    fvalue[0] = f(x[i],y[0])
    fvalue[n] = f(x[i],y[n])
    ystep = (ylim2)/float(n)
    for j in range(1,n):
        y[j] = j*ystep
        fvalue[j] = f(x[i],y[j])       
    Integrand[i] = trapz(fvalue,y,ystep)

Ivalue = trapz(Integrand,x,step)

print 'The integral is',round(Ivalue, 6),'for', n,'intervals.'

print 'Q3 pt. 1: code'

# Defining derivatives of f 
def f(x,y):
    f = (x**2-y**2+4*x+6*y+1)/(x**2+y**2+6)
    return f
    
def fx(x,y):
    fx = -4*((x**2)-x*((y**2)-3*y+2.5)-y**2-6)/((x**2+y**2+6)**2)
    return fx
    
def fy(x,y):
    fy = -6*((y**2)+(2/3)*y*((x**2)+2*x+3.5)-(x**2)-6)/((y**2+x**2+6)**2)
    return fy
    
def fxx(x,y):
    fxx = 8*(x**3-1.5*x**2*(y**2-3*y+2.5)-3*x*(y**2+6)+0.5*(y**2+6)*(y**2-3*y+2.5))/((x**2+y**2+6)**3)
    return fxx

def fyy(x,y):
    fyy = (12*(y**3+y**2*(x**2+2*x+3.5)-3*y*(x**2+6)-(1/3)*(x**2+6)*(x**2+2*x+3.5)))/((y**2+x**2+6)**3)
    return fyy

def fxy(x,y):
    fxy = 8*((x**3)*(y-1.5)+(3*(x**2)*y)-x*((y**3)-4.5*(y**2)-y+9)-y*((y**2)+6))/((y**2+x**2+6)**3)
    return fxy

# Defining grad (delta)
def triangle(x, y):
    triangle = fxx(x,y)*fyy(x,y) - (fxy(x,y)**2)
    return round(triangle, 6) 
    
# Nature of stationary points based on conditions
def determine_nature(x,y):
   
    if triangle(x,y) > 0.0 and fxx(x,y) > 0.0:
        print 'minimum.'
    if triangle(x,y) > 0.0 and fxx(x,y) < 0.0:
        print 'maximum.'
    if triangle(x,y) < 0.0: 
        print 'saddlepoint.'
    if triangle(x,y) == 0.0:
        print 'No information is given for this point.'
        
print 'Q3 pt. 2: '

# Creating a contour graph of F
xp = arange(-5.0, 5.01, 0.01)
yp = arange(-5.0, 5.01, 0.01)

X,Y = meshgrid(xp,yp)

F = (X**2 - Y**2 +4*X +6*Y +1)/(X**2 +Y**2 +6)
fig = figure(1)

fig = figure(1)
ax = fig.add_subplot(111)
ax.contour(X,Y,F, levels = [-5.0, -4.9, -4.8, -4.7, -4.6, -4.5, -4.4, -4.3, -4.2, -4.1, -4.0, -3.9, -3.8, -3.7, -3.6, -3.5, -3.4, -3.3, -3.2, -3.1, -3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0, 3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0], cmap = 'hsv')

ax.set_title('$F = f(x,y) = (x^2 - y^2 +4x +6y +1)/(x^2 +y^2 +6)$', size = 'large')
ax.set_xlabel('$x$', size = 'large')
ax.set_ylabel('$y$', size = 'large')

show()

print 'Looking at the graph, I see that the function centres around approx. (-0.5, -3.5) as a minimum, (-2.5, 2.5) as a saddlepoint and (3.5, 1) as a maximum.'

print 'Q3 pt. 3: '

# Choose estimates for stationary points based on contour graph
x3 = -0.5
y3 = -3.5

x4 = 3.5
y4 = 1

x5 = -2.5
y5 = 2.5 

# Find roots of f based on 3 estimated pairs of coordinates
def froot(x,y,err,nmax):
    errorx = 1.0
    errory = 1.0
    n = 0
    while (((abs(errorx) > abs(err)) or (abs(errory) > abs(err))) and n < nmax):
        n = n+1
        errorx = (fxy(x,y)*fy(x,y)-fx(x,y)*fyy(x,y))/(fyy(x,y)*fxx(x,y)-fxy(x,y)*fxy(x,y)) 
        errory = (fxx(x,y)*fy(x,y)-fx(x,y)*fxy(x,y))/(fxy(x,y)*fxy(x,y)-fxx(x,y)*fyy(x,y)) 
        x = x + errorx
        y = y + errory
        return round(x,6), round(y,6)

print 'The location of the first stationary point is', froot(x3, y3, err, nmax),'and is a', determine_nature(x3,y3)
print 'The location of the second stationary point is', froot(x4, y4, err, nmax),'and is a', determine_nature(x4,y4)
print 'The location of the third stationary point is', froot(x5, y5, err, nmax),'and is a', determine_nature(x5,y5)

# The exact locations of stationary points are wrong, I think it is because of the froot() errorx and errory
# However, I know that these errors in froot() are exactly the same as in Question 2, which is why I cannot solve this problem
# My predictions of what nature the points had was confirmed by determine_nature()

