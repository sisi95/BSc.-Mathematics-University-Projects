# -*- coding: utf-8 -*-
"""
Created on Sat Nov  5 17:10:05 2016

@author: sjoukjeijlstra
"""

from numpy import *
from pylab import *
import numpy.linalg
from numpy import array, zeros, diag, diagflat, dot
import math

def f(x):
    f = (1-x**2)*sin(2*pi*x)
    return f
    
def g(x):
    g = sin(x)
    return g

x = linspace(-1.5, 2.5)
plot(x, f(x))

# The data points are equally spaced on [-1, 2] for polynomials degree n = 6 and n = 10

h = []  
for i in xrange(0, 7, 1):
    a = f(-1 + 0.5*i)
    h.append([-1+ (0.5*i), round(a, 3)])
plot(*zip(*h), marker='o', color='r', ls='')
xlim([-1.5,2.5])
#show()
    
j = []
for i in xrange(0, 11, 1):
    a = f(-1 + 0.3*i)
    j.append([round(-1+ (0.3*i), 3), round(a, 3)])
plot(*zip(*j), marker='o', color='b', ls='')
xlim([-1.5,2.5])
show()

# Definition for interpolant

def p(x):
    w1 = []
    #w2 = []
    
    for i in xrange(0, 6):
        p = ((x - h[i+1][0])/(h[i][0] - h[i+1][0]))*f(h[i][0]) + ((x - h[i][0])/(h[i+1][0] - h[i][0]))*f(h[i+1][0])
        #x = linspace(h[i][0], h[i+1][0])
        w1.append( p)
    #for i in xrange(0, 10):
      #  p = ((x - j[i+1][0])/(j[i][0] - j[i+1][0]))*f(j[i][0]) + ((x - j[i][0])/(j[i+1][0] - j[i][0]))*f(j[i+1][0])
      #  w2.append(p)
    plot(w1)
    #plot(w2)
    #xlim([-10, 2])
    show()
