# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 10:32:58 2016

@author: sjoukjeijlstra
"""

from numpy import *
from pylab import *
import numpy.linalg
from numpy import array, zeros, diag, diagflat, dot

# Defining the given matrix 

a = array([[9.0, 1.0, 2.0, 1.0, 0.0, 4.0], [0.0, 8.0, 1.0, 2.0, 1.0, 0.0], [1.0, 0.0, 9.0, 1.0, 2.0, 1.0], [0.0, 1.0, 0.0, 7.0, 1.0, 2.0], [0.0, 0.0, 1.0, 0.0, 9.0, 1.0], [4.0, 0.0, 0.0, 1.0, 0.0, 6.0]])
b = array([27.0, 22.0, 17.0, 23.0, 13.0, 24.0])

print 'Question 1:'

# Initial guess for x0

x0 = array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

# Defining the Jacobi Method

def jacobi_method(x0, n):
    
    if x0 is None: 
        x0 = zeros(len(a[0]))
    
    # Create Tj = D^(-1)*(L+U)
    d = diag(a) # Diagonal of a in its own 1D array
    l = a - diagflat(d) # Creates L+U with main diagonal 0's. Makes diagonal of a into a 2D array (same type as a) then its subtracted from a

    for i in range(0, n):
        x0 = (b - dot(l, x0))/d # Iterates Jacobi by setting each row, x1, x2, ..., xn, to scheme we want. Then does it n times. 
    return x0

# Define the Gauss-Siedel Method

def gs_method(x0, n):
    
    if x0 is None:
        x0 = zeros(len(a[0]))
    
    # Create Tsg = (D-L)^(-1)*U
    l = tril(a, 0) # Lower triangular matrix of a with elements above diagonal zero
    u = -triu(a, 1) # Upper triangular matrix of a with elements under 1st diagonal zero
    
    for i in range(0, n):
        x0 = dot(inv(l), b + dot(u, x0)) # Using defintion of Gauss-Siedel method
    return x0

# Define the Successive Relaxation method

def sr_method(x0, w, n):
    
    if x0 is None:
        x0 = zeros(len(a[0]))
        
    l = tril(a, 0) # Lower triangular matrix of a with elements above -1st diagonal zero
    u = -triu(a, 1) # Upper triangular matrix of a with elements under 1st diagonal zero
    d = diag(a) # Diagonal of a 
    
    for i in range(0, n):
        x0 = dot(inv((w*l)), dot(((1-w)*d + u), x0) + w*b)
    return x0

print 'Question 2:'

#print sr_method(x0, -0.5, 10) # w<0 does not converge 
                              # w = 0 gives a singular matrix 
#print sr_method(x0, 0.1, 10)  # 0<w<1 gives successive under relaxation for systems convergent if GS method doesn't converge 
#print sr_method(x0, 1.0, 10)  # w = 1 is the Gauss-Siedel method
#print sr_method(x0, 1.5, 10)  # 1<w<2 gives successive over relaxation method for systems convergent using GS method
#print sr_method(x0, 2.0, 10)  # w>=2 does not converge

print 'Question 3:' 

# The spectral radius of a matrix is the sum of its modulus eigenvalues

# Spectral radius of Tj

def specradius_Tj():
    
    d = diag(a)
    l = -tril(a, -1)
    u = -a + tril(a, 0)
    
    # Create iteration matrix for Tj: inv(D)*(L+U)
    iteration_matrix_Tj = dot(inv(diagflat(d)), l + u)
    
    inf_norm_Tj = amax(abs(iteration_matrix_Tj)) # Infinity norm
    
    spec_radius_Tj = max(eigvals(iteration_matrix_Tj)) # Spectral radius of Tj
    return real(spec_radius_Tj) # Returns real part of spectral radius
    
print 'p(Tj) =', specradius_Tj()

# Spectral radius of Tgs

def specradius_Tgs():
    
    d = diag(a) 
    l = -tril(a, -1)
    u = -a + tril(a, 0)
    
    # Create iteration matrix for Tgs: inv(D-L)*U 
    iteration_matrix_Tgs = dot(inv(d-l), u)
    
    inf_norm_Tgs = amax(abs(iteration_matrix_Tgs)) # Infinity norm 

    spec_radius_Tgs = max(eigvals(iteration_matrix_Tgs))
    return round(real(spec_radius_Tgs), 8) # Real part of spectral radius to 8 dp
    
print 'p(Tgs) =', specradius_Tgs()

# Spectral radius of Tsor

def specradius_Tsor(w):
    
    d = diag(a) 
    l = -tril(a, -1)
    u = -a + tril(a, 0)
    
    # Create iteration matrix for Tsor: inv(D-wL)*[(1-w)*D + U]
    iteration_matrix_Tsor = dot(inv(d-w*l), dot((1-w), d) + u)
    
    inf_norm_Tsor = amax(abs(iteration_matrix_Tsor)) # Infinity matrix for Tsor
    
    spec_radius_Tsor = max(eigvals(iteration_matrix_Tsor))
    return round(real(spec_radius_Tsor), 8)
    
print 'p(Tsor) =', specradius_Tsor(1.1)

