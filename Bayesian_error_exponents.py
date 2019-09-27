#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 08:47:38 2018

@author: sarath

comment: plot error exponents in the Bayesian formulation.
values of the exponents is stored in the array error_exponents, it will be printed 
when you run this code.
"""

import numpy
import operator as op
from scipy.optimize import linprog
import matplotlib
import matplotlib.pyplot as plt

# set cost function in line 79, and range of n in line 64,
# set other parameters below:
gamma = 1
qstar = 0.8 # minimum point of cost function
d = 2 # number of elements in the alphabet set {1,2,\cdots, d}
p = 0.5 #H0
q1=0.7 #Leftmost point of H1
q2=0.9 #Rightmost point of H1
k = 100 #number of discrete points in the set Q

# compute probability of a sequence given probability of ones
def prob(n,no_ones,p):
    return (p**no_ones)*((1-p)**(n-no_ones))

def my_range(start, end, step):
    while start <= end:
        yield start
        start += step

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, xrange(n, n-r, -1), 1)
    denom = reduce(op.mul, xrange(1, r+1), 1)
    return numer//denom

# given q and decision rule, compute the utility of zero-sum equivalent game
# if no of ones \geq threshold, then accept H1. istar is the threshold.
def type2error(n,q,istar): 
    temp = 0
    for i in range (0,istar):
        comb = ncr(n,i)
        temp = temp + comb*prob(n,i,q)
#    temp = temp - (abs(q-qstar))
    return temp
def type1error(n,istar):
    if istar == n+1:
        return 0
    else:
        temp = 0
        for i in range (istar,n+1):
            temp = temp + ncr(n,i)*prob(n,i,p)
    return temp

def relative_entropy(p,q):
    return p*numpy.log(p/q)+(1-p)*numpy.log((1-p)/(1-q))

error_exponent = [0]*(50)
l = 0
for n in range(100,110,10):    
    #print n
    # solve the attacker's LP
    A_ub = numpy.zeros((n+3,k+1))
    b_ub = numpy.zeros((n+3,1))
    error_matrix = numpy.zeros((n+3,k+1))
    c = [0]*(k+1)
    B_ublb = numpy.zeros((k+1,2))
    B_ublb[0,0]= -1
    B_ublb[:,1] = 100

    for i in range(0,n+2):
        #print i
        A_ub[i,0] = 1
        for j in range(0,k):
            q = q1+j*(q2-q1)/k
            A_ub[i,j+1] = -1*(type2error(n,q,i)+gamma*type1error(n,i)-(abs(q-qstar)))
            error_matrix[i,j+1] = (type2error(n,q,i)+gamma*type1error(n,i))
    for j in range(1,k+1):
        A_ub[n+2,j] = 1
    b_ub[n+2] = 1
    c[0] = -1 # becasue LP function does minimisation
    A_eq = numpy.ones((1,k+1))
    A_eq[0,0] = 0
    b_eq = 1
    soln_attackerLP = linprog(c,A_ub,b_ub,A_eq,b_eq, method='simplex')
    #print soln_attackerLP
    
    A_ub_def = numpy.zeros((k+1,n+3))
    b_ub_def = numpy.zeros((k+1,1))
    error_matrix_transpose = numpy.zeros((k+1,n+3))
    b_ub_def[k,0] = -1
    A_ub_def[0:k,1:n+3] = -1*numpy.transpose(A_ub[0:n+2,1:k+1])
    error_matrix_transpose[0:k,1:n+3] = numpy.transpose(error_matrix[0:n+2,1:k+1])
    for i in range(0,k):
        A_ub_def[i,0] = -1
    for j in range(1,n+3):
        A_ub_def[k,j] = -1
    c_def = [0]*(n+3)
    c_def[0] =1
    A_eq_def = numpy.ones((1,n+3))
    A_eq_def[0,0] = 0
    b_eq_def = 1
        
    B_ublb_def = numpy.zeros((n+3,2))
    B_ublb_def[0,0]= -1
    B_ublb_def[:,1] = 100
    soln_defLP = linprog(c_def,A_ub_def,b_ub_def,A_eq_def,b_eq_def,method = 'simplex')
    #print soln_defLP
    
    x = soln_attackerLP.x
    x = numpy.asmatrix(x[1:k+1])
    y = soln_defLP.x
    y = numpy.asmatrix(y[1:n+3])
    error_exponent [l] = -1*numpy.log(numpy.matrix.sum(numpy.multiply(numpy.matmul(numpy.transpose(x),y), error_matrix_transpose[0:k,1:n+3])))/n
    print error_exponent[l]
    l = l+1
#find actual error exponent analytically
array1 = [0]*(100)
array2 = [0]*(100)
for i in range (0,100):
    q = 0.5+i*(qstar-0.5)/100
    array1[i] = relative_entropy(q,0.5)
    array2[i] = relative_entropy(q,qstar)
    if abs(array1[i] - array2[i])<0.005:
        istar = i
#        print istar --this will print the index and fetch array1 or array2 
#        for error exponent
        
