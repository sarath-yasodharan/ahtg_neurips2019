#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 10:02:24 2018

@author: sarath
computes error exponents in Neyman-Pearson formulation. 
type2errors from index 9 onwards prints the exponent.
"""

import numpy
import operator as op
from scipy.optimize import linprog
import matplotlib.pyplot as plt

# set cost function of the attacker in line 44, range of n in line 64
# set other parameters here:
qstar = 0.8 # minimum point of cost function
gamma = 1 # scaling factor in the defender's revenue expression 
d = 2 # number of elements in the alphabet set {1,2,\cdots, d}
p = 0.5 #H0
q1=0.7 #Leftmost point of H1
q2=0.8 #Rightmost point of H1
#n = 50 # number of observations
epsilon = 0.1 #bound on false alarm probability
k = 100; # number of types for discretisation

def prob(no_ones,p):
    prob = (p**no_ones)*((1-p)**(n-no_ones))
    return prob

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, xrange(n, n-r, -1), 1)
    denom = reduce(op.mul, xrange(1, r+1), 1)
    return numer//denom

# given q and decision rule, compute the utility of zero-sum equivalent game
# if no of ones \geq threshold, then accept H1
def u(q,rejection_threshold): 
    temp = 0
    for i in range (0,rejection_threshold):
        comb = ncr(n,i)
        temp = temp + comb*prob(i,q)
    temp = temp - (0.01*abs(q-qstar))**2
    return temp
def type2error(q,rejection_threshold):
    temp = 0
    for i in range (0,rejection_threshold):
        comb = ncr(n,i)
        temp = temp + comb*prob(i,q)
    return temp
    

def false_alarm_prob(rejection_threshold):
    temp = 0
    for i in range(rejection_threshold,n+1):
        comb = ncr(n,i)
        temp = temp + gamma*comb*prob(i,p)
    return temp

def relative_entropy(p,q):
    return p*numpy.log(p/q)+(1-p)*numpy.log((1-p)/(1-q))

iterations=100
type2errors = [0]*(iterations)
attacker_equilibrium = [0]*(iterations)
for n in range(10,iterations+1):
    print n
    #A contains the values of u for all possible q in grid Q and index-based decision rule
    A = numpy.zeros((n+3,k+1))


    for i in range(0,n+2):
        A[i,0] = 1
        for j in range(0,k):
            q = q1+j*(q2-q1)/k
            A[i,j+1] = -1*u(q,i)
            for j in range(0,k):
                A[n+2,j+1] = 1

    # compute false alarm probabilities for index-based decision rules
    false_alarm = [0]*(n+2)
    threshold_star = 0 # to indicate the possible set of decision rules that meets
                    # the false alarm constraint
    for i in range(0,n+2):
        false_alarm[i] = false_alarm_prob(i)

    for i in range(0,n+2):
        if false_alarm[i] <= epsilon:
            threshold_star = i
            break

    A_ub = A[threshold_star:n+3,:] # truncate the matrix A for LP
    c = [0]*(k+1)
    c[0] = -1 # linprog does minimisation
    b_ub = numpy.zeros((A_ub.shape[0],1))
    b_ub[A_ub.shape[0]-1]=1
    A_eq = numpy.zeros((1,k+1))
    for i in range(0,k):
        A_eq[0,i+1] = 1
        b_eq = 1

    soln = linprog(c,A_ub,b_ub,A_eq,b_eq,method='simplex')
    i_n = numpy.where(soln.x)[0][1]-1
    q_n = q1+i_n*(q2-q1)/k
    attacker_equilibrium[n-1] = q_n
    type2errors[n-1] = -1*numpy.log(type2error(q_n,threshold_star))/n
    


