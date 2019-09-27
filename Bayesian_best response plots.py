#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 13:40:42 2018

@author: sarath

comment: generate best response plots for Bayesnain formulation in d=2.
"""

import numpy
import math
import matplotlib.pyplot as plt
import operator as op
from decimal import Decimal

d = 2 # number of elements in the alphabet set {1,2,\cdots, d}

# set parameters here. cost function needs to be set on line 
qstar = 0.8
gamma = 1 # scaling factor in the defender's revenue expression 
p = 0.5 #H_0
q1=0.7 #Leftmost point of H1
q2=0.9 #Rightmost point of H1
n = 200 # number of observations

# compute probability of a sequence given probability of ones
def prob(no_ones,p):
    return (p**no_ones)*((1-p)**(n-no_ones))

# likelihood ratio test, given q and an n-length sequence
def lrtest(no_ones,q):
    threshold = ((numpy.log((1-p)/(1-q)))*n+numpy.log(gamma))/numpy.log((q*(1-p))/(p*(1-q)))
    if no_ones >= threshold:
       return 1
    else:
       return 0

# compute cCr       
def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, xrange(n, n-r, -1), 1)
    denom = reduce(op.mul, xrange(1, r+1), 1)
    return numer//denom

threshold = [0]*(100)
for i in range(0,100):
	q = q1+i*(q2-q1)/100
	threshold [i]= int(math.ceil(((numpy.log((1-p)/(1-q)))*n+numpy.log(gamma))/numpy.log((q*(1-p))/(p*(1-q)))))
threshold_min = threshold[0]
threshold_max = threshold[99]
#print threshold

#utility of the attacker for a given q, when the defender plays a LRT with threshold 'threshold'	
def ua_given_threshold (q,threshold):
	temp = 0
	j = 0
	for i in range(0,threshold):
		comb = Decimal(ncr(n, j))
		j = j+1
		temp = temp + Decimal(prob(i,q))*comb
	temp = temp -Decimal(abs(q-qstar))
	return temp

# best response of the attacker for a given threshold
def attacker_best_response (threshold):
	output = [0]*(20)
	#print threshold
	for i in range (0,20):
		q = q1+i*(q2-q1)/20
		output[i] = ua_given_threshold(q,threshold)
	istar = numpy.argmax(output)
	q = q1+istar*(q2-q1)/20
	return q

no_thresholds = 21
output = [0]*(no_thresholds)
q = qstar
threshold_min = int(math.ceil(((numpy.log((1-p)/(1-q)))*n+numpy.log(gamma))/numpy.log((q*(1-p))/(p*(1-q))))) - 10
print threshold_min
threshold_max = threshold_min+20
threshold_running = threshold_min
for i in range (0,no_thresholds):
	output[i]=attacker_best_response(threshold_running)
	threshold_running = threshold_running + 1
x = numpy.linspace(q1,0.898,100)
y = output
plt.gcf().subplots_adjust(bottom=0.15)
#plt.gcf().subplots_adjust(left=0.20)
plt.scatter(numpy.linspace(q1,0.898,100),threshold,marker='o')
plt.scatter(y,numpy.linspace(threshold_min,threshold_max,no_thresholds),c='r',marker='x')
plt.xlabel('q',fontsize=14)
plt.ylabel('Threshold', fontsize=14)
plt.savefig('test3.png')
plt.show()