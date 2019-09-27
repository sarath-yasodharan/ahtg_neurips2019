#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 13:46:20 2018

@author: sarath
plots the utility of the attacker in a finer scale when defender plays a 
fixed threshold strategy
"""

import numpy
import math
import matplotlib.pyplot as plt
import operator as op
from decimal import Decimal

# set cost function of the attacker in line 57, set the x and y limits of plots
# in line 85 and 86
# set other parameters here:
qstar = 0.8
def_br = 133 # set best response threshold of the defender
gamma = 1 # scaling factor in the defender's revenue expression 
d = 2 # number of elements in the alphabet set {1,2,\cdots, d}
p = 0.5 #H_0
q1=0.7 #Leftmost point of H1
q2=0.9 #Rightmost point of H1
n = 200 # number of observations

# compute probability of a sequence given probability of ones
def prob(no_ones,p):
    prob = (p**no_ones)*((1-p)**(n-no_ones))
    return prob

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

#attacker best response for a given threshold -- very fine search
def attacker_best_response_finer (threshold):
	output = [0]*(1001)
	#print threshold
	for i in range (0,1001):
		q = qstar - 10/(2000*float(n)) +i*10*2/(2000*1000*float(n))
		output[i] = ua_given_threshold(q,threshold)
	return output


att_revenue = attacker_best_response_finer(def_br)
x = numpy.linspace(qstar - 10/(2000*float(n)),qstar - 10/(2000*float(n)) +1000*10*2/(2000*1000*float(n)),1001)
#print(x,att_revenue)
istar = numpy.argmax(att_revenue)
q =  qstar - 10/(2000*float(n)) +istar*10*2/(2000*1000*float(n))
print istar
print q
#plt.scatter(x,att_revenue)
#plt.show()
#plt.ylim(-1E-6,1E-6)
#
#plt.xticks(['1 2 3 4 5'])
#plt.plot(att_revenue)

plt.ylim(-2E-6,1E-5)
plt.xlim(0.79999375, 0.80000625)
plt.scatter(x,att_revenue)
plt.xticks(rotation=90)
plt.xlabel('q',fontsize=14)
plt.ylabel('Attacker revenue',fontsize=14)
plt.gcf().subplots_adjust(bottom=0.30)
plt.gcf().subplots_adjust(left=0.20)
plt.savefig('test6.png')