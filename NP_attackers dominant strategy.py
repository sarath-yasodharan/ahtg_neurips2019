#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 10:15:51 2018

@author: sarath
plots the attacker's threshold rule as a function of n (illustrates lemma A3 in the paper)
"""

import numpy
import operator as op
import matplotlib
import matplotlib.pyplot as plt

d = 2 # number of elements in the alphabet set {1,2,\cdots, d}

# set the cost function of the attacker in line 75 and range of n in line 60

# set the other parameters here:
qstar = 0.8 # minimum point of cost function
p = 0.5 #H0
q1=0.7 #Leftmost point of H1
q2=0.8 #Rightmost point of H1
epsilon = 0.1 #bound on false alarm probability

def my_range(start, end, step):
    while start <= end:
        yield start
        start += step

def prob(n,no_ones,p):
    return (p**no_ones)*((1-p)**(n-no_ones))
    

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, xrange(n, n-r, -1), 1)
    denom = reduce(op.mul, xrange(1, r+1), 1)
    return numer//denom

# given q and decision rule, compute the utility of zero-sum equivalent game
# if no of ones \geq threshold, then accept H1
def type2error(n,q,istar,pi): 
    temp = 0
    for i in range (0,istar+1):
        comb = ncr(n,i)
        temp = temp + comb*prob(n,i,q)
    temp = temp + ncr(n,istar+1)*prob(n,istar+1,q)*pi
#    temp = temp - (abs(q-qstar))
    return temp


def relative_entropy(p,q):
    return p*numpy.log(p/q)+(1-p)*numpy.log((1-p)/(1-q))

u_attacker=[0]*(100)
qhat=[0]*(50)
exponent =[0]*(50)
k=0
for n in my_range(10,500,10):
    #print n
# find istar and pi for the decision rule that has FA rate exactly equal to epsilon
    temp = 1
    for i in range(0,n+1):
        temp = temp - ncr(n,i)*prob(n,i,p)
        if temp <= epsilon:
            break
    istar = i-1 #randomisation is at threshold istar+1
    temp = 1
    for i in range(0,istar+2):
        temp = temp - ncr(n,i)*prob(n,i,p)
    pi = (epsilon - temp) /(ncr(n,istar+1)*prob(n,istar+1,p))
    
    for j in range(0,100):
        q = q1+j*(q2-q1)/100
        u_attacker[j] = type2error(n,q,istar,pi) - (abs(q-qstar))
    jstar = numpy.argmax(u_attacker)
    qhat[k] = q1+jstar*(q2-q1)/100
    exponent[k] = -1*numpy.log(type2error(n,qhat[k],istar,pi))/n
    k = k+1    

matplotlib.rcParams.update({'font.size': 12})
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.20)
plt.plot(numpy.linspace(10,500,50),exponent)
plt.xlabel('n',fontsize=14)
plt.ylabel(r'$-\frac{1}{n} \log\, e_n(\hat{q}_n, \hat{\varphi}_n) $', fontsize=14)
plt.ylim([0,0.25])
#plt.savefig('test1.png')
plt.close()
plt.plot(numpy.linspace(10,500,50),qhat)
plt.gcf().subplots_adjust(bottom=0.15)
plt.gcf().subplots_adjust(left=0.2)
plt.xlabel('n',fontsize=14)
plt.ylabel(r'$\hat{q}_n$', fontsize=14)
#plt.savefig('test2.png')
