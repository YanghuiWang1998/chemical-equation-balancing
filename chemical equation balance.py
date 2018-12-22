#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 16:50:00 2018

@author: yanghuiwang
"""


#1. turn a string into a chemical equation (regular expression)

#2. turn a chemical equation into a system of equations

#3. solve linear system (matrix algebra)

import sympy as sp
sp.init_printing(use_unicode=True)
import re
from collections import OrderedDict


def balance(eq):
    s=str(eq) # convert the equation into a string
    n=len(re.findall(r'\W',s))+2 # width of the matrix
    x=sp.symbols('x0:%d'%(n-1))
    left=re.findall(r'(.+)=',s)[0] # the left side of the equation
    leftelem=left.split('+') # the elements on the left side of the equation
    right=re.findall(r'=(.+)',s)[0] # the right side of the equation
    rightelem=right.split('+') # the elements on the right side of the equation
    # find all unique elements and preserve the order 
    elem=re.findall(r'([A-Z](?:[a-z]+)?)[0-9]?',left)
    seen = set()
    uniqueelem=[i for i in elem if not (i in seen or seen.add(i))]
    # construct an ordered dictionary with the elements and the coefficients
    emptylist=[[0]*n for i in range(len(uniqueelem))]
    dictionary=OrderedDict(zip(uniqueelem,emptylist))

    # fill the matrix of the left side
    index=0
    for i in leftelem:
        k=re.findall(r'[A-Z][a-z]?[0-9]?',i) # get the list of elements in the component 
        for j in range(len(k)):  # for each element in that component, 
            if re.findall(r'[0-9]',k[j]): # if a number follows the element, extract that number to be the coefficent
                word=re.findall(r'[A-z]',k[j])[0]
                dictionary[word][index]=int(re.findall(r'[0-9]',k[j])[0])
            else: # if there is no number following the component
                dictionary[k[j]][index]+=1 # fill number 1 into the matrix
        index+=1
    
    # fill the matrix of the right side (with negative signs)
    indexright=len(leftelem) # start counting after the left side ends
    for i in rightelem:
        k=re.findall(r'[A-Z][a-z]?[0-9]?',i) # get the list of elements in the component 
        for j in range(len(k)): # for each element in that component, 
            if re.findall(r'[0-9]',k[j]): # if a number follows the element
                word=re.findall(r'[A-z]',k[j])[0]
                dictionary[word][indexright]=-int(re.findall(r'[0-9]',k[j])[0])
            else: # if there is no number following the component
                dictionary[k[j]][indexright]+=-1 # fill number -1 into the matrix
        indexright+=1

    M=sp.Matrix(dictionary.values()) # get the matrix

    sols=sp.solve_linear_system(M,*x) # solve the matrix

    for j in range(n-2):
        sols[x[j]]=sols[x[j]].subs({'x%d'%(n-2):1}) # substitute the last element of sols dictionary with 1 with a for loop
    denom=sp.lcm([sp.fraction(sols[i])[1] for i in sols]) # get the largest common factor
    for i in sols:
        sols[i]*=denom # update the value of dictionary with a for loop
    sols.update({sp.symbols('x%d'%(n-2)):denom}) # add the last element into the sols dictionary

    #concagate left elements with coefficients
    leftelem=[str(sols[sp.symbols('x%d'%i)])+leftelem[i] if sols[sp.symbols('x%d'%i)]!=1 else leftelem[i] for i in range(len(leftelem)) ]
    #concagate right elements with coefficients
    rightelem=[str(sols[sp.symbols('x%d'%(i+len(leftelem)))])+rightelem[i] if sols[sp.symbols('x%d'%(i+len(leftelem)))]!=1 else rightelem[i] for i in range(len(rightelem))  ]
    # concagate left and right sides
    return '+'.join(leftelem)+'='+'+'.join(rightelem)

print balance("H2+O2=H2O")
print balance("PhCH3+KMnO4+H2SO4=PhCOOH+K2SO4+MnSO4+H2O")