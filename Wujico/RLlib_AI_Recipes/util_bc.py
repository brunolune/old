#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 16:07:21 2019

@author: Bruno & Eric
"""
#from itertools import intersection

flat_list = lambda l: [item for sublist in l for item in sublist]

import math

def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

def intersec(a, b): 
    a_set = set(a) 
    b_set = set(b)      
    return len(a_set.intersection(b_set))

def flat_tuple(T):
    if not isinstance(T, tuple): return(T,)
    elif len(T) == 0: return ()
    else: return flat_tuple(T[0]) + flat_tuple(T[1:])