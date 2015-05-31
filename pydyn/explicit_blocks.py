#!python3
#
# Copyright (C) 2014-2015 Julius Susanto. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""
PYPOWER-Dynamics
Functions for standard blocks (solves a step)

"""

import numpy as np

# Gain block    
# yo = p * yi
# p is a scalar gain coefficient
def gain_block(yi,p):
    yo = p * yi
    
    return yo
    
# Integrator block    
# K / sT    
# p = [K, T]
def int_block(h,x0,yi,p):
    f = yi * p[0] / p[1]
    x1 = x0 + h * f
    yo = x1
            
    return yo, x1, f

# Lag block
# K / (1 + sT)
# p = [K, T]
def lag_block(h,x0,yi,p):   
    f = (yi - x0) / p[1]
    x1 = x0 + h * f
    yo = p[0] * x1
            
    return yo, x1, f
    
# Lead-Lag block
# (1 + sTa) / (1 + sTb)
# p = [Ta, Tb]
def leadlag_block(h,x0,yi,p):  
    f = (yi - x0) / p[1]
    x1 = x0 + h * f
    yo = x1 + p[0] * (yi - x0) / p[1]
    
    return yo, x1, f

# Limiter block    
# yo = min_lim, if yi < min_lim
# yo = max_lim, if yi > max_lim
# yo = yi, min_lim <= yi <= max_lim
# p = [min_lim, max_lim]
def lim_block(yi,p):
    min_lim = p[0]
    max_lim = p[1]
    if yi < min_lim:
        yo = min_lim
    elif yi > max_lim:
        yo = max_lim
    else:
        yo = yi
    
    return yo
    
# Multiplication block    
# yo = yi1 * yi2 * ... * yin
# yi = [yi1, yi2, ... yin]
def mult_block(yi):
    yo = np.prod(yi)
    
    return yo
    
# Summation block    
# yo = yi1 + yi2 + ... + yin
# yi = [yi1, yi2, ... yin]
def sum_block(yi):
    yo = sum(yi)
    
    return yo
         
# Washout block
# (s / (1 + sT)
# p is the time constant T
def wout_block(h,x0,yi,p):  
    f = (yi - x0) / p
    x1 = x0 + h * f
    yo = (yi - x1) / p
    
    return yo, x1, f
    
