#!python3
#
# Copyright (C) 2014-2015 Julius Susanto
#
# PYPOWER-Dynamics is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# PYPOWER-Dynamics is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PYPOWER-Dynamics. If not, see <http://www.gnu.org/licenses/>.

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
    
