#!python3
#
# Copyright (C) 2014 Julius Susanto
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
Functions for standard blocks
(using modified Euler method for integration)

"""

import numpy as np
from integrators import mod_euler, runge_kutta

# Lag block
# K / (1 + sT)
# p = [K, T]
def lag_block(h,x0,yi,p):   
    f = '(yi - x) / p[1]'
    x1 = mod_euler(x0,h,f,yi,p)
    yo = p[0] * x1
            
    return yo, x1

# Integrator block    
# K / sT    
# p = [K, T]
def int_block(h,x0,yi,p):
    f = 'yi * p[0] / p[1]'
    x1 = mod_euler(x0,h,f,yi,p)
    yo = x1
            
    return yo, x1

# Lead-Lag block
# (1 + sTa) / (1 + sTb)
# p = [Ta, Tb]
def leadlag_block(h,x0,yi,p):  
    f = '(yi - x) / p[1]'
    x1 = mod_euler(x0,h,f,yi,p)
    yo = x1 + p[0] * (yi - x0) / p[1]
    
    return yo, x1
  
# Summation block    
# yo = yi1 + yi2 + ... + yin
# yi = [yi1, yi2, ... yin]
def sum_block(yi):
    yo = sum(yi)
    
    return yo
       
# Gain block    
# yo = p * yi
# p is a scalar gain coefficient
def gain_block(yi,p):
    yo = p * yi
    
    return yo
  
# Multiplication block    
# yo = yi1 * yi2 * ... * yin
# yi = [yi1, yi2, ... yin]
def mult_block(yi):
    yo = np.prod(yi)
    
    return yo

    
