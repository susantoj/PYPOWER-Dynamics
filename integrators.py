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
Generic explicit numerical integrators

"""

import numpy as np

def integrate(x0,h,f,yi,p,opt):
    """
    Integrator function wrapper
    """

    if opt == 'mod_euler':
        x1 = mod_euler(x0,h,f,yi,p)
    elif opt == 'runge_kutta':    
        x1 = runge_kutta(x0,h,f,yi,p)
    else:
        print('WARNING: integration algorithm ' + str(opt) + ' not found...')
    
    return x1
    
def mod_euler(x0,h,f,yi,p):
    """
    Modified Euler integration step
    
        x1 = mod_euler(x0,h,f,yi,p)
    
    x0 is the initial value of the state variable at time t0
    h is the time step
    f is the ordinary differential equation to be evaluated and is expressed in
      the form of a string, e.g. 'x * np.cos(p[0] * t) + y[0]'
    yi is a vector of inputs that can be used in f
    p is a vector of coefficients or parameters that can be used in f
    
    Returns x1, the computed state variable at time t0 + h
    """
    
    x = x0
    K1 = h * eval(f)
    
    x = x0 + K1
    K2 = h * eval(f)
    
    x1 = x0 + (K1 + K2) / 2
    
    return x1

def runge_kutta(x0,h,f,yi,p):
    """
    4th order Runge-Kutta integration step (based on Simpson's rule)
    
        x1 = runge_kutta(x0,h,f,yi,p)
    
    x0 is the initial value of the state variable at time t0
    h is the time step
    f is the ordinary differential equation to be evaluated and is expressed in
      the form of a string, e.g. 'x * np.cos(p[0] * t) + y[0]'
    yi is a vector of inputs that can be used in f
    p is a vector of coefficients or parameters that can be used in f
    
    Returns x1, the computed state variable at time t0 + h
    """
    
    x = x0
    K1 = h * eval(f)
    
    y = x0 + K1/2
    K2 = h * eval(f)
    
    x = x0 + K2/2
    K3 = h * eval(f)
    
    x = x0 + K3
    K4 = h * eval(f)
    
    x1 = x0 + (K1 + 2 * K2 + 2 * K3 + K4) / 6
    
    return x1