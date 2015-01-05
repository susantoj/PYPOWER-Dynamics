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
External Grid Model Class
Grid is modelled as a constant voltage behind a transient reactance
and two differential equations representing the swing equations.
"""

import numpy as np
from pydyn.integrators import integrate

class ext_grid:
    def __init__(self, ID, gen_no, Xdp, H, iopt):
        self.id = ID
        self.gen_no = gen_no
        self.opt = iopt
        
        self.signals = {}
        self.states = {}
        
        self.params = {}          
        self.params['Xdp'] = Xdp
        self.params['H'] = H        
        
    def initialise(self,vt0,S0):
        """
        Initialise grid emf based on load flow voltage and grid current injection
        """
        # Calculate initial armature current
        Ia0 = np.conj(S0 / vt0)
        phi0 = np.angle(Ia0)
        
        # Calculate steady state machine emf (i.e. voltage behind synchronous reactance)
        Eq0 = vt0 + np.complex(0,self.params['Xdp']) * Ia0
        delta0 = np.angle(Eq0)
        
        p0 = 1 / self.params['Xdp'] * np.abs(vt0) * np.abs(Eq0) * np.sin(delta0 - np.angle(vt0))
        
        # Initialise signals, states and parameters
        self.signals['Vt'] = np.abs(vt0)
        self.signals['P'] = p0
        self.signals['Pm'] = p0
        
        self.states['Eq'] = np.abs(Eq0)       
        self.states['omega'] = 1
        self.states['delta'] = delta0
        
    def solve_step(self,h):
        """
        Solve differential equations for the next time step
        """
        # State variables
        omega_0 = self.states['omega']
        delta_0 = self.states['delta']
        
        # Solve swing equation
        p = self.params['H']
        yi = [self.signals['Pm'], self.signals['P']]
        f = '1/(2 * p) * (yi[0] - yi[1])'
        omega_1 = integrate(omega_0,h,f,yi,p,self.opt)
        
        p = self.params['H']
        yi = omega_0
        f = '314.16 * (yi - 1)'
        delta_1 = integrate(delta_0,h,f,yi,p,self.opt)
        
        # Update state variables
        self.states['omega'] = omega_1
        self.states['delta'] = delta_1
    
    def calc_currents(self, vt):
        """
        Solve grid current injections (in network reference frame)
        """
        delta = self.states['delta']
        Eq = self.states['Eq']
        
        p = 1 / self.params['Xdp'] * np.abs(vt) * np.abs(Eq) * np.sin(delta - np.angle(vt))
        
        # Update signals
        self.signals['P'] = p
        self.signals['Vt'] = np.abs(vt)
        
        i_grid = Eq * np.exp(1j * delta) / np.complex(0,self.params['Xdp'])
        
        return i_grid
    
    