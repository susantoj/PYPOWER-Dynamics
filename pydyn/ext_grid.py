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

class ext_grid:
    def __init__(self, ID, gen_no, Xdp, H, dynopt):
        self.id = ID
        self.gen_no = gen_no
        self.opt = dynopt['iopt']
        
        self.signals = {}
        self.states = {}
        self.states0 = {}
        self.dsteps = {}
        
        self.params = {}          
        self.params['Xdp'] = Xdp
        self.params['H'] = H
        self.params['fn'] = dynopt['fn']
        
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

    def calc_currents(self, vt):
        """
        Solve grid current injections (in network reference frame)
        """
        delta = self.states['delta']
        Eq = self.states['Eq']
        Xdp = self.params['Xdp']
        
        p = np.abs(vt) * Eq * np.sin(delta - np.angle(vt)) / Xdp
        
        # Update signals
        self.signals['P'] = p
        self.signals['Vt'] = np.abs(vt)
        
        i_grid = Eq * np.exp(1j * delta) / np.complex(0,Xdp)
        
        return i_grid
        
    def solve_step(self,h,dstep):
        """
        Solve machine differential equations for the next stage in the integration step
        """
        # State variables
        omega_0 = self.states['omega']
        delta_0 = self.states['delta']
        
        # Solve swing equation
        f1 = 1/(2 * self.params['H']) * (self.signals['Pm'] / omega_0 - self.signals['P'])
        k_omega = h * f1
        
        f2 = 2 * np.pi * self.params['fn'] * (omega_0 - 1)
        k_delta = h * f2
        
        if self.opt == 'mod_euler':
            # Modified Euler
            # Update state variables
            if dstep == 0:
                self.states['omega'] = omega_0 + k_omega
                self.dsteps['omega'] = [k_omega]            
                self.states['delta'] = delta_0 + k_delta
                self.dsteps['delta'] = [k_delta]
            elif dstep == 1:
                self.states['omega'] = omega_0 + 0.5 * (k_omega - self.dsteps['omega'][0])     
                self.states['delta'] = delta_0 + 0.5 * (k_delta - self.dsteps['delta'][0]) 
        
        elif self.opt == 'runge_kutta':
            # 4th Order Runge-Kutta Method
            # Update state variables
            if dstep == 0:
                # Save initial states
                self.states0['omega'] = omega_0
                self.states0['delta'] = delta_0
                
                self.states['omega'] = omega_0 + 0.5 * k_omega
                self.dsteps['omega'] = [k_omega]            
                self.states['delta'] = delta_0 + 0.5 * k_delta
                self.dsteps['delta'] = [k_delta]
            elif dstep == 1:
                self.states['omega'] = omega_0 + 0.5 * k_omega
                self.dsteps['omega'].append(k_omega)           
                self.states['delta'] = delta_0 + 0.5 * k_delta
                self.dsteps['delta'].append(k_delta)
            elif dstep == 2:
                self.states['omega'] = omega_0 + k_omega
                self.dsteps['omega'].append(k_omega)           
                self.states['delta'] = delta_0 + k_delta
                self.dsteps['delta'].append(k_delta)
            elif dstep == 3:
                self.states['omega'] = self.states0['omega'] + 1/6 * (self.dsteps['omega'][0] + 2*self.dsteps['omega'][1] + 2*self.dsteps['omega'][2] + k_omega)
                self.states['delta'] = self.states0['delta'] + 1/6 * (self.dsteps['delta'][0] + 2*self.dsteps['delta'][1] + 2*self.dsteps['delta'][2] + k_delta)
        
    
    