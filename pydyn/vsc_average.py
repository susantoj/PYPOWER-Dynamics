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
Voltage Source Converter Model Class
Average model of a VSC in voltage-control mode (i.e. controlled voltage source behind an impedance). 
"""

import numpy as np

class vsc_average:
    def __init__(self, ID, gen_no, Rl, Xl, dynopt):
        self.id = ID
        self.gen_no = gen_no
        self.opt = dynopt['iopt']
        
        self.signals = {}
        self.states = {}
        self.states0 = {}
        self.dsteps = {}
        
        self.params = {}      
        self.params['Rl'] = Rl
        self.params['Xl'] = Xl
        self.params['fn'] = dynopt['fn']
        
        # Equivalent Norton impedance for Ybus modification
        self.Yg = 1 / (Rl + 1j * Xl)
        
    def initialise(self,vt0,S0):
        """
        Initialise converter emf based on load flow voltage and grid current injection
        """
        # Calculate initial armature current
        Ia0 = np.conj(S0 / vt0)
        phi0 = np.angle(Ia0)
        
        # Calculate steady state machine emf (i.e. voltage behind synchronous reactance)
        Edq0 = vt0 + (self.params['Rl'] + 1j * self.params['Xl']) * Ia0
        delta0 = np.angle(Edq0)
        
        # Convert currents to rotor reference frame
        Id0 = np.abs(Ia0) * np.sin(delta0 - phi0)
        Iq0 = np.abs(Ia0) * np.cos(delta0 - phi0)
               
        # Initialise signals, states and parameters
        self.signals['Vt'] = np.abs(vt0)
        
        self.signals['Edq'] = Edq0
        self.signals['Ed'] = np.real(Edq0)
        self.signals['Eq'] = np.imag(Edq0)
        self.signals['Id'] = Id0
        self.signals['Iq'] = Iq0    

        self.states['delta'] = delta0
        self.states['omega'] = 1

    def calc_currents(self, vt):
        """
        Solve grid current injections (in network reference frame)
        """
        Edq = self.signals['Ed'] + 1j * self.signals['Eq']
        delta = np.angle(Edq)
                
        # Calculate terminal voltage in dq reference frame
        Vd = np.abs(vt) * np.sin(self.states['delta'] - np.angle(vt))
        Vq = np.abs(vt) * np.cos(self.states['delta'] - np.angle(vt))
        
        # Calculate Id and Iq (Norton equivalent current injection in dq frame)
        Ia = (Edq - vt) / (self.params['Rl'] + 1j * self.params['Xl'])
        phi = np.angle(Ia)
        Id = np.abs(Ia) * np.sin(delta - phi)
        Iq = np.abs(Ia) * np.cos(delta - phi)
        
        # Calculate machine current injection (Norton equivalent current injection in network frame)
        In = (Iq - 1j * Id) * np.exp(1j * (self.states['delta']))
        Im = In + self.Yg * vt
        
        # Update signals
        self.signals['Edq'] = Edq
        self.signals['Vd'] = Vd
        self.signals['Vq'] = Vq
        self.signals['Id'] = Id
        self.signals['Iq'] = Iq        
        self.signals['Vt'] = np.abs(vt) 

        self.states['delta'] = delta
        
        return Im
        
    def solve_step(self,h,dstep):
        """
        Solve machine differential equations for the next stage in the integration step
        """
        
        # State variables do not change in this model
        pass
        
    
    