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
4th Order Synchronous Machine Model

"""

import numpy as np

class sym_order4:
    def __init__(self, filename, iopt):
        self.signals = {}
        self.states = {}
        self.dsteps = {}
        self.params = {}
        self.opt = iopt
        
        self.parser(filename)
    
    def parser(self, filename):
        """
        Parse a machine file (*.mach) and populate dictionary of parameters
        """
        f = open(filename, 'r')
        
        for line in f:
            if line[0] != '#' and line.strip() != '':   # Ignore comments and blank lines
                tokens = line.strip().split('=')
                if tokens[0].strip() == 'ID':
                    self.id = tokens[1].strip()
                elif tokens[0].strip() == 'GEN_NO':
                    self.gen_no = int(tokens[1].strip())
                else:
                    self.params[tokens[0].strip()] = float(tokens[1].strip())
                
        f.close()
    
    def initialise(self, vt0, S0):
        """
        Initialise machine signals and states based on load flow voltage and complex power injection
        """

        # Calculate initial armature current
        Ia0 =  np.conj(S0 / vt0)
        phi0 = np.angle(Ia0)
        
        # Calculate steady state machine emf (i.e. voltage behind synchronous reactance)
        Eq0 = vt0 + np.complex(self.params['Ra'],self.params['Xq']) * Ia0
        delta0 = np.angle(Eq0)
        
        # Convert currents to rotor reference frame
        Id0 = np.abs(Ia0) * np.sin(delta0 - phi0)
        Iq0 = np.abs(Ia0) * np.cos(delta0 - phi0)
        
        # Calculate machine state variables and Vfd
        Vfd0 = np.abs(Eq0) + (self.params['Xd'] - self.params['Xq']) * Id0
        
        # Initial transient EMF
        Eqp0 = Vfd0 - (self.params['Xd'] - self.params['Xdp']) * Id0
        Edp0 = (self.params['Xq'] - self.params['Xqp']) * Iq0   
        
        # Initial Vd, Vq
        Vd0 = Edp0 + self.params['Xqp'] * Iq0
        Vq0 = Eqp0 - self.params['Xdp'] * Id0
        
        # Calculate active and reactive power
        p0 = Vd0 * Id0 + Vq0 * Iq0
        q0 = Vq0 * Id0 - Vd0 * Iq0
        
        # Initialise signals, states and parameters        
        self.signals['Vfd'] = Vfd0
        self.signals['Id'] = Id0
        self.signals['Iq'] = Iq0
        self.signals['Vd'] = Vd0
        self.signals['Vq'] = Vq0
        self.signals['Vt'] = np.abs(vt0)
        self.signals['P'] = p0
        self.signals['Q'] = q0
        self.signals['Pm'] = p0
        
        self.states['omega'] = 1
        self.states['delta'] = delta0
        self.states['Eqp'] = Eqp0
        self.states['Edp'] = Edp0
        
        self.check_diffs()
    
    def check_diffs(self):
        """
        Check if differential equations are zero (on initialisation)
        """
    
        # State variables
        Eqp_0 = self.states['Eqp']
        Edp_0 = self.states['Edp']
        
        Vfd = self.signals['Vfd']
        Id = self.signals['Id']
        Iq = self.signals['Iq']
        
        Xd = self.params['Xd']
        Xdp = self.params['Xdp']
        Td0p = self.params['Td0p']
        
        Xq = self.params['Xq']
        Xqp = self.params['Xqp']
        Tq0p = self.params['Tq0p']
        
        dEqp = (Vfd - (Xd - Xdp) * Id - Eqp_0) / Td0p
        dEdp = ((Xq - Xqp) * Iq - Edp_0) / Tq0p
        
        if dEdp != 0 or dEqp != 0:
            print('Differential equations not zero on initialisation...')
    
    def solve_step(self,h,dstep):
        """
        Solve machine differential equations for the next time step
        """
        
        # Initial state variables
        omega_0 = self.states['omega']
        delta_0 = self.states['delta']
        Eqp_0 = self.states['Eqp']
        Edp_0 = self.states['Edp']
        
        # Electrical differential equations
        p = [self.params['Xd'], self.params['Xdp'], self.params['Td0p']]
        yi = [self.signals['Vfd'], self.signals['Id']]
        f1 = (yi[0] - (p[0] - p[1]) * yi[1] - Eqp_0) / p[2]
        Eqp_1 = Eqp_0 + h * f1
        
        p = [self.params['Xq'], self.params['Xqp'], self.params['Tq0p']]
        yi = self.signals['Iq']
        f2 = ((p[0] - p[1]) * yi - Edp_0) / p[2]
        Edp_1 = Edp_0 + h * f2
        
        # Swing equation
        f3 = 1/(2 * self.params['H']) * (self.signals['Pm'] - self.signals['P'])
        omega_1 = omega_0 + h * f3
        
        f4 = 314.16 * (omega_0 - 1)
        delta_1 = delta_0 + h * f4
        
        # Update state variables
        if dstep == 0:
            # Predictor step
            self.states['Eqp'] = Eqp_1
            self.dsteps['Eqp'] = f1
            self.states['Edp'] = Edp_1
            self.dsteps['Edp'] = f2
            self.states['omega'] = omega_1
            self.dsteps['omega'] = f3
            self.states['delta'] = delta_1
            self.dsteps['delta'] = f4
        else:
            # Corrector step
            self.states['Eqp'] = Eqp_1 - h/2 * self.dsteps['Eqp']
            self.states['Edp'] = Edp_1 - h/2 * self.dsteps['Edp']
            self.states['omega'] = omega_1 - h/2 * self.dsteps['omega']
            self.states['delta'] = delta_1 - h/2 * self.dsteps['delta']
    
    def calc_currents(self,vt):
        """
        Calculate machine current injections (in network reference frame)
        """
        # Calculate terminal voltage in dq reference frame
        Vd = np.abs(vt) * np.sin(self.states['delta'] - np.angle(vt))
        Vq = np.abs(vt) * np.cos(self.states['delta'] - np.angle(vt))
        
        # Calculate Id and Iq
        Id = (self.states['Eqp'] - Vq) / self.params['Xdp']
        Iq = (Vd - self.states['Edp']) / self.params['Xqp']
        
        # Calculate power output
        p = Vd * Id + Vq * Iq        
        # Equivalent formulation
        #p = self.states['Eqp'] * Iq + self.states['Edp'] * Id + (self.params['Xdp'] - self.params['Xqp']) * Id * Iq
        
        q = Vq * Id - Vd * Iq
        S = p - 1j * q
        
        # Calculate machine current injection (Norton equivalent current injection in network frame)
        Im = (self.states['Eqp'] - 1j * self.states['Edp']) * np.exp(1j * (self.states['delta'])) / (1j * self.params['Xdp'])
        
        # Update signals
        self.signals['Id'] = Id
        self.signals['Iq'] = Iq
        self.signals['Vd'] = Vd
        self.signals['Vq'] = Vq
        self.signals['P'] = p
        self.signals['Q'] = q
        self.signals['Vt'] = np.sqrt(Vd**2 + Vq**2)
        
        return Im