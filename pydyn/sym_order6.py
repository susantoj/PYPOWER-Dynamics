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
6th Order Synchronous Machine Model
Based on Anderson-Fouad model
Anderson, P. M., Fouad, A. A., "Power System Control and Stability", Wiley-IEEE Press, New York, 2002
"""

import numpy as np

class sym_order6:
    def __init__(self, filename, iopt):
        self.id = ''
        self.gen_no = 0
        self.signals = {}
        self.states = {}
        self.states0 = {}
        self.dsteps = {}
        self.params = {}
        self.opt = iopt
        self.omega_n = 2 * np.pi * 50
        
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
        Ia0 = np.conj(S0 / vt0)
        phi0 = np.angle(Ia0)
        
        # Calculate steady state machine emf (i.e. voltage behind synchronous reactance)
        Eq0 = vt0 + np.complex(self.params['Ra'],self.params['Xq']) * Ia0
        delta0 = np.angle(Eq0)
        
        # Convert currents to rotor reference frame
        Id0 = np.abs(Ia0) * np.sin(delta0 - phi0)
        Iq0 = np.abs(Ia0) * np.cos(delta0 - phi0)
        
        # Calculate machine state variables and Vfd
        Vfd0 = np.abs(Eq0) + (self.params['Xd'] - self.params['Xq']) * Id0
        Eqp0 = Vfd0 - (self.params['Xd'] - self.params['Xdp']) * Id0
        Eqpp0 = Eqp0 - (self.params['Xdp'] - self.params['Xdpp']) * Id0
        
        Edp0 = (self.params['Xq'] - self.params['Xqp']) * Iq0 
        Edpp0 = Edp0 + (self.params['Xqp'] - self.params['Xqpp']) * Iq0     
        
        Vd0 = Edpp0 + self.params['Xqpp'] * Iq0 - self.params['Ra'] * Id0
        Vq0 = Eqpp0 - self.params['Xdpp'] * Id0 - self.params['Ra'] * Iq0
        
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
        self.states['Eqpp'] = Eqpp0
        self.states['Edp'] = Edp0
        self.states['Edpp'] = Edpp0
    
    def solve_step(self,h,dstep):
        """
        Solve machine differential equations for the next step in modified Euler method iteration
        """
        
        # Initial state variables
        omega_0 = self.states['omega']
        delta_0 = self.states['delta']
        Eqp_0 = self.states['Eqp']
        Edp_0 = self.states['Edp']
        Edpp_0 = self.states['Edpp']
        Eqpp_0 = self.states['Eqpp']
        
        # Electrical differential equations
        p = [self.params['Xd'], self.params['Xdp'], self.params['Td0p']]
        yi = [self.signals['Vfd'], self.signals['Id']]
        f1 = (yi[0] - (p[0] - p[1]) * yi[1] - Eqp_0) / p[2]
        k_Eqp = h * f1
        
        p = [self.params['Xq'], self.params['Xqp'], self.params['Tq0p']]
        yi = self.signals['Iq']
        f2 = ((p[0] - p[1]) * yi - Edp_0) / p[2]
        k_Edp = h * f2
        
        p = [self.params['Xdp'], self.params['Xdpp'], self.params['Td0pp']]
        yi = [Eqp_0, self.signals['Id']]
        f3 = (yi[0] - (p[0] - p[1]) * yi[1] - Eqpp_0) / p[2]
        k_Eqpp = h * f3
        
        p = [self.params['Xqp'], self.params['Xqpp'], self.params['Tq0pp']]
        yi = [Edp_0, self.signals['Iq']]
        f4 = (yi[0] + (p[0] - p[1]) * yi[1] - Edpp_0) / p[2]
        k_Edpp = h * f4
        
        # Swing equation
        f5 = 1/(2 * self.params['H']) * (self.signals['Pm'] - self.signals['P'])
        k_omega = h * f5
        
        f6 = self.omega_n * (omega_0 - 1)
        k_delta = h * f6
        
        if self.opt == 'mod_euler':
            # Modified Euler
            # Update state variables
            if dstep == 0:
                # Predictor step
                self.states['Eqp'] = Eqp_0 + k_Eqp
                self.dsteps['Eqp'] = [k_Eqp]
                self.states['Edp'] = Edp_0 + k_Edp
                self.dsteps['Edp'] = [k_Edp]
                self.states['Eqpp'] = Eqpp_0 + k_Eqpp
                self.dsteps['Eqpp'] = [k_Eqpp]
                self.states['Edpp'] = Edpp_0 + k_Edpp
                self.dsteps['Edpp'] = [k_Edpp]
                self.states['omega'] = omega_0 + k_omega
                self.dsteps['omega'] = [k_omega]
                self.states['delta'] = delta_0 + k_delta
                self.dsteps['delta'] = [k_delta]
            else:
                # Corrector step
                self.states['Eqp'] = Eqp_0 + 0.5 * (k_Eqp - self.dsteps['Eqp'][0])
                self.states['Edp'] = Edp_0 + 0.5 * (k_Edp - self.dsteps['Edp'][0])
                self.states['Eqpp'] = Eqpp_0 + 0.5 * (k_Eqpp - self.dsteps['Eqpp'][0])
                self.states['Edpp'] = Edpp_0 + 0.5 * (k_Edpp - self.dsteps['Edpp'][0])
                self.states['omega'] = omega_0 + 0.5 * (k_omega - self.dsteps['omega'][0])     
                self.states['delta'] = delta_0 + 0.5 * (k_delta - self.dsteps['delta'][0])
        
        elif self.opt == 'runge_kutta':
            # 4th Order Runge-Kutta Method
            # Update state variables
            if dstep == 0:
                # Save initial states
                self.states0['omega'] = omega_0
                self.states0['delta'] = delta_0
                self.states0['Eqp'] = Eqp_0 
                self.states0['Edp'] = Edp_0
                self.states0['Eqpp'] = Eqpp_0 
                self.states0['Edpp'] = Edpp_0
                
                self.states['Eqp'] = Eqp_0 + 0.5 * k_Eqp
                self.dsteps['Eqp'] = [k_Eqp]
                self.states['Edp'] = Edp_0 + 0.5 * k_Edp
                self.dsteps['Edp'] = [k_Edp]
                self.states['Eqpp'] = Eqpp_0 + 0.5 * k_Eqpp
                self.dsteps['Eqpp'] = [k_Eqpp]
                self.states['Edpp'] = Edpp_0 + 0.5 * k_Edpp
                self.dsteps['Edpp'] = [k_Edpp]
                self.states['omega'] = omega_0 + 0.5 * k_omega
                self.dsteps['omega'] = [k_omega]            
                self.states['delta'] = delta_0 + 0.5 * k_delta
                self.dsteps['delta'] = [k_delta]
            elif dstep == 1:
                self.states['Eqp'] = Eqp_0 + 0.5 * k_Eqp
                self.dsteps['Eqp'].append(k_Eqp)
                self.states['Edp'] = Edp_0 + 0.5 * k_Edp
                self.dsteps['Edp'].append(k_Edp)
                self.states['Eqpp'] = Eqpp_0 + 0.5 * k_Eqpp
                self.dsteps['Eqpp'].append(k_Eqpp)
                self.states['Edpp'] = Edpp_0 + 0.5 * k_Edpp
                self.dsteps['Edpp'].append(k_Edpp)
                self.states['omega'] = omega_0 + 0.5 * k_omega
                self.dsteps['omega'].append(k_omega)           
                self.states['delta'] = delta_0 + 0.5 * k_delta
                self.dsteps['delta'].append(k_delta)
            elif dstep == 2:
                self.states['Eqp'] = Eqp_0 + k_Eqp
                self.dsteps['Eqp'].append(k_Eqp)
                self.states['Edp'] = Edp_0 + k_Edp
                self.dsteps['Edp'].append(k_Edp)
                self.states['Eqpp'] = Eqpp_0 + k_Eqpp
                self.dsteps['Eqpp'].append(k_Eqpp)
                self.states['Edpp'] = Edpp_0 + k_Edpp
                self.dsteps['Edpp'].append(k_Edpp)
                self.states['omega'] = omega_0 + k_omega
                self.dsteps['omega'].append(k_omega)           
                self.states['delta'] = delta_0 + k_delta
                self.dsteps['delta'].append(k_delta)
            elif dstep == 3:
                self.states['Eqp'] = self.states0['Eqp'] + 1/6 * (self.dsteps['Eqp'][0] + 2*self.dsteps['Eqp'][1] + 2*self.dsteps['Eqp'][2] + k_Eqp)
                self.states['Edp'] = self.states0['Edp'] + 1/6 * (self.dsteps['Edp'][0] + 2*self.dsteps['Edp'][1] + 2*self.dsteps['Edp'][2] + k_Edp)
                self.states['Eqpp'] = self.states0['Eqpp'] + 1/6 * (self.dsteps['Eqpp'][0] + 2*self.dsteps['Eqpp'][1] + 2*self.dsteps['Eqpp'][2] + k_Eqpp)
                self.states['Edpp'] = self.states0['Edpp'] + 1/6 * (self.dsteps['Edpp'][0] + 2*self.dsteps['Edpp'][1] + 2*self.dsteps['Edpp'][2] + k_Edpp)
                self.states['omega'] = self.states0['omega'] + 1/6 * (self.dsteps['omega'][0] + 2*self.dsteps['omega'][1] + 2*self.dsteps['omega'][2] + k_omega)
                self.states['delta'] = self.states0['delta'] + 1/6 * (self.dsteps['delta'][0] + 2*self.dsteps['delta'][1] + 2*self.dsteps['delta'][2] + k_delta)
            
    def calc_currents(self,vt):
        """
        Calculate machine current injections (in network reference frame)
        """
        
        # Calculate terminal voltage in dq reference frame
        Vd = np.abs(vt) * np.sin(self.states['delta'] - np.angle(vt))
        Vq = np.abs(vt) * np.cos(self.states['delta'] - np.angle(vt))
        
        # For machines with no saliency
        # Calculate Id and Iq (Norton equivalent current injection in dq frame)
        if self.params['Ra'] > 0:
            Iq = (-self.params['Ra'] * (Vq-self.states['Eqpp']) + self.params['Xdpp'] * (Vd - self.states['Edpp'])) / \
                    (self.params['Xdpp'] * self.params['Xqpp'] + self.params['Ra'] ** 2)
            Id = -(Vd - self.states['Edpp'] - self.params['Xqpp'] * Iq) / self.params['Ra']
        else:
            # Ra = 0 (or if Ra is negative, Ra is ignored)
            Id = (self.states['Eqpp'] - Vq) / self.params['Xdpp']
            Iq = (Vd - self.states['Edpp']) / self.params['Xqpp']
        
        # Formulation for salient pole machines
        # (only for machines with Ra = 0 at the moment)
        if self.params['Ra'] == 0:
            # Unadjusted current injection (in dq frame)
            Y0 = 0.5 * (self.params['Xdpp'] + self.params['Xqpp']) / (self.params['Xdpp'] * self.params['Xqpp'])
            Id_u = (self.states['Eqpp'] - Vq) * Y0
            Iq_u = -(self.states['Edpp'] - Vd) * Y0
            
            # Adjusted current injection (in dq frame)
            Y1 = 0.5 * (self.params['Xdpp'] - self.params['Xqpp']) / (self.params['Xdpp'] * self.params['Xqpp'])
            Id_a = Y1 * (self.states['Eqpp'] - Vq)
            Iq_a = Y1 * (self.states['Edpp'] - Vd)
            
            # Total current injection
            Id = Id_u + Id_a
            Iq = Iq_u + Iq_a
        
        # Calculate power output
        p = Vd * Id + Vq * Iq
        q = Vq * Id - Vd * Iq
        S = np.complex(p,q)
        
        # Calculate machine current injection
        Im = (self.states['Eqpp'] - 1j * self.states['Edpp']) * np.exp(1j * (self.states['delta'])) / (1j * self.params['Xdpp'])
        
        # Update signals
        self.signals['Id'] = Id
        self.signals['Iq'] = Iq
        self.signals['Vd'] = Vd
        self.signals['Vq'] = Vq
        self.signals['P'] = p
        self.signals['Q'] = q
        self.signals['Vt'] = np.sqrt(self.signals['Vq'] **2 + self.signals['Vd'] **2)
        
        return Im