#!python3
#
# Copyright (C) 2014-2015 Julius Susanto. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""
PYPOWER-Dynamics
4th Order Synchronous Machine Model

"""

import numpy as np

class sym_order4:
    def __init__(self, filename, dynopt):
        self.signals = {}
        self.states = {}
        self.states0 = {}
        self.dsteps = {}
        self.params = {}
        self.opt = dynopt['iopt']
        self.omega_n = 2 * np.pi * dynopt['fn']
        
        # Check for speed-voltage term option 
        if 'speed_volt' in dynopt:
            self.speed_volt = dynopt['speed_volt']
        else:
            self.speed_volt = False
        
        self.parser(filename)
        
        # Convert impedances and H to 100MVA base
        if 'MVA_Rating' in self.params.keys():
            base_mva = self.params['MVA_Rating']
            self.params['H'] = self.params['H'] * base_mva / 100
            self.params['Ra'] = self.params['Ra'] * 100 / base_mva
            self.params['Xd'] = self.params['Xd'] * 100 / base_mva
            self.params['Xdp'] = self.params['Xdp'] * 100 / base_mva
            self.params['Xdpp'] = self.params['Xdpp'] * 100 / base_mva
            self.params['Xq'] = self.params['Xq'] * 100 / base_mva
            self.params['Xqp'] = self.params['Xqp'] * 100 / base_mva
            self.params['Xqpp'] = self.params['Xqpp'] * 100 / base_mva
        
        # Equivalent Norton impedance for Ybus modification
        self.Yg = (self.params['Ra'] - 1j * 0.5 * (self.params['Xdp'] + self.params['Xqp'])) / (self.params['Ra'] **2 + (self.params['Xdp'] * self.params['Xqp']))
    
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
        
        Vd0 = np.abs(vt0) * np.sin(delta0 - np.angle(vt0))
        Vq0 = np.abs(vt0) * np.cos(delta0 - np.angle(vt0))
        
        # Calculate machine state variables and Vfd
        Eqp0 = Vq0 + self.params['Ra'] * Iq0 + self.params['Xdp'] * Id0
        Edp0 = Vd0 + self.params['Ra'] * Id0 - self.params['Xqp'] * Iq0
        Vfd0 = np.abs(Eqp0) + (self.params['Xd'] - self.params['Xdp']) * Id0
        
        # Calculate active and reactive power
        p0 = (Vd0 + self.params['Ra']*Id0) * Id0 + (Vq0  + self.params['Ra']*Iq0) * Iq0
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
        
    def calc_currents(self,vt):
        """
        Calculate machine current injections (in network reference frame)
        """
        # Calculate terminal voltage in dq reference frame
        Vd = np.abs(vt) * np.sin(self.states['delta'] - np.angle(vt))
        Vq = np.abs(vt) * np.cos(self.states['delta'] - np.angle(vt))
        
        # Calculate Id and Iq (Norton equivalent current injection in dq frame)
        Eqp = self.states['Eqp']
        Edp = self.states['Edp']
        Ra = self.params['Ra']
        Xdp = self.params['Xdp']
        Xqp = self.params['Xqp']
        
        # Check if speed-voltage term should be included
        if self.speed_volt:
            omega = self.states['omega']
        else:
            omega = 1
        
        Id = (Eqp - Ra / (Xqp * omega) * (Vd - Edp) - Vq / omega) / (Xdp + Ra ** 2 / (omega * omega * Xqp))
        Iq = (Vd / omega + Ra * Id / omega - Edp) / Xqp
        
        # Calculate power output    
        p = (Vd + self.params['Ra']*Id) * Id + (Vq  + self.params['Ra']*Iq) * Iq             
        q = Vq * Id - Vd * Iq
        
        # Calculate machine current injection (Norton equivalent current injection in network frame)
        delta = self.states['delta']
        In = (Iq - 1j * Id) * np.exp(1j * (self.states['delta']))
        Im = In + self.Yg * vt
        
        # Update signals
        self.signals['Id'] = Id
        self.signals['Iq'] = Iq
        self.signals['Vd'] = Vd
        self.signals['Vq'] = Vq
        self.signals['P'] = p
        self.signals['Q'] = q
        self.signals['Vt'] = np.abs(vt)
        self.signals['Vang'] = np.angle(vt)
        
        return Im
        
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
        
        if round(dEdp,6) != 0 or round(dEqp,6) != 0:
            print('Warning: differential equations not zero on initialisation...')
            print('dEdp = ' + str(dEdp) + ', dEqp = ' + str(dEqp))
    
    def solve_step(self,h,dstep):
        """
        Solve machine differential equations for the next stage in the integration step
        """
        
        # Initial state variables
        omega_0 = self.states['omega']
        delta_0 = self.states['delta']
        Eqp_0 = self.states['Eqp']
        Edp_0 = self.states['Edp']
        
        Xd = self.params['Xd']
        Xdp = self.params['Xdp']
        Xq = self.params['Xq']
        Xqp = self.params['Xqp']
        Td0p = self.params['Td0p']
        Tq0p = self.params['Tq0p']
        
        Vfd = self.signals['Vfd']
        Id = self.signals['Id']
        Iq = self.signals['Iq']
        
        # Electrical differential equations
        f1 = (Vfd - (Xd - Xdp) * Id - Eqp_0) / Td0p
        k_Eqp = h * f1
        
        f2 = ((Xq - Xqp) * Iq - Edp_0) / Tq0p
        k_Edp = h * f2
        
        # Swing equation
        f3 = 1/(2 * self.params['H']) * (self.signals['Pm'] / omega_0 - self.signals['P'])
        k_omega = h * f3
        
        f4 = self.omega_n * (omega_0 - 1)
        k_delta = h * f4
        
        if self.opt == 'mod_euler':
            # Modified Euler
            # Update state variables
            if dstep == 0:
                # Predictor step
                self.states['Eqp'] = Eqp_0 + k_Eqp
                self.dsteps['Eqp'] = [k_Eqp]
                self.states['Edp'] = Edp_0 + k_Edp
                self.dsteps['Edp'] = [k_Edp]
                self.states['omega'] = omega_0 + k_omega
                self.dsteps['omega'] = [k_omega]
                self.states['delta'] = delta_0 + k_delta
                self.dsteps['delta'] = [k_delta]
            elif dstep == 1:
                # Corrector step
                self.states['Eqp'] = Eqp_0 + 0.5 * (k_Eqp - self.dsteps['Eqp'][0])
                self.states['Edp'] = Edp_0 + 0.5 * (k_Edp - self.dsteps['Edp'][0])
                self.states['omega'] = omega_0 + 0.5 * (k_omega - self.dsteps['omega'][0])     
                self.states['delta'] = delta_0 + 0.5 * (k_delta - self.dsteps['delta'][0])
                self.signals['Tm'] = self.signals['Pm'] / omega_0
        
        elif self.opt == 'runge_kutta':
            # 4th Order Runge-Kutta Method
            # Update state variables
            if dstep == 0:
                # Save initial states
                self.states0['omega'] = omega_0
                self.states0['delta'] = delta_0
                self.states0['Eqp'] = Eqp_0 
                self.states0['Edp'] = Edp_0
                
                self.states['Eqp'] = Eqp_0 + 0.5 * k_Eqp
                self.dsteps['Eqp'] = [k_Eqp]
                self.states['Edp'] = Edp_0 + 0.5 * k_Edp
                self.dsteps['Edp'] = [k_Edp]
                self.states['omega'] = omega_0 + 0.5 * k_omega
                self.dsteps['omega'] = [k_omega]            
                self.states['delta'] = delta_0 + 0.5 * k_delta
                self.dsteps['delta'] = [k_delta]
            elif dstep == 1:
                self.states['Eqp'] = Eqp_0 + 0.5 * k_Eqp
                self.dsteps['Eqp'].append(k_Eqp)
                self.states['Edp'] = Edp_0 + 0.5 * k_Edp
                self.dsteps['Edp'].append(k_Edp)
                self.states['omega'] = omega_0 + 0.5 * k_omega
                self.dsteps['omega'].append(k_omega)           
                self.states['delta'] = delta_0 + 0.5 * k_delta
                self.dsteps['delta'].append(k_delta)
            elif dstep == 2:
                self.states['Eqp'] = Eqp_0 + k_Eqp
                self.dsteps['Eqp'].append(k_Eqp)
                self.states['Edp'] = Edp_0 + k_Edp
                self.dsteps['Edp'].append(k_Edp)
                self.states['omega'] = omega_0 + k_omega
                self.dsteps['omega'].append(k_omega)           
                self.states['delta'] = delta_0 + k_delta
                self.dsteps['delta'].append(k_delta)
            elif dstep == 3:
                self.states['Eqp'] = self.states0['Eqp'] + 1/6 * (self.dsteps['Eqp'][0] + 2*self.dsteps['Eqp'][1] + 2*self.dsteps['Eqp'][2] + k_Eqp)
                self.states['Edp'] = self.states0['Edp'] + 1/6 * (self.dsteps['Edp'][0] + 2*self.dsteps['Edp'][1] + 2*self.dsteps['Edp'][2] + k_Edp)
                self.states['omega'] = self.states0['omega'] + 1/6 * (self.dsteps['omega'][0] + 2*self.dsteps['omega'][1] + 2*self.dsteps['omega'][2] + k_omega)
                self.states['delta'] = self.states0['delta'] + 1/6 * (self.dsteps['delta'][0] + 2*self.dsteps['delta'][1] + 2*self.dsteps['delta'][2] + k_delta)
                self.signals['Tm'] = self.signals['Pm'] / omega_0
   