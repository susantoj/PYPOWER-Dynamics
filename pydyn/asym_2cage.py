#!python3
#
# Copyright (C) 2014-2015 Julius Susanto. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""
PYPOWER-Dynamics
Double Cage Asynchronous Machine Model

Model equations based on section 15.2.5 of:
Milano, F., "Power System Modelling and Scripting", Springer-Verlag, 2010

"""

import numpy as np

class asym_2cage:
    def __init__(self, filename, dynopt):
        self.signals = {}
        self.states = {}
        self.states0 = {}
        self.dsteps = {}
        self.params = {}
        self.opt = dynopt['iopt']
        self.omega_n = 2 * np.pi * dynopt['fn']
        
        self.parser(filename)
        
        # Convert parameters to 100MVA base
        if 'MVA_Rating' in self.params.keys():
            self.base_mva = self.params['MVA_Rating']
        else:
            self.base_mva = 100
        
        # Calculate H from J (kg.m2) and pf (no of poles)
        if 'J' in self.params and 'pf' in self.params:
            self.params['H'] = 0.5 * self.params['J'] / (self.base_mva * 1e6) * (self.omega_n * 2 / self.params['pf']) ** 2
        
        # Calculate internal parameters       
        self.params['X0'] = (self.params['Xs'] + self.params['Xm'])
        self.params['Xp'] = (self.params['Xs'] + self.params['Xr'] * self.params['Xm'] / (self.params['Xr'] + self.params['Xm']))
        self.params['Xpp'] = (self.params['Xs'] + self.params['Xr'] * self.params['Xr2'] * self.params['Xm'] / (self.params['Xr'] * self.params['Xr2'] + self.params['Xm'] * self.params['Xr'] + self.params['Xm'] * self.params['Xr2']))
        self.params['T0p'] = (self.params['Xr'] + self.params['Xm']) / (self.params['Rr'])
        self.params['T0pp'] = (self.params['Xr2'] + (self.params['Xr'] * self.params['Xm']) / (self.params['Xr'] + self.params['Xm'])) / (self.params['Rr2'])
        
        # Motor start signal
        self.signals['start'] = 0
        
        # Equivalent Norton impedance for Ybus modification (NOTE: currently not used)
        self.Ym = self.params['Rs'] - 1j * self.params['Xs']
    
    def parser(self, filename):
        """
        Parse a machine file (*.mot) and populate dictionary of parameters
        """
        f = open(filename, 'r')
        
        for line in f:
            if line[0] != '#' and line.strip() != '':   # Ignore comments and blank lines
                tokens = line.strip().split('=')
                if tokens[0].strip() == 'ID':
                    self.id = tokens[1].strip()
                elif tokens[0].strip() == 'BUS_NO':
                    self.bus_no = int(tokens[1].strip())
                else:
                    self.params[tokens[0].strip()] = float(tokens[1].strip())
                
        f.close()
    
    def initialise(self, vt0, S0):
        """
        Initialise machine signals and states based on load flow voltage and complex power injection
        NOTE: currently only initialised at standstill
        """
        
        # Initialisations for locked rotor machine
        Id0 = 0
        Iq0 = 0
        Vd0 = 0
        Vq0 = 0
        Edp0 = 0
        Eqp0 = 0
        Edpp0 = 0
        Eqpp0 = 0
        slip0 = 1
        p0 = 0
        q0 = 0
        Te0 = 0
        
        """
        # Placeholder for initialising a running motor
        """

        # Initialise signals, states and parameters        
        self.signals['Id'] = Id0
        self.signals['Iq'] = Iq0
        self.signals['Vd'] = Vd0
        self.signals['Vq'] = Vq0
        self.signals['Vt'] = np.abs(vt0)
        self.signals['P'] = p0
        self.signals['Q'] = q0
        self.signals['Te'] = Te0
        self.signals['omega'] = 1 - slip0
        self.signals['Im'] = 0
        
        self.states['s'] = slip0
        self.states['Eqp'] = Eqp0
        self.states['Edp'] = Edp0
        self.states['Eqpp'] = Eqpp0
        self.states['Edpp'] = Edpp0
        
        self.check_diffs()
    
    def calc_tmech(self,s):
        """
        Calculate mechanical load torque (with a quadratic load model)
        """
        
        Tm = self.params['a'] * (1-s) ** 2
        
        return Tm
    
    def calc_currents(self,vt):
        """
        Calculate machine current injections (in network reference frame)
        """
        
        if self.signals['start'] == 1:
            # Calculate terminal voltage in dq reference frame (set to rotate with q-axis)
            Vd = -np.abs(vt) * np.sin(np.angle(vt))
            Vq = np.abs(vt) * np.cos(np.angle(vt))
            
            # Calculate Id and Iq (Norton equivalent current injection in dq frame)
            Eqpp = self.states['Eqpp']
            Edpp = self.states['Edpp']
            s = self.states['s']
            Rs = self.params['Rs']
            X0 = self.params['X0']
            Xpp = self.params['Xpp']
            T0pp = self.params['T0pp']
            
            Iq = (Rs / Xpp * (Vq - Eqpp) - Vd + Edpp) / (Xpp + Rs ** 2 / Xpp)
            Id = (Vq - Eqpp - Rs * Iq) / Xpp
            
            # Calculate power output and electrical torque
            p = -(Vd * Id + Vq * Iq)             
            q = -(Vq * Id - Vd * Iq)
            Te = (Edpp * Id + Eqpp * Iq)
            
            # Calculate machine current injection (Norton equivalent current injection in network frame)
            In = (Id + 1j * Iq) * np.exp(1j * (-np.pi / 2))
            Im = -In * self.base_mva / 100 #+ self.Ym * vt
            
            # Update signals
            self.signals['Id'] = Id
            self.signals['Iq'] = Iq
            self.signals['Vd'] = Vd
            self.signals['Vq'] = Vq
            self.signals['Te'] = Te
            self.signals['P'] = p
            self.signals['Q'] = q
            self.signals['Im'] = np.abs(Im)
            self.signals['Vt'] = np.abs(vt)
            self.signals['Vang'] = np.angle(vt)
            self.signals['omega'] =  1 - s
        
        else:
            Im = 0
        
        return Im
    
    def solve_step(self,h,dstep):
        """
        Solve machine differential equations for the next stage in the integration step
        """
        
        if self.signals['start'] == 1:
        
            # Initial state variables
            s_0 = self.states['s']
            Eqp_0 = self.states['Eqp']
            Edp_0 = self.states['Edp']
            Eqpp_0 = self.states['Eqpp']
            Edpp_0 = self.states['Edpp']
            
            Rs = self.params['Rs']
            X0 = self.params['X0']
            Xp = self.params['Xp']
            Xpp = self.params['Xpp']
            T0p = self.params['T0p']
            T0pp = self.params['T0pp']
            H = self.params['H']
            
            Id = self.signals['Id']
            Iq = self.signals['Iq']
            Te = self.signals['Te']
            
            # Electrical differential equations
            f1 = self.omega_n / np.pi * (-s_0 * Edp_0 - (Eqp_0 - (X0 - Xp) * Id) / T0p ) 
            k_Eqp = h * f1 
            
            f2 = self.omega_n / np.pi * (s_0 * Eqp_0 - (Edp_0 + (X0 - Xp) * Iq) / T0p ) 
            k_Edp = h * f2
            
            f3 = f1 + self.omega_n / np.pi * (s_0 * (Edp_0 - Edpp_0) + (Eqp_0 - Eqpp_0 + (Xp - Xpp) * Id) / T0pp ) 
            k_Eqpp = h * f3
            
            f4 = f2 + self.omega_n / np.pi * (-s_0 * (Eqp_0 - Eqpp_0) + (Edp_0 - Edpp_0 - (Xp - Xpp) * Iq) / T0pp ) 
            k_Edpp = h * f4
            
            # Mechanical equation
            Tm = self.calc_tmech(s_0)
            f5 = (Tm - Te) / (2 * H)
            k_s = h * f5
            
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
                    self.states['s'] = s_0 + k_s
                    self.dsteps['s'] = [k_s]
                elif dstep == 1:
                    # Corrector step
                    self.states['Eqp'] = Eqp_0 + 0.5 * (k_Eqp - self.dsteps['Eqp'][0])
                    self.states['Edp'] = Edp_0 + 0.5 * (k_Edp - self.dsteps['Edp'][0])
                    self.states['Eqpp'] = Eqpp_0 + 0.5 * (k_Eqpp - self.dsteps['Eqpp'][0])
                    self.states['Edpp'] = Edpp_0 + 0.5 * (k_Edpp - self.dsteps['Edpp'][0])
                    self.states['s'] = s_0 + 0.5 * (k_s - self.dsteps['s'][0]) 
            
            elif self.opt == 'runge_kutta':
                # 4th Order Runge-Kutta Method
                # Update state variables
                if dstep == 0:
                    # Save initial states
                    self.states0['s'] = s_0
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
                    self.states['s'] = s_0 + 0.5 * k_s
                    self.dsteps['s'] = [k_s]            
                elif dstep == 1:
                    self.states['Eqp'] = Eqp_0 + 0.5 * k_Eqp
                    self.dsteps['Eqp'].append(k_Eqp)
                    self.states['Edp'] = Edp_0 + 0.5 * k_Edp
                    self.dsteps['Edp'].append(k_Edp)
                    self.states['Eqpp'] = Eqpp_0 + 0.5 * k_Eqpp
                    self.dsteps['Eqpp'].append(k_Eqpp)
                    self.states['Edpp'] = Edpp_0 + 0.5 * k_Edpp
                    self.dsteps['Edpp'].append(k_Edpp)
                    self.states['s'] = s_0 + 0.5 * k_s
                    self.dsteps['s'].append(k_s)           
                elif dstep == 2:
                    self.states['Eqp'] = Eqp_0 + k_Eqp
                    self.dsteps['Eqp'].append(k_Eqp)
                    self.states['Edp'] = Edp_0 + k_Edp
                    self.dsteps['Edp'].append(k_Edp)
                    self.states['Eqpp'] = Eqpp_0 + k_Eqpp
                    self.dsteps['Eqpp'].append(k_Eqpp)
                    self.states['Edpp'] = Edpp_0 + k_Edpp
                    self.dsteps['Edpp'].append(k_Edpp)
                    self.states['s'] = s_0 + k_s
                    self.dsteps['s'].append(k_s)           
                elif dstep == 3:
                    self.states['Eqp'] = self.states0['Eqp'] + 1/6 * (self.dsteps['Eqp'][0] + 2*self.dsteps['Eqp'][1] + 2*self.dsteps['Eqp'][2] + k_Eqp)
                    self.states['Edp'] = self.states0['Edp'] + 1/6 * (self.dsteps['Edp'][0] + 2*self.dsteps['Edp'][1] + 2*self.dsteps['Edp'][2] + k_Edp)
                    self.states['Eqpp'] = self.states0['Eqpp'] + 1/6 * (self.dsteps['Eqpp'][0] + 2*self.dsteps['Eqpp'][1] + 2*self.dsteps['Eqpp'][2] + k_Eqpp)
                    self.states['Edpp'] = self.states0['Edpp'] + 1/6 * (self.dsteps['Edpp'][0] + 2*self.dsteps['Edpp'][1] + 2*self.dsteps['Edpp'][2] + k_Edpp)
                    self.states['s'] = self.states0['s'] + 1/6 * (self.dsteps['s'][0] + 2*self.dsteps['s'][1] + 2*self.dsteps['s'][2] + k_s)

    def check_diffs(self):
        """
        Check if differential equations are zero (on initialisation)
        """
    
        # State variables
        Eqp_0 = self.states['Eqp']
        Edp_0 = self.states['Edp']
        s_0 = self.states['s']
        
        Id = self.signals['Id']
        Iq = self.signals['Iq']
        Te = self.signals['Te']
        
        Rs = self.params['Rs']
        X0 = self.params['X0']
        T0p = self.params['T0p']
        Xp = self.params['Xp']
        
        dEdp = self.omega_n * s_0 * Eqp_0 - (Edp_0 + (X0 - Xp) * Iq) / T0p
        dEqp = -self.omega_n * s_0 * Edp_0 - (Eqp_0 - (X0 - Xp) * Id) / T0p
        ds = self.calc_tmech(1) - Te
        
        if round(dEdp,6) != 0 or round(dEqp,6) != 0 or round(ds,6) != 0:
            print('Warning: differential equations not zero on initialisation...')
            print('dEdp = ' + str(dEdp) + ', dEqp = ' + str(dEqp) + ', ds = ' + str(ds))