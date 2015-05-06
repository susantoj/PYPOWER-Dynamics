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
Single Cage Asynchronous Machine Model

Model equations based on section 15.2.4 of:
Milano, F., "Power System Modelling and Scripting", Springer-Verlag, 2010

"""

import numpy as np

class asym_1cage:
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
            self.params['H'] = self.params['H'] * self.base_mva / 100
            self.params['a'] = self.params['a'] * self.base_mva / 100
            self.params['Rs'] = self.params['Rs'] * 100 / self.base_mva
            self.params['Xs'] = self.params['Xs'] * 100 / self.base_mva
            self.params['Xm'] = self.params['Xm'] * 100 / self.base_mva
            self.params['Rr'] = self.params['Rr'] * 100 / self.base_mva
            self.params['Xr'] = self.params['Xr'] * 100 / self.base_mva
        else:
            self.base_mva = 100
        
        # Calculate internal parameters       
        self.params['X0'] = self.params['Xs'] + self.params['Xm']
        self.params['Xp'] = self.params['Xs'] + self.params['Xr'] * self.params['Xm'] / (self.params['Xr'] + self.params['Xm'])
        self.params['T0p'] = (self.params['Xr'] + self.params['Xm']) / (self.omega_n * self.params['Rr'])
        
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
        slip0 = 1
        p0 = 0
        q0 = 0
        Te0 = 0
        
        """
        # Placeholder code for initialising a running motor
        Rs = self.params['Rs']
        X0 = self.params['X0']
        T0p = self.params['T0p']
        Xp = self.params['Xp']
        
        # Calculate initial armature current
        Ia0 =  np.conj(S0 / vt0)
        phi0 = np.angle(Ia0)
        
        # Convert currents to rotor reference frame
        Id0 = np.abs(Ia0) * np.sin(-np.pi/2 - phi0)
        Iq0 = np.abs(Ia0) * np.cos(-np.pi/2 - phi0)
        
        Vd0 = -np.abs(vt0) * np.sin(np.angle(vt0))
        Vq0 = np.abs(vt0) * np.cos(np.angle(vt0))
        
        # Calculate active and reactive power
        p0 = -(Vd0 * Id0 + Vq0 * Iq0)             
        q0 = -(Vq0 * Id0 - Vd0 * Iq0)
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
        
        self.states['s'] = slip0
        self.states['Eqp'] = Eqp0
        self.states['Edp'] = Edp0
        
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
            Eqp = self.states['Eqp']
            Edp = self.states['Edp']
            s = self.states['s']
            Rs = self.params['Rs']
            X0 = self.params['X0']
            T0p = self.params['T0p']
            Xp = self.params['Xp']
            
            Iq = (Rs / Xp * (Vq - Eqp) - Vd + Edp) / (Xp + Rs ** 2 / Xp)
            Id = (Vq - Eqp - Rs * Iq) / Xp
            
            # Calculate power output and electrical torque
            p = -(Vd * Id + Vq * Iq)             
            q = -(Vq * Id - Vd * Iq)
            Te = (Edp * Id + Eqp * Iq) #/ self.omega_n
            
            # Calculate machine current injection (Norton equivalent current injection in network frame)
            In = (Id + 1j * Iq) * np.exp(1j * (-np.pi / 2))
            Im = -In #+ self.Ym * vt
            
            # Update signals
            self.signals['Id'] = Id
            self.signals['Iq'] = Iq
            self.signals['Vd'] = Vd
            self.signals['Vq'] = Vq
            self.signals['Te'] = Te
            self.signals['P'] = p
            self.signals['Q'] = q
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
            
            Rs = self.params['Rs']
            X0 = self.params['X0']
            T0p = self.params['T0p']
            Xp = self.params['Xp']
            H = self.params['H']
            
            Id = self.signals['Id']
            Iq = self.signals['Iq']
            Te = self.signals['Te']
            
            # Electrical differential equations
            f1 = (-self.omega_n * s_0 * Edp_0 - (Eqp_0 - (X0 - Xp) * Id) / T0p ) * self.base_mva / 100
            k_Eqp = h * f1
            
            f2 = (self.omega_n * s_0 * Eqp_0 - (Edp_0 + (X0 - Xp) * Iq) / T0p ) * self.base_mva / 100
            k_Edp = h * f2
            
            # Mechanical equation
            Tm = self.calc_tmech(s_0)
            f3 = (Tm - Te) / (2 * H)
            k_s = h * f3
            
            if self.opt == 'mod_euler':
                # Modified Euler
                # Update state variables
                if dstep == 0:
                    # Predictor step
                    self.states['Eqp'] = Eqp_0 + k_Eqp
                    self.dsteps['Eqp'] = [k_Eqp]
                    self.states['Edp'] = Edp_0 + k_Edp
                    self.dsteps['Edp'] = [k_Edp]
                    self.states['s'] = s_0 + k_s
                    self.dsteps['s'] = [k_s]
                elif dstep == 1:
                    # Corrector step
                    self.states['Eqp'] = Eqp_0 + 0.5 * (k_Eqp - self.dsteps['Eqp'][0])
                    self.states['Edp'] = Edp_0 + 0.5 * (k_Edp - self.dsteps['Edp'][0])
                    self.states['s'] = s_0 + 0.5 * (k_s - self.dsteps['s'][0]) 
            
            elif self.opt == 'runge_kutta':
                # 4th Order Runge-Kutta Method
                # Update state variables
                if dstep == 0:
                    # Save initial states
                    self.states0['s'] = s_0
                    self.states0['Eqp'] = Eqp_0 
                    self.states0['Edp'] = Edp_0
                    
                    self.states['Eqp'] = Eqp_0 + 0.5 * k_Eqp
                    self.dsteps['Eqp'] = [k_Eqp]
                    self.states['Edp'] = Edp_0 + 0.5 * k_Edp
                    self.dsteps['Edp'] = [k_Edp]
                    self.states['s'] = s_0 + 0.5 * k_s
                    self.dsteps['s'] = [k_s]            
                elif dstep == 1:
                    self.states['Eqp'] = Eqp_0 + 0.5 * k_Eqp
                    self.dsteps['Eqp'].append(k_Eqp)
                    self.states['Edp'] = Edp_0 + 0.5 * k_Edp
                    self.dsteps['Edp'].append(k_Edp)
                    self.states['s'] = s_0 + 0.5 * k_s
                    self.dsteps['s'].append(k_s)           
                elif dstep == 2:
                    self.states['Eqp'] = Eqp_0 + k_Eqp
                    self.dsteps['Eqp'].append(k_Eqp)
                    self.states['Edp'] = Edp_0 + k_Edp
                    self.dsteps['Edp'].append(k_Edp)
                    self.states['s'] = s_0 + k_s
                    self.dsteps['s'].append(k_s)           
                elif dstep == 3:
                    self.states['Eqp'] = self.states0['Eqp'] + 1/6 * (self.dsteps['Eqp'][0] + 2*self.dsteps['Eqp'][1] + 2*self.dsteps['Eqp'][2] + k_Eqp)
                    self.states['Edp'] = self.states0['Edp'] + 1/6 * (self.dsteps['Edp'][0] + 2*self.dsteps['Edp'][1] + 2*self.dsteps['Edp'][2] + k_Edp)
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