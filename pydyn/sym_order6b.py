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
Based on Sauer-Pai model
Sauer, P.W., Pai, M. A., "Power System Dynamics and Stability", Stipes Publishing, 2006 
"""

import numpy as np

class sym_order6b:
    def __init__(self, filename, dynopt):
        self.id = ''
        self.gen_no = 0
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
            self.params['Xa'] = self.params['Xa'] * 100 / base_mva
            self.params['Xd'] = self.params['Xd'] * 100 / base_mva
            self.params['Xdp'] = self.params['Xdp'] * 100 / base_mva
            self.params['Xdpp'] = self.params['Xdpp'] * 100 / base_mva
            self.params['Xq'] = self.params['Xq'] * 100 / base_mva
            self.params['Xqp'] = self.params['Xqp'] * 100 / base_mva
            self.params['Xqpp'] = self.params['Xqpp'] * 100 / base_mva
        
        # Internal variables
        self.params['gamma_d1'] = (self.params['Xdpp'] - self.params['Xa']) / (self.params['Xdp'] - self.params['Xa'])
        self.params['gamma_d2'] = (1 - self.params['gamma_d1']) / (self.params['Xdp'] - self.params['Xa'])        
        self.params['gamma_q1'] = (self.params['Xqpp'] - self.params['Xa']) / (self.params['Xqp'] - self.params['Xa'])
        self.params['gamma_q2'] = (1 - self.params['gamma_q1']) / (self.params['Xqp'] - self.params['Xa'])
        
        # Equivalent Norton impedance for Ybus modification
        self.Yg = (self.params['Ra'] - 1j * 0.5 * (self.params['Xdpp'] + self.params['Xqpp'])) / (self.params['Ra'] **2 + (self.params['Xdpp'] * self.params['Xqpp']))
    
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
        
        Vd0 = np.abs(vt0) * np.sin(delta0 - np.angle(vt0))
        Vq0 = np.abs(vt0) * np.cos(delta0 - np.angle(vt0))
        
        # Calculate machine state variables and Vfd
        Ra = self.params['Ra']
        Xa = self.params['Xa']
        Xd = self.params['Xd']
        Xdp = self.params['Xdp']
        Xdpp = self.params['Xdpp']
        Xqp = self.params['Xqp']
        Xqpp = self.params['Xqpp']     
        gamma_d1 = self.params['gamma_d1']
        gamma_d2 = self.params['gamma_d2']
        gamma_q1 = self.params['gamma_q1']
        gamma_q2 = self.params['gamma_q2']
        
        Edp0 = Vd0 - Xqpp * Iq0 + Ra * Id0 - (1 - gamma_q1) * (Xqp - Xa) * Iq0
        Eqp0 = Vq0 + Xdpp * Id0 + Ra * Iq0 + (1 - gamma_d1) * (Xdp - Xa) * Id0
        phid_pp0 = Eqp0 - (Xdp - Xa) * Id0
        phiq_pp0 = -Edp0 - (Xqp - Xa) * Iq0
        Vfd0 = Eqp0 + (Xd - Xdp) * (Id0 - gamma_d2 * phid_pp0 - (1 - gamma_d1) * Id0 + gamma_d2 * Eqp0)        
        
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
        self.signals['Tm'] = p0
        
        self.states['omega'] = 1
        self.states['delta'] = delta0
        self.states['Eqp'] = Eqp0
        self.states['phiq_pp'] = phiq_pp0
        self.states['Edp'] = Edp0
        self.states['phid_pp'] = phid_pp0
        
        self.check_diffs()
        
    def check_diffs(self):
        """
        Check if differential equations are zero (on initialisation)
        """
    
        # Initial state variables
        Eqp_0 = self.states['Eqp']
        Edp_0 = self.states['Edp']
        phid_pp_0 = self.states['phid_pp']
        phiq_pp_0 = self.states['phiq_pp']
        
        Xa = self.params['Xa']
        Xd = self.params['Xd']
        Xdp = self.params['Xdp']
        Xdpp = self.params['Xdpp']
        Xq = self.params['Xq']
        Xqp = self.params['Xqp']
        Xqpp = self.params['Xqpp']
        Td0p = self.params['Td0p']
        Td0pp = self.params['Td0pp']
        Tq0p = self.params['Tq0p']
        Tq0pp = self.params['Tq0pp']
        
        gamma_d1 = self.params['gamma_d1']
        gamma_d2 = self.params['gamma_d2']
        gamma_q1 = self.params['gamma_q1']
        gamma_q2 = self.params['gamma_q2']
        
        Vfd = self.signals['Vfd']
        Id = self.signals['Id']
        Iq = self.signals['Iq']
        
        dEqp = (Vfd - (Xd - Xdp) * (Id - gamma_d2 * phid_pp_0 - (1 - gamma_d1) * Id + gamma_d2 * Eqp_0) - Eqp_0) / Td0p
        dEdp = ((Xq - Xqp) * (Iq - gamma_q2 * phiq_pp_0 - (1 - gamma_q1) * Iq - gamma_q2 * Edp_0) - Edp_0) / Tq0p
        dphid_pp = (Eqp_0 - (Xdp - Xa) * Id - phid_pp_0) / Td0pp
        dphiq_pp = (-Edp_0 - (Xqp - Xa) * Iq - phiq_pp_0) / Tq0pp
        
        if round(dEdp,6) != 0 or round(dEqp,6) != 0 or round(dphid_pp,6) != 0 or round(dphiq_pp,6) != 0:
            print('Warning: differential equations not zero on initialisation...')
            print('dEdp = ' + str(dEdp) + ', dEqp = ' + str(dEqp) + ', dphid_pp = ' + str(dphid_pp) + ', dphiq_pp = ' + str(dphiq_pp))
    
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
        phid_pp = self.states['phid_pp']
        phiq_pp = self.states['phiq_pp']
        
        Ra = self.params['Ra']
        Xa = self.params['Xa']
        Xdp = self.params['Xdp']
        Xdpp = self.params['Xdpp']
        Xqp = self.params['Xqp']
        Xqpp = self.params['Xqpp']     
        gamma_d1 = self.params['gamma_d1']
        gamma_q1 = self.params['gamma_q1']
        
        # Check if speed-voltage term should be included
        if self.speed_volt:
            omega = self.states['omega']
        else:
            omega = 1
        
        Id = (-Vq / omega + gamma_d1 * Eqp + (1 - gamma_d1) * phid_pp - Ra / (omega * Xqpp) * (Vd - gamma_q1 * Edp + (1 - gamma_q1) * phiq_pp)) / (Xdpp + Ra ** 2 / (omega * omega * Xqpp))
        Iq = (Vd / omega + (Ra * Id / omega) - gamma_q1 * Edp + (1 - gamma_q1) * phiq_pp) / Xqpp
        
        # Calculate power output
        p = (Vd + self.params['Ra']*Id) * Id + (Vq  + self.params['Ra']*Iq) * Iq
        q = Vq * Id - Vd * Iq
        S = np.complex(p,q)
        
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
        
    def solve_step(self,h,dstep):
        """
        Solve machine differential equations for the next stage in the integration step
        """
        
        # Initial state variables
        omega_0 = self.states['omega']
        delta_0 = self.states['delta']
        Eqp_0 = self.states['Eqp']
        Edp_0 = self.states['Edp']
        phid_pp_0 = self.states['phid_pp']
        phiq_pp_0 = self.states['phiq_pp']
        
        Xa = self.params['Xa']
        Xd = self.params['Xd']
        Xdp = self.params['Xdp']
        Xdpp = self.params['Xdpp']
        Xq = self.params['Xq']
        Xqp = self.params['Xqp']
        Xqpp = self.params['Xqpp']
        Td0p = self.params['Td0p']
        Td0pp = self.params['Td0pp']
        Tq0p = self.params['Tq0p']
        Tq0pp = self.params['Tq0pp']
        
        gamma_d1 = self.params['gamma_d1']
        gamma_d2 = self.params['gamma_d2']
        gamma_q1 = self.params['gamma_q1']
        gamma_q2 = self.params['gamma_q2']
                
        Vfd = self.signals['Vfd']
        Id = self.signals['Id']
        Iq = self.signals['Iq']
        
        # Electrical differential equations
        f1 = (Vfd - (Xd - Xdp) * (Id - gamma_d2 * phid_pp_0 - (1 - gamma_d1) * Id + gamma_d2 * Eqp_0) - Eqp_0) / Td0p
        k_Eqp = h * f1
        
        f2 = ((Xq - Xqp) * (Iq - gamma_q2 * phiq_pp_0 - (1 - gamma_q1) * Iq - gamma_q2 * Edp_0) - Edp_0) / Tq0p
        k_Edp = h * f2
        
        f3 = (Eqp_0 - (Xdp - Xa) * Id - phid_pp_0) / Td0pp
        k_phid_pp = h * f3
        
        f4 = (-Edp_0 - (Xqp - Xa) * Iq - phiq_pp_0) / Tq0pp
        k_phiq_pp = h * f4
        
        # Swing equation
        f5 = 0.5 / self.params['H'] * (self.signals['Pm'] / omega_0 - self.signals['P'])
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
                self.states['phiq_pp'] = phiq_pp_0 + k_phiq_pp
                self.dsteps['phiq_pp'] = [k_phiq_pp]
                self.states['phid_pp'] = phid_pp_0 + k_phid_pp
                self.dsteps['phid_pp'] = [k_phid_pp]
                self.states['omega'] = omega_0 + k_omega
                self.dsteps['omega'] = [k_omega]
                self.states['delta'] = delta_0 + k_delta
                self.dsteps['delta'] = [k_delta]
            elif dstep == 1:
                # Corrector step
                self.states['Eqp'] = Eqp_0 + 0.5 * (k_Eqp - self.dsteps['Eqp'][0])
                self.states['Edp'] = Edp_0 + 0.5 * (k_Edp - self.dsteps['Edp'][0])
                self.states['phiq_pp'] = phiq_pp_0 + 0.5 * (k_phiq_pp - self.dsteps['phiq_pp'][0])
                self.states['phid_pp'] = phid_pp_0 + 0.5 * (k_phid_pp - self.dsteps['phid_pp'][0])
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
                self.states0['phiq_pp'] = phiq_pp_0 
                self.states0['phid_pp'] = phid_pp_0
                
                self.states['Eqp'] = Eqp_0 + 0.5 * k_Eqp
                self.dsteps['Eqp'] = [k_Eqp]
                self.states['Edp'] = Edp_0 + 0.5 * k_Edp
                self.dsteps['Edp'] = [k_Edp]
                self.states['phiq_pp'] = phiq_pp_0 + 0.5 * k_phiq_pp
                self.dsteps['phiq_pp'] = [k_phiq_pp]
                self.states['phid_pp'] = phid_pp_0 + 0.5 * k_phid_pp
                self.dsteps['phid_pp'] = [k_phid_pp]
                self.states['omega'] = omega_0 + 0.5 * k_omega
                self.dsteps['omega'] = [k_omega]            
                self.states['delta'] = delta_0 + 0.5 * k_delta
                self.dsteps['delta'] = [k_delta]
            elif dstep == 1:
                self.states['Eqp'] = Eqp_0 + 0.5 * k_Eqp
                self.dsteps['Eqp'].append(k_Eqp)
                self.states['Edp'] = Edp_0 + 0.5 * k_Edp
                self.dsteps['Edp'].append(k_Edp)
                self.states['phiq_pp'] = phiq_pp_0 + 0.5 * k_phiq_pp
                self.dsteps['phiq_pp'].append(k_phiq_pp)
                self.states['phid_pp'] = phid_pp_0 + 0.5 * k_phid_pp
                self.dsteps['phid_pp'].append(k_phid_pp)
                self.states['omega'] = omega_0 + 0.5 * k_omega
                self.dsteps['omega'].append(k_omega)           
                self.states['delta'] = delta_0 + 0.5 * k_delta
                self.dsteps['delta'].append(k_delta)
            elif dstep == 2:
                self.states['Eqp'] = Eqp_0 + k_Eqp
                self.dsteps['Eqp'].append(k_Eqp)
                self.states['Edp'] = Edp_0 + k_Edp
                self.dsteps['Edp'].append(k_Edp)
                self.states['phiq_pp'] = phiq_pp_0 + k_phiq_pp
                self.dsteps['phiq_pp'].append(k_phiq_pp)
                self.states['phid_pp'] = phid_pp_0 + k_phid_pp
                self.dsteps['phid_pp'].append(k_phid_pp)
                self.states['omega'] = omega_0 + k_omega
                self.dsteps['omega'].append(k_omega)           
                self.states['delta'] = delta_0 + k_delta
                self.dsteps['delta'].append(k_delta)
            elif dstep == 3:
                self.states['Eqp'] = self.states0['Eqp'] + 1/6 * (self.dsteps['Eqp'][0] + 2*self.dsteps['Eqp'][1] + 2*self.dsteps['Eqp'][2] + k_Eqp)
                self.states['Edp'] = self.states0['Edp'] + 1/6 * (self.dsteps['Edp'][0] + 2*self.dsteps['Edp'][1] + 2*self.dsteps['Edp'][2] + k_Edp)
                self.states['phiq_pp'] = self.states0['phiq_pp'] + 1/6 * (self.dsteps['phiq_pp'][0] + 2*self.dsteps['phiq_pp'][1] + 2*self.dsteps['phiq_pp'][2] + k_phiq_pp)
                self.states['phid_pp'] = self.states0['phid_pp'] + 1/6 * (self.dsteps['phid_pp'][0] + 2*self.dsteps['phid_pp'][1] + 2*self.dsteps['phid_pp'][2] + k_phid_pp)
                self.states['omega'] = self.states0['omega'] + 1/6 * (self.dsteps['omega'][0] + 2*self.dsteps['omega'][1] + 2*self.dsteps['omega'][2] + k_omega)
                self.states['delta'] = self.states0['delta'] + 1/6 * (self.dsteps['delta'][0] + 2*self.dsteps['delta'][1] + 2*self.dsteps['delta'][2] + k_delta)
                self.signals['Tm'] = self.signals['Pm'] / omega_0
