#!python3
"""
PYPOWER-Dynamics
External Grid Model Class
Grid is modelled as a constant voltage behind a transient reactance
"""

import numpy as np
#import matplotlib.pyplot as plt
from integrators import mod_euler, runge_kutta

class ext_grid:
    def __init__(self, Xdp):
        self.Xdp = Xdp
        self.Ed = 1
        self.delta = 0
        self.Eqp = 0
        self.Vfd = 1
        self.Id = 0
        self.Iq = 0
        
    def initialise(self,vt0,S0):
        """
        Initialise grid emf based on load flow voltage and grid current injection
        """
        # Calculate initial grid current
        Ig0 = np.conj(S0 / vt0)
        phi0 = np.angle(Ig0)
        
        # Calculate steady state grid emf (i.e. voltage behind synchronous reactance)
        Eq0 = vt0 + np.complex(0,self.Xdp) * Ig0
        self.delta = np.angle(Eq0)
        
        # Convert currents to rotor reference frame
        Id0 = np.abs(Ig0) * np.sin(self.delta - phi0)
        Iq0 = np.abs(Ig0) * np.cos(self.delta - phi0)
        
        self.Vfd = np.abs(Eq0)
        self.Eqp = self.Vfd + self.Xdp * Id0
        
        Vdq = np.abs(vt0) * np.exp(np.complex(0,self.delta - np.angle(vt0) - np.pi/2))
        
        self.Ed = np.real(Vdq) + self.Xdp * Iq0
        
        Eq = np.imag(Vdq) - self.Xdp * Id0
        #self.E = v_grid + np.complex(0,self.Xdp) * i_grid
        
        self.Id = Id0
        self.Iq = Iq0
    
    def solve_step(self,h):
        p = self.Xdp
        yi = [self.Vfd, self.Id]
        f = '(yi[0] - p * yi[1] - x)'
        Eqp_1 = runge_kutta(self.Eqp,h,f,yi,p)
    
    def calc_currents(self, vt):
        """
        Solve grid equations
        """
        Vdq = np.abs(vt) * np.exp(np.complex(0,np.angle(vt) + np.pi/2 - self.delta))
        Vd = np.real(Vdq)
        Vq = np.imag(Vdq)
        
        self.Id = (self.Ed -Vd) / self.Xdp
        self.Iq = -(self.Eqp - Vq) / self.Xdp
        
        # Calculate power output
        p = Vd * self.Id + Vq * self.Iq
        q = Vq * self.Id - Vd * self.Iq
        S = np.complex(p,q)
        
        i_grid = np.conj(S / vt)
        
        return i_grid
    
    