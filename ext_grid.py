#!python3
#
# Copyright (C) 2014 Julius Susanto
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
and with infinite inertia.
"""

import numpy as np
from integrators import mod_euler, runge_kutta

class ext_grid:
    def __init__(self, Xdp):
        self.Xdp = Xdp
        self.E = 1
        
    def initialise(self,vt0,S0):
        """
        Initialise grid emf based on load flow voltage and grid current injection
        """
        # Calculate initial grid current
        Ig0 = np.conj(S0 / vt0)       
        self.E = vt0 + np.complex(0,self.Xdp) * Ig0
    
    def calc_currents(self, vt):
        """
        Solve grid equations
        """
        
        i_grid = self.E * np.exp(-1j * np.angle(vt)) / np.complex(0,self.Xdp)
        
        return i_grid
    
    