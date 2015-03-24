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
Build modified Ybus matrix

"""

import numpy as np
from pypower.idx_bus import BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, \
    VM, VA, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, REF


def mod_Ybus(Ybus, elements, bus, gen, baseMVA):
    # Add equivalent generator and grid admittances to Ybus matrix
    for element in elements.values():
        Ye = 0
        
        # 4th/6th order machines and converters
        if element.__module__ in ['pydyn.sym_order4', 'pydyn.sym_order6a', 'pydyn.sym_order6b', 'pydyn.vsc_average']:
            i = gen[element.gen_no,0]
            Ye = element.Yg
        
        # External grid
        if element.__module__ == 'pydyn.ext_grid':
            i = gen[element.gen_no,0]
            Ye = 1 / (1j * element.params['Xdp'])
        
        if Ye != 0:
            Ybus[i,i] = Ybus[i,i] + Ye

    # Add equivalent load admittance to Ybus matrix    
    Pl, Ql = bus[:, PD], bus[:, QD]
    for i in range(len(Pl)):
        S_load = (Pl[i] - 1j * Ql[i]) / baseMVA
        y_load = S_load / (bus[i, VM] ** 2)
        Ybus[i,i] = Ybus[i,i] + y_load
    
    return Ybus