#!python3
#
# Copyright (C) 2014-2015 Julius Susanto. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

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