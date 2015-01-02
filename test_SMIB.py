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
Single Machine Infinite Bus (SMIB) Test

"""

from controller import controller
from interface import init_interfaces
from sym_order6 import sym_order6
from sym_order4 import sym_order4
from ext_grid import ext_grid

from mod_Ybus import mod_Ybus

from scipy.sparse.linalg import splu
import numpy as np
import matplotlib.pyplot as plt

from pypower.runpf import runpf
from pypower.ext2int import ext2int
from pypower.makeYbus import makeYbus
from pypower.loadcase import loadcase
from pypower.idx_bus import BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, \
    VM, VA, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, REF
    
if __name__ == '__main__':
    
    #########
    # SETUP #
    #########
    
    print('----------------------------')
    print('PYPOWER-Dynamics - SMIB Test')
    print('----------------------------')
    
    # Open output file
    f = open('output.csv', 'w')
    
    # Load PYPOWER case
    ppc = loadcase('smib_case.py')
    
    # Program options
    h = 0.01                # step length (s)
    t_sim = 15              # simulation time (s)
    max_err = 0.001        # Maximum error in network iteration (voltage mismatches)
    max_iter = 25           # Maximum number of network iterations
    iopt = 'mod_euler'      # integrator option
    #iopt = 'runge_kutta'
    
    # Create dynamic model objects
    oCtrl = controller('smib.dyn', iopt)
    #oMach = sym_order4('smib_round.mach', iopt)
    oMach = sym_order6('smib_round.mach', iopt)     
    oGrid = ext_grid(0.1, 0)
    
    # Create dictionary of elements
    # Hard-coded placeholder (to be replaced by a more generic loop)
    elements = {}
    elements[oCtrl.id] = oCtrl
    elements[oMach.id] = oMach
    elements['grid'] = oGrid
    
    # Make list of current injection sources (generators, external grids, etc)
    sources = []
    for source in elements.values():
        if source.__module__ in ['sym_order6', 'sym_order4', 'ext_grid']:
            sources.append(source)
    
    # Set up interfaces
    interfaces = init_interfaces(elements)
    
    ##################
    # INITIALISATION #
    ##################
    
    print('Initialising models...')
    
    # Run power flow and update bus voltages and angles in PYPOWER case object
    results, success = runpf(ppc) 
    ppc["bus"][:, VM] = results["bus"][:, VM]
    ppc["bus"][:, VA] = results["bus"][:, VA]
    
    # Build Ybus matrix
    ppc_int = ext2int(ppc)
    baseMVA, bus, branch = ppc_int["baseMVA"], ppc_int["bus"], ppc_int["branch"]
    Ybus, Yf, Yt = makeYbus(baseMVA, bus, branch)
    
    # Build modified Ybus matrix
    Ybus = mod_Ybus(Ybus, elements, bus, ppc_int['gen'], baseMVA)
    
    # Calculate initial voltage phasors
    v0 = bus[:, VM] * (np.cos(np.radians(bus[:, VA])) + 1j * np.sin(np.radians(bus[:, VA])))
    
    # Initialise sources from load flow
    for source in sources:
        source_bus = ppc_int['gen'][source.gen_no,0]
        S_source = np.complex(results["gen"][source.gen_no, 1] / baseMVA, results["gen"][source.gen_no, 2] / baseMVA)
        v_source = v0[source_bus]
        source.initialise(v_source,S_source)
    
    # Interface controllers and machines (for initialisation)
    for intf in interfaces:
        int_type = intf[0]
        var_name = intf[1]
        if int_type == 'OUTPUT':
            # If an output, interface in the reverse direction for initialisation
            intf[2].signals[var_name] = intf[3].signals[var_name]
        else:
            # Inputs are interfaced in normal direction during initialisation
            intf[3].signals[var_name] = intf[2].signals[var_name]
    
    # Initialise controller
    oCtrl.initialise()
    
    #############
    # MAIN LOOP #
    #############
    
    # Factorise Ybus matrix
    Ybus_inv = splu(Ybus)
    
    y1 = []
    t_axis = []
    v_prev = v0
    f.write('time,Vref,Vt,Vfd,Id,Iq,Vd,Vq,P,Q,Pm,omega,delta\n')
    print('Simulating...')
    for t in range(int(t_sim / h) + 1):
        if np.mod(t,1/h) == 0:
            print('t=' + str(t*h) + 's')
            
        # Interface controllers and machines
        for intf in interfaces:
            var_name = intf[1]
            intf[3].signals[var_name] = intf[2].signals[var_name]
        
        # Solve differential equations
        oCtrl.solve_step(h)
        oMach.solve_step(h)  
        
        # Solve network equations
        verr = 1
        i = 1
        # Iterate until network voltages in successive iterations are within tolerance
        while verr > max_err and i < max_iter:        
            # Update current injections for sources
            I = np.zeros(len(bus), dtype='complex')
            for source in sources:
                source_bus = ppc_int['gen'][source.gen_no,0]
                I[source_bus] = source.calc_currents(v_prev[source_bus])
            
            # Solve for network voltages
            vtmp = Ybus_inv.solve(I) 
            verr = np.abs(np.dot((vtmp-v_prev),np.transpose(vtmp-v_prev)))
            v_prev = vtmp
            i = i + 1
        
        # Signal to plot
        y1.append(oMach.signals['Vt'])
        t_axis.append(t*h)
        
        # Write to output file
        f.write(str(t*h) + ',' + str(oCtrl.signals['Vref']) + ',' + str(oMach.signals['Vt']) + ',' + \
            str(oMach.signals['Vfd']) + ',' + str(oMach.signals['Id'])  + ',' + str(oMach.signals['Iq']) \
             + ',' + str(oMach.signals['Vd']) + ',' + str(oMach.signals['Vq']) + ',' + str(oMach.signals['P']) \
             + ',' + str(oMach.signals['Q']) + ',' + str(oMach.signals['Pm']) + ',' + str(oMach.states['omega']) \
             + ',' + str(oMach.states['delta']) + '\n')
        
        # Events (placeholder)
        if t*h == 1:
            oCtrl.signals['Vref'] = oCtrl.signals['Vref'] + 0.05
            #oMach.signals['Pm'] = oMach.signals['Pm'] - 0.05

        if t*h == 8:
            oCtrl.signals['Vref'] = oCtrl.signals['Vref'] - 0.05

    plt.plot(t_axis,y1)
    plt.show()
    f.close()