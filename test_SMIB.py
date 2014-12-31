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
from sym_order6 import sym_order6
from sym_order4 import sym_order4
from ext_grid import ext_grid

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
    
    # Integrator option
    iopt = 'mod_euler'
    #iopt = 'runge_kutta'
    
    # Create dynamic model objects
    oCtrl = controller('smib.dyn', iopt)
    #oMach = sym_order4('smib_round.mach', iopt)
    oMach = sym_order6('smib_round.mach', iopt)     
    oGrid = ext_grid(0.01)
    
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
    
    # Calculate initial voltage phasors
    v0 = results["bus"][:, VM] * (np.cos(np.radians(results["bus"][:, VA])) + 1j * np.sin(np.radians(results["bus"][:, VA])))
    
    # Set up augmented Ybus matrix
    #y_gen = np.complex(oMach.params['Ra'], -0.5 * (oMach.params['Xdpp'] + oMach.params['Xqpp'])) / \
    #        np.complex(oMach.params['Ra'] ** 2, oMach.params['Xdpp'] * oMach.params['Xqpp'])
    #y_grid = np.complex(oGrid.params['Ra'], -0.5 * (oGrid.params['Xdpp'] + oGrid.params['Xqpp'])) / \
    #        np.complex(oGrid.params['Ra'] ** 2, oGrid.params['Xdpp'] * oGrid.params['Xqpp'])
    Ybus[0,0] = Ybus[0,0] + 1 / (1j * oGrid.Xdp)
    Ybus[1,1] = Ybus[1,1] + 1 / (oMach.params['Ra'] + 1j * 0.5 * (oMach.params['Xdpp'] + oMach.params['Xqpp']))
    
    # Add equivalent load admittance to Ybus matrix
    Pl, Ql = ppc["bus"][:, PD], ppc["bus"][:, QD]
    for i in range(len(Pl)):
        S_load = np.complex(Pl[i],Ql[i])
        y_load = S_load / results["bus"][i, VM] ** 2
        Ybus[i,i] = Ybus[i,i] + y_load
    
    # Calculate grid current
    grid_bus = 0    # Bus index for grid (placeholder)
    S_grid = np.complex(results["gen"][grid_bus, 1] / baseMVA, results["gen"][grid_bus, 2] / baseMVA)
    v_ang = np.radians(bus[grid_bus, VA])
    v_grid = bus[grid_bus, VM] * np.complex(np.cos(v_ang), np.sin(v_ang)) 

    # Initialise machine from load flow
    # Get voltage and complex power injection at generator terminals
    gen_bus = 1 # Bus index for generator (placeholder)
    S_gen = np.complex(results["gen"][gen_bus, 1] / baseMVA, results["gen"][gen_bus, 2] / baseMVA)
    v_ang = np.radians(bus[gen_bus, VA])
    v_gen = v0[1] #bus[gen_bus, VM] * np.complex(np.cos(v_ang), np.sin(v_ang))    
    
    # Initialise machine and grid emf
    oMach.initialise(v_gen,S_gen) 
    oGrid.initialise(v_grid,S_grid)     
    
    # Machine to controller interfacing
    oCtrl.signals['Vt'] = oMach.signals['Vt']
    oCtrl.signals['Vfd'] = oMach.signals['Vfd']
    
    # Initialise controller
    oCtrl.initialise()
    
    #############
    # MAIN LOOP #
    #############
    
    # Factorise Ybus matrix
    Ybus_inv = splu(Ybus)
    
    y1 = []
    t_axis = []
    h = 0.01
    vt = v_gen
    vg = v_grid
    v_prev = v0
    f.write('time,Vref,Vt,Vfd,Id,Iq,Vd,Vq,P,Q,Pm,omega,delta\n')
    print('Simulating...')
    for t in range(1501):
        if np.mod(t,100) == 0:
            print('t=' + str(t*h) + 's')
            
        # Controller and machine interfacing
        oMach.signals['Vfd'] = oCtrl.signals['Vfd']
        oCtrl.signals['Vt'] = oMach.signals['Vt']
        
        # Solve differential equations
        oCtrl.solve_step(h)
        oMach.solve_step(h)  
        
        # Solve network equations
        verr = 1
        i = 1
        # Iterate until network voltages in successive iterations are within tolerance
        while verr > 0.001 and i <25:        
            # Update generator current injections
            Im = oMach.calc_currents(vt) 
            Ig = oGrid.calc_currents(vg)
            I = np.array([Ig, Im])
            #print(I)
            
            # Solve for network voltages
            vtmp = Ybus_inv.solve(I) 
            #print(vtmp)
            verr = np.abs(np.dot((vtmp-v_prev),np.transpose(vtmp-v_prev)))
            #print(verr)
            v_prev = vtmp
            vt = vtmp[1]
            vg = vtmp[0]
            i = i + 1
        
        #y1.append(Ig)
        y1.append(oMach.signals['Vt'])
        #y1.append(oCtrl.signals['Vfd'])
        t_axis.append(t*h)
        
        # Write to output file
        f.write(str(t*h) + ',' + str(oCtrl.signals['Vref']) + ',' + str(oMach.signals['Vt']) + ',' + \
            str(oMach.signals['Vfd']) + ',' + str(oMach.signals['Id'])  + ',' + str(oMach.signals['Iq']) \
             + ',' + str(oMach.signals['Vd']) + ',' + str(oMach.signals['Vq']) + ',' + str(oMach.signals['P']) \
             + ',' + str(oMach.signals['Q']) + ',' + str(oMach.signals['Pm']) + ',' + str(oMach.states['omega']) \
             + ',' + str(oMach.states['delta']) + '\n')
        
        # Events
        if t == 100:
            oCtrl.signals['Vref'] = oCtrl.signals['Vref'] + 0.05
            #oMach.signals['Pm'] = oMach.signals['Pm'] - 0.05
            pass
            
        if t == 800:
            oCtrl.signals['Vref'] = oCtrl.signals['Vref'] - 0.05
            pass
        
    plt.plot(t_axis,y1)
    plt.show()
    f.close()