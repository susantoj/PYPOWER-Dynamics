#!python3
#
# Copyright (C) 2014-2015 Julius Susanto. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""
PYPOWER-Dynamics
Open Loop Test

"""

from pydyn.controller import controller
from pydyn.sym_order6a import sym_order6a
from pydyn.sym_order4 import sym_order4

from scipy.sparse.linalg import splu
import numpy as np
import matplotlib.pyplot as plt
import os, sys
    
if __name__ == '__main__':
    
    #########
    # SETUP #
    #########
    print('---------------------------------')
    print('PYPOWER-Dynamics - Open Loop Test')
    print('---------------------------------')
    
    # Open output file
    f = open('output.csv', 'w')
    
    # Program options
    dynopt = {}
    h = 0.01                # step length (s)
    t_sim = 15              # simulation time (s)
    dynopt['fn'] = 50       # Nominal system frequency (Hz)
    
    # Integrator option
    #dynopt['iopt'] = 'mod_euler'
    dynopt['iopt'] = 'runge_kutta'
    
    # Create dynamic model objects
    oCtrl = controller('smib.dyn',dynopt)
    #oMach = sym_order4('smib_round.mach',dynopt)
    oMach = sym_order6a('smib_round.mach',dynopt)     
    
    ##################
    # INITIALISATION #
    ##################
    
    print('Initialising models...')
    
    # Initialise machine
    # Get voltage and complex power injection at generator terminals
    gen_bus = 1 # Bus index for generator
    S_gen = 0
    v_gen = 1    

    # Initialise machine and grid emf
    oMach.initialise(v_gen,S_gen)   
    
    # Machine to controller interfacing
    oCtrl.signals['Vt'] = oMach.signals['Vt']
    oCtrl.signals['Vfd'] = oMach.signals['Vfd']
    
    # Initialise controller
    oCtrl.initialise()
    
    #############
    # MAIN LOOP #
    #############

    y1 = []
    t_axis = []
    vt = v_gen
    f.write('time,Vref,Vt,Vfd,Id,Iq,Vd,Vq,P,Q,Pm,omega,delta\n')
    print('Simulating...')
    for t in range(int(t_sim / h) + 1):
        if np.mod(t,1/h) == 0:
            print('t=' + str(t*h) + 's')
            
        # Controller and machine interfacing
        oMach.signals['Vfd'] = oCtrl.signals['Vfd']
        oCtrl.signals['Vt'] = oMach.signals['Vt']
        
        # Solve differential equations
        oCtrl.solve_step(h, 0)
        oMach.solve_step(h, 0)  
        oMach.signals['Vt'] = oMach.states['Eqpp']
        
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
                       
        if t == 800:
            oCtrl.signals['Vref'] = oCtrl.signals['Vref'] - 0.05
        
    plt.plot(t_axis,y1)
    plt.xlabel('Time (s)')
    plt.ylabel('GEN1:Vt (pu)')
    plt.show()
    f.close()