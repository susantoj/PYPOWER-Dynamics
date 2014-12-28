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
Open Loop Test

"""

from controller import controller
from sym_order6 import sym_order6
from sym_order4 import sym_order4
from ext_grid import ext_grid

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
    
    # Create dynamic model objects
    oCtrl = controller('smib.dyn')
    #oMach = sym_order4('smib_round.mach')
    oMach = sym_order6('smib_round.mach')     
    
    ##################
    # INITIALISATION #
    ##################
    
    print('Initialising models...')
    
    # Initialise machine
    # Get voltage and complex power injection at generator terminals
    gen_bus = 1 # Bus index for generator (placeholder)
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
    h = 0.01
    vt = v_gen
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
    plt.show()
    f.close()