#!python3
#
# Copyright (C) 2014-2015 Julius Susanto. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""
PYPOWER-Dynamics
Motor Startup Test Case

"""
# Dynamic model classes
from pydyn.asym_1cage import asym_1cage
from pydyn.ext_grid import ext_grid

# Simulation modules
from pydyn.events import events
from pydyn.recorder import recorder
from pydyn.run_sim import run_sim

# External modules
from pypower.loadcase import loadcase
import matplotlib.pyplot as plt
import numpy as np
    
if __name__ == '__main__':
    
    #########
    # SETUP #
    #########
    
    print('-----------------------------')
    print('PYPOWER-Dynamics - Motor Test')
    print('-----------------------------')

    # Load PYPOWER case
    ppc = loadcase('test_case.py')
    
    # Program options
    dynopt = {}
    dynopt['h'] = 0.01                # step length (s)
    dynopt['t_sim'] = 10.0              # simulation time (s)
    dynopt['max_err'] = 0.0001        # Maximum error in network iteration (voltage mismatches)
    dynopt['max_iter'] = 25           # Maximum number of network iterations
    dynopt['verbose'] = False         # option for verbose messages
    dynopt['fn'] = 50                 # Nominal system frequency (Hz)
    
    # Integrator option
    dynopt['iopt'] = 'mod_euler'
    #dynopt['iopt'] = 'runge_kutta'
    
    # Create dynamic model objects
    oMach = asym_1cage('motor.mach', dynopt) 
    oGrid = ext_grid('GRID1', 0, 0.1, 99999, dynopt)
    
    # Create dictionary of elements
    elements = {}
    elements[oMach.id] = oMach
    elements[oGrid.id] = oGrid
    
    # Create event stack
    oEvents = events('test_events.evnt')
    
    # Create recorder object
    oRecord = recorder('recorder.rcd')
    
    # Run simulation
    oRecord = run_sim(ppc,elements,dynopt,oEvents,oRecord)

    # Plot variables
    plt.plot(oRecord.t_axis,oRecord.results['MOT1:Vt'])
    plt.xlabel('Time (s)')
    plt.ylabel('MOT1:Vt (pu)')
    plt.show()
    
    # Write recorded variables to output file
    oRecord.write_output('output.csv')