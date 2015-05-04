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
Single Machine Infinite Bus (SMIB) Test

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