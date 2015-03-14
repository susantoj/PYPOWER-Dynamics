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
from pydyn.controller import controller
from pydyn.sym_order6a import sym_order6a
from pydyn.sym_order4 import sym_order4
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
    
    print('----------------------------')
    print('PYPOWER-Dynamics - SMIB Test')
    print('----------------------------')

    # Load PYPOWER case
    ppc = loadcase('smib_case.py')
    
    # Program options
    dynopt = {}
    dynopt['h'] = 0.01                # step length (s)
    dynopt['t_sim'] = 15              # simulation time (s)
    dynopt['max_err'] = 0.0001        # Maximum error in network iteration (voltage mismatches)
    dynopt['max_iter'] = 25           # Maximum number of network iterations
    dynopt['verbose'] = False         # option for verbose messages
    dynopt['fn'] = 50                 # Nominal system frequency (Hz)
    
    # Integrator option
    #dynopt['iopt'] = 'mod_euler'
    dynopt['iopt'] = 'runge_kutta'
    
    # Create dynamic model objects
    oCtrl = controller('smib.dyn', dynopt)
    oMach = sym_order6a('smib_round.mach', dynopt)
    #oMach = sym_order4('smib_round.mach', iopt) 
    oGrid = ext_grid('GRID1', 0, 0.1, 99999, dynopt)
    
    # Create dictionary of elements
    # Hard-coded placeholder (to be replaced by a more generic loop)
    elements = {}
    elements[oCtrl.id] = oCtrl
    elements[oMach.id] = oMach
    elements[oGrid.id] = oGrid
    
    # Create event stack
    oEvents = events('smib_events.evnt')
    
    # Create recorder object
    oRecord = recorder('smib_recorder.rcd')
    
    # Run simulation
    oRecord = run_sim(ppc,elements,dynopt,oEvents,oRecord)
    
    # Calculate relative rotor angles
    rel_delta = np.array(oRecord.results['GEN1:delta']) - np.array(oRecord.results['GRID1:delta'])
    
    # Plot variables
    #plt.plot(oRecord.t_axis,rel_delta)
    plt.plot(oRecord.t_axis,oRecord.results['GEN1:Vt'])
    plt.xlabel('Time (s)')
    plt.ylabel('GEN1:Vt (pu)')
    plt.show()
    
    # Write recorded variables to output file
    oRecord.write_output('output.csv')