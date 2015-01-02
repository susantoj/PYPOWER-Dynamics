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
# Dynamic model classes
from controller import controller
from sym_order6 import sym_order6
from sym_order4 import sym_order4
from ext_grid import ext_grid

# Simulation modules
from events import events
from recorder import recorder
from run_sim import run_sim

# External modules
from pypower.loadcase import loadcase
import matplotlib.pyplot as plt
    
if __name__ == '__main__':
    
    #########
    # SETUP #
    #########
    
    print('----------------------------')
    print('PYPOWER-Dynamics - SMIB Test')
    print('----------------------------')
    
    # Open output file
    #f = open('output.csv', 'w')
    
    # Load PYPOWER case
    ppc = loadcase('smib_case.py')
    
    # Integrator option
    iopt = 'mod_euler'      
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
    
    # Create event stack
    oEvents = events('smib_events.evnt')
    
    # Create recorder object
    oRecord = recorder('smib_recorder.rcd')
    
    # Run simulation
    oRecord = run_sim(ppc,elements,oEvents,oRecord)
    
    # Plot variables
    plt.plot(oRecord.t_axis,oRecord.results['GEN1:Vt'])
    plt.show()
    
    # Write recorded variables to output file
    #oRecord.write_output('output.csv')