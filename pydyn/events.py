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
Events Class
Sets up and handles events in the simulation
"""

import numpy as np
from pypower.idx_bus import BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, \
    VM, VA, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, REF

class events:
    def __init__(self, filename):
        self.event_stack = []
        self.parser(filename) 
            
    def parser(self, filename):
        """
        Parse an event file (*.evnt) and populate event stack
        """
        f = open(filename, 'r')
        
        for line in f:
            if line[0] != '#' and line.strip() != '':   # Ignore comments and blank lines
                tokens = line.strip().split(',')
                
                # Parse signal events
                if tokens[1].strip() in ['SIGNAL', 'FAULT', 'LOAD', 'STATE']:
                    self.event_stack.append([float(tokens[0].strip()), tokens[1].strip(), tokens[2].strip(), tokens[3].strip(), tokens[4].strip()])
                
                elif tokens[1].strip() in ['CLEAR_FAULT', 'TRIP_BRANCH']:
                    self.event_stack.append([float(tokens[0].strip()), tokens[1].strip(), tokens[2].strip()])
                    
        f.close()
        
    def handle_events(self, t, elements, ppc, baseMVA):
        """
        Checks and handles the event stack during a simulation time step
        """
        refactorise = False
        
        if self.event_stack:
            if self.event_stack[0][0] < t:
                print('Event missed at t=' + str(self.event_stack[0][0]) + 's... Check simulation time step!')
                del self.event_stack[0]
            
            # Event exists at time step
            while self.event_stack and self.event_stack[0][0] == t:
                event_type = self.event_stack[0][1]
                
                # Handle signal events
                if event_type == 'SIGNAL':
                    obj_id = self.event_stack[0][2]
                    sig_id = self.event_stack[0][3]
                    value = float(self.event_stack[0][4])
                    elements[obj_id].signals[sig_id] = value
                    
                    print('SIGNAL event at t=' + str(t) + 's on element "' + obj_id + '". ' + sig_id + ' = ' + str(value) + '.')
                
                if event_type == 'STATE':
                    obj_id = self.event_stack[0][2]
                    sig_id = self.event_stack[0][3]
                    value = float(self.event_stack[0][4])
                    elements[obj_id].states[sig_id] = value
                    
                    print('STATE event at t=' + str(t) + 's on element "' + obj_id + '". ' + sig_id + ' = ' + str(value) + '.')
                
                if event_type == 'FAULT':
                    bus_id = int(self.event_stack[0][2])
                    Rf = float(self.event_stack[0][3])
                    Xf = float(self.event_stack[0][4])
                    
                    if Rf == 0:
                        ppc["bus"][bus_id, GS] = 1e6
                    elif Rf < 0:
                        ppc["bus"][bus_id, GS] = 0
                        Rf = 'Inf'
                    else:
                        ppc["bus"][bus_id, GS] = 1 / Rf * baseMVA
                    
                    if Xf == 0:
                        ppc["bus"][bus_id, BS] = -1e6
                    elif Xf < 0:
                        ppc["bus"][bus_id, BS] = 0
                        Xf = 'Inf'
                    else:
                        ppc["bus"][bus_id, BS] = -1 / Xf * baseMVA
                    
                    refactorise = True
                    
                    print('FAULT event at t=' + str(t) + 's on bus at row "' + str(bus_id) + '" with fault impedance Zf = ' + str(Rf) + ' + j' + str(Xf) + ' pu.')
                
                if event_type == 'CLEAR_FAULT':
                    bus_id = int(self.event_stack[0][2])
                    ppc["bus"][bus_id, BS] = 0
                    ppc["bus"][bus_id, GS] = 0
                    refactorise = True
                    
                    print('CLEAR_FAULT event at t=' + str(t) + 's on bus at row "' + str(bus_id) + '".')
                
                if event_type == 'TRIP_BRANCH':
                    branch_id = int(self.event_stack[0][2])
                    ppc["branch"] = np.delete(ppc["branch"],branch_id, 0)
                    refactorise = True
                    
                    print('TRIP_BRANCH event at t=' + str(t) + 's on branch "' + str(branch_id) + '".')
                
                if event_type == 'LOAD':
                    bus_id = int(self.event_stack[0][2])
                    Pl = float(self.event_stack[0][3])
                    Ql = float(self.event_stack[0][4])
                    
                    ppc["bus"][bus_id, PD] = Pl
                    ppc["bus"][bus_id, QD] = Ql
                    
                    refactorise = True
                    
                    print('LOAD event at t=' + str(t) + 's on bus at row "' + str(bus_id) + '" with S = ' + str(Pl) + ' MW + j' + str(Ql) + ' MVAr.')
                    
                del self.event_stack[0]
                
        return ppc, refactorise