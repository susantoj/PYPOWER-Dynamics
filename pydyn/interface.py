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
Functions for the interface between controller and machine variables

"""

def init_interfaces(elements):
    ints_list = []
    
    for element in elements.values():
        if element.__module__ == 'pydyn.controller':
            for line in element.equations:
                if line[1] == 'INPUT':
                    new_int = [line[1],line[2],elements[line[3]],element]
                    ints_list.append(new_int)
                
                if line[1] == 'OUTPUT':
                    new_int = [line[1],line[0],element,elements[line[3]]]
                    ints_list.append(new_int)
    
    return ints_list