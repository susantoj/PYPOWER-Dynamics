#!python3
#
# Copyright (C) 2014-2015 Julius Susanto. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

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