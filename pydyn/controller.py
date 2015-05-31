#!python3
#
# Copyright (C) 2014-2015 Julius Susanto. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""
PYPOWER-Dynamics
Controller Model Class
Parses, initialises and solves a dynamic controller model file (*.dyn)

"""

import pydyn.explicit_blocks as blocks
import numpy as np

class controller:
    def __init__(self, filename, dynopt):
        self.id = ''
        self.signals = {}
        self.states = {}
        self.states0 = {}
        self.dsteps = {}
        self.equations = []
        self.init = []
        self.opt = dynopt['iopt']
        
        self.parser(filename)
    
    def parser(self, filename):
        """
        Dynamic model file parser that returns the controller ID, an ordered list of model equations (as tuples),
        a dictionary of signals and state variables initialised to zero and a list of initialisation
        equations.        
        """
        
        f = open(filename, 'r')
        
        init_flag = False
        for line in f:
            if line[0] != '#' and line.strip() != '':   # Ignore comments and blank lines
                tokens1 = line.strip().split('=')                
                if tokens1[0].strip() == 'INIT':
                    init_flag = True
                
                if init_flag == False:
                    if tokens1[0].strip() == 'ID':
                        # Controller ID
                        self.id = tokens1[1].strip()
                    else:
                        # Controller definition equations
                        equation = [tokens1[0].strip()]
                        tokens2 = tokens1[1].split('(')
                        equation.append(tokens2[0].strip())          
                        tokens3 = tokens2[1].strip(')').split(',')
                        for token in tokens3:
                            if token.strip() != '':
                                equation.append(token.strip())
                        
                        self.signals[tokens1[0].strip()] = 0
                        self.states[tokens1[0].strip()] = 0
                        self.equations.append(equation)
                
                elif tokens1[0].strip() != 'INIT':
                    # Initialisation equations
                    init = [tokens1[0].strip(), tokens1[1].strip()]
                    tokens2 = tokens1[2].split('(')
                    init.append(tokens2[0].strip())          
                    tokens3 = tokens2[1].strip(')').split(',')
                    for token in tokens3:
                        if token.strip() != '':
                            init.append(token.strip())
                    
                    self.init.append(init)
        
        f.close()
    
    def initialise(self):
        """
        Initialise controller
        """
        for line in self.init:
            type = line[0]
            var = line[1]
            block = line[2]
            
            if block == 'SUM':
                yi = self.neg_token(line[3:])
                yo = sum(yi)
                
            elif block == 'MULT':
                yi = self.neg_token(line[3:])
                yo = np.prod(yi)
                
            elif block == 'CONST':
                yo = float(line[3])
                
            if type == 'SIGNAL':
                self.signals[var] = yo
            elif type == 'STATE':
                self.states[var] = yo
    
    def solve_step(self,h,dstep):
        """
        Solve controller for the next time step
        """
        if self.opt == 'runge_kutta' and dstep in [0,1]:
            # Halve the step size for 1st and 2nd steps of 4th order Runge-Kutta method
            h = h / 2
                
        for line in self.equations:
            signal = line[0]
            block = line[1]
                       
            # Current state variable(s)
            x0 = self.states[signal]
            
            yo = None
            x1 = None

            if block == 'CONST':
                yo = float(line[2])
            
            elif block == 'GAIN':
                p = float(line[3])
                yo = blocks.gain_block(yi,p)
            
            elif block == 'INT':
                yi = self.signals[line[2]]
                p = [float(x) for x in line[3:]]
                yo, x1, f = blocks.int_block(h,x0,yi,p)
            
            elif block == 'LAG':
                yi = self.signals[line[2]]
                p = [float(x) for x in line[3:]]
                yo, x1, f = blocks.lag_block(h,x0,yi,p) 
            
            elif block == 'LDLAG':
                yi = self.signals[line[2]]
                p = [float(x) for x in line[3:]]
                yo, x1, f = blocks.leadlag_block(h,x0,yi,p)
            
            elif block == 'LIM':
                yi = self.neg_token(line[2:])
                p = [float(x) for x in line[3:]]
                yo = blocks.lim_block(yi,p)
                
            elif block == 'MULT':
                yi = self.neg_token(line[2:])
                yo = blocks.mult_block(yi)
            
            elif block == 'OUTPUT':
                self.signals[signal] = self.signals[line[2]]
            
            elif block == 'SUM':
                yi = self.neg_token(line[2:])
                yo = blocks.sum_block(yi)
            
            elif block == 'WOUT':
                yi = self.signals[line[2]]
                p = float(line[3])
                yo, x1, f = blocks.wout_block(h,x0,yi,p)                
            
            if yo:
                self.signals[signal] = yo
            
            if x1:
                if self.opt == 'mod_euler':
                    if dstep == 0:
                        self.states[signal] = x1
                        self.dsteps[signal] = [f]
                    elif dstep == 1:
                        self.states[signal] = x1 - h * self.dsteps[signal][0]
                        
                elif self.opt == 'runge_kutta':
                    if dstep == 0:
                        self.states0[signal] = x0
                        self.states[signal] = x1
                        self.dsteps[signal] = [f * 2 * h]
                    elif dstep == 1:
                        self.states[signal] = x1
                        self.dsteps[signal].append(f * 2 * h)
                    elif dstep == 2:
                        self.states[signal] = x1
                        self.dsteps[signal].append(f * h)
                    elif dstep == 3:
                        self.states[signal] = self.states0[signal] + 1/6 * (self.dsteps[signal][0] + 2*self.dsteps[signal][1] + 2*self.dsteps[signal][2] + f * h)
                
    def neg_token(self, tokens):
        """
        Consider negative sign in list of tokens
        Returns a list of token values
        """
        yi = []        
        for x in tokens:
            if x.replace('.','0').isnumeric() == True:
                yi.append(float(x))
            elif x[0] == '-':
                yi.append(-self.signals[x[1:]])
            else:
                yi.append(self.signals[x])
        
        return yi