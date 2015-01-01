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
Controller Model Class
Parses, initialises and solves a dynamic controller model file (*.dyn)

"""

import explicit_blocks as blocks
import numpy as np

class controller:
    def __init__(self, filename, iopt):
        self.id = ''
        self.signals = {}
        self.states = {}
        self.equations = []
        self.init = []
        self.opt = iopt
        
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
                
            if type == 'SIGNAL':
                self.signals[var] = yo
            elif type == 'STATE':
                self.states[var] = yo
    
    def solve_step(self,h):
        """
        Solve controller for the next time step
        """
        for line in self.equations:
            signal = line[0]
            block = line[1]
            
            # Current state variable(s)
            x0 = self.states[signal]
                        
            if block == 'LAG':
                yi = self.signals[line[2]]
                p = [float(x) for x in line[3:]]
                yo, x1 = blocks.lag_block(h,x0,yi,p,self.opt)
                self.states[signal] = x1
                self.signals[signal] = yo
                
            elif block == 'INT':
                yi = self.signals[line[2]]
                p = [float(x) for x in line[3:]]
                yo, x1 = blocks.int_block(h,x0,yi,p,self.opt)
                self.states[signal] = x1
                self.signals[signal] = yo
            
            elif block == 'LDLAG':
                yi = self.signals[line[2]]
                p = [float(x) for x in line[3:]]
                yo, x1 = blocks.leadlag_block(h,x0,yi,p,self.opt)
                self.states[signal] = x1
                self.signals[signal] = yo
            
            elif block == 'GAIN':
                yo = blocks.gain_block(yi,p)
                self.signals[signal] = yo
            
            elif block == 'SUM':
                yi = self.neg_token(line[2:])
                yo = blocks.sum_block(yi)
                self.signals[signal] = yo
            
            elif block == 'MULT':
                yi = self.neg_token(line[2:])
                yo = blocks.mult_block(yi)
                self.signals[signal] = yo
    
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