#!/bin/python3
'''
    @release_date  : $release_date
    @version       : $release_version
    @author        : Sarath Chandra Dantu


     Copyright (C) 2020 Sarath Chandra Dantu & Alessandro Pandini

     This file is part of: Residue ranking Python (ResPy)

     ResPy is a free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     ResPy is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with ResPy.  If not, see <http://www.gnu.org/licenses/>.

'''
'''
    python imports
'''
import os,logging,collections
'''
    generic tools
'''
import dynoutil.hash_maps as hmaps
import dynoio.fileio as fileIO
'''
    class specific tools
'''
from dynolib.residuelib import ResidueLib

class Residue(object):
    '''
        classdocs
    '''
    def __init__(self):
        '''
            Constructor
        '''
        self._logger        =   logging.getLogger('Residue')
        self._initialize();
        
    def _initialize(self):
        self._notset    =   "NOTSET"
        self._resname   =   self._notset;   
        self._resid     =   -1;   
        self._restype   =   
    def resid(self,dict_params):
        return self._resid

    def resname(self,dict_params):
        return self._resname
    
    def properties(self):
        return 0



    
