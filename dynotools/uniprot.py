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
import dynoio.fileio as fileIO
import dynoutil.uniutils as utils
import dynoutil.hash_maps as hmaps
import os,logging,collections

from dynolib.uniprotlib import UniLib
class Uniprot(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
            Constructor
        '''
        self._logger        =   logging.getLogger('Dyno UP')
        self._initialize();
        
    def _initialize(self):
        self._uniprot_id="";

        self._unilib=UniLib();
        
        #print(hmaps.get_ft_types())
        self._fmt   =   "txt";
        self.MONOMER_DICT={};
        self.DICT_RES_PROP={};
        self._num_aa=0;     self._ec_id="";     self._subunit="NO PDB";
        self._resprop=collections.namedtuple('resprop', 'resid aa ss func mut');
        self._list_sequence=[];
        self._list_features=[];
        self._fmt   =   'xml';

    def get_uniprot_record(self):
        if(os.path.isfile(self._uni_record)==False):
            utils.get_uniprot_record(self._uniprot_id,fmt=self._fmt);

    def process_uniprot_record(self):
        spdata=""; http_error=False
        if(os.path.isfile(self._uni_record)==False):
            self.get_uniprot_record();

        if(os.path.isfile(self._uni_record)==True):
            self._unilib.getuniprotdata(self._uni_record,self._fmt);

    def manager(self,dict_params):
        self._uniprot_id    =   dict_params['uniprotid']
        self._uni_record    =   self._uniprot_id+"."+self._fmt;  
        self.process_uniprot_record()
