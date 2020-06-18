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
import dynoIO.fileIO as fileIO
import dynoutil.uniutils as utils
import dynoutil.hash_maps as hmaps
import os,logging,collections

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
        print(hmaps.get_ft_types())
        self._fmt   =   "txt";
        self.MONOMER_DICT={};
        self.DICT_RES_PROP={};
        self._num_aa=0;     self._ec_id="";     self._subunit="NO PDB";
        self._resprop=collections.namedtuple('resprop', 'resid aa ss func mut');
        self._list_sequence=[];
        self._list_features=[];

    def get_uniprot_record(self):
        if(os.path.isfile(self._uni_record)==False):
            utils.get_uniprot_record(self._uniprot_id,fmt="txt");

    def process_uniprot_record(self):
        spdata=""; http_error=False
        
        if(os.path.isfile(self._uni_record)==False):
            self.get_uniprot_record();

        if(os.path.isfile(self._uni_record)==True):
            _unidata    =   fileIO.read_file(self._uni_record);
            _splituni   =   _unidata.split("\n");
            
            _full_name="";  _sname="";  _cofactor="";   
            
            self._dict_pdbs =   {};   
            
            cofac=False;    _got_name=False;
            _seq_save=False;    self._sequence="";
            for line in _splituni:
                more_split=line.split();
                if(len(more_split)>=1)and(more_split[0]!="//"):
                    if(self._num_aa==0):
                        if(more_split[0]=="ID"):
                            self._num_aa=more_split[3]
                    if(more_split[0]=="DE")and(more_split[1]=="RecName:")and(_got_name==False):
                        
                        a=more_split[2].split("=")[1];
                        _full_name=a;
                        for i in range(3,len(more_split)):
                            if("|" in more_split[i])==False:
                                _full_name+=" %s"%(more_split[i]);
                        
                        _full_name=utils.remove_grammar(_full_name);
                        _got_name=True;
                        
                    if(more_split[0]=="DE")and(more_split[1]!="RecName:"):
                        spt=more_split[1];
                        spt=spt[0:2]
                        
                        if(spt=="EC"):
                            self._ec_id=more_split[1].split("=")[1]
                            if(self._ec_id.find(";")>-1):
                                self._ec_id=utils.remove_grammar(self._ec_id)
                                    
                    if(more_split[0]=="GN")and(more_split[1][0]=="N"):
                        _sname=utils.remove_grammar(more_split[1].split("=")[1])
                    if(more_split[0]=="CC")and(len(more_split)>2):
                        if(more_split[2]=="COFACTOR:"):
                            cofac=True;
                        if(cofac==True)and(more_split[1]!="-!-"):
                            _cofactor=more_split[1].split("=")[1]
                            cofac=False;
                        if(more_split[2]=="SUBUNIT:"):
                            self._subunit=utils.remove_grammar(more_split[3])
                    if(more_split[0]=="DR")and(more_split[1]=="PDB;"):
                        _pdbid      =   utils.remove_grammar(more_split[2])
                        _method     =   utils.remove_grammar(more_split[3])
                        _res        =   0; 
                        _nchains    =   0;
                        
                        if(_method=="X-ray"):
                            _res        =   float(more_split[4]);
                            _nchains    =   utils.get_chains(more_split[6])
                        
                        if(_method=="NMR"):
                            _res        =   0.0;
                            _nchains    =   utils.get_chains(more_split[5])

                        self._dict_pdbs[_pdbid]=[_res,_method,_nchains];
                    if(more_split[0]=="FT")and(more_split[1]=="PDB;"):
                        _pdbid      =   utils.remove_grammar(more_split[2])
                        _method     =   utils.remove_grammar(more_split[3])
                        _res        =   0;
                        _nchains    =   0;

                        if(_method=="X-ray"):
                            _res        =   float(more_split[4]);
                            _nchains    =   utils.get_chains(more_split[6])

                        if(_method=="NMR"):
                            _res        =   0.0;
                            _nchains    =   utils.get_chains(more_split[5])
                    
                        self._dict_pdbs[_pdbid]=[_res,_method,_nchains];
                    if(more_split[0]=="SQ"):
                        _seq_save=True;
                    if(_seq_save==True)and(len(more_split)<=6):
                        self._sequence+="".join(more_split)
                    if(more_split[0]=="FT"):
                        '''
                            just get all the FT data
                        '''
                        sw  =   more_split[1]
                        #self._list_features.append(more_split)

        if(len(self._sequence)  ==  self._num_aa):
            self._list_sequence=list(self._sequence)
        #print(self._list_features)

    def manager(self,dict_params):
        self._uniprot_id    =   dict_params['uniprotid']
        self._uni_record    =   self._uniprot_id+"."+self._fmt;  
        self.process_uniprot_record()
