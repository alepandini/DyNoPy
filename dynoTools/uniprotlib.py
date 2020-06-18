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
import xmlschema as xmls

class UniLib(object):
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
        
        _url_schema     =   "https://www.uniprot.org/docs/uniprot.xsd";
        self._schema    =   xmls.XMLSchema(_url_schema);

        #print(hmaps.get_ft_types())
        self._fmt   =   "xml";
        self._list_sequence=[];
        self._list_features=[];

        self._resprop   =   collections.namedtuple('resprop', 'resid aa ss func mut');
        self._pdbprop   =   collections.namedtuple('pdbprop','pdbid method resolution numchains chains resrange')

        self.NUM_AA =   0;  self.RES_FIRST  =   1;  self.RES_LAST=1;
        
        self.SEQUENCE   =   "";
        
        self.PDB_LIST   =   [];
        
        self.EC_ID="";  self.ACCESSION_ID="";   self.NAME_PROTEIN="";
        

    def getuniprotdata(self,uni_file,fmt="xml"):
        self._uni_record    =   uni_file;
        self._fmt           =   fmt;
        if(self._fmt    ==  'xml'):
            self._read_xml();
        if(self._fmt    ==  'txt'):
            self._read_txt();

    def _read_xml(self):
        self._entry_dict    =   self._schema.to_dict(self._uni_record);
        '''
        for key in self._entry_dict.keys():
            print(key)
        print('^^^'*20)
        '''
        self._content =  self._entry_dict['entry'][0]
        self._keys    =  list(self._content.keys())
        
        self.ACCESSION_ID   =   self._content[self._keys[4]];
        self.NAME_PROTEIN   =   self._content[self._keys[5]];
        
        '''
            i=6 full details on name
            i=8     organism
            i=10    references for all 
        '''
        self._get_ec_id();
        self._get_sequence_details();
        self._get_pdb_list();
        self._get_residue_annotations();
        exit()
    def _get_residue_annotations(self):
        _data   =   self._get_data(i=15);
        for d in _data:
            _ft_type    =   d['@type']
            _ft_mod     =   _ft_type.replace(' ','_').upper()
            if(_ft_type ==   'chain'):
                _res_d  = d['location'];
                self.RES_FIRST  =   int(_res_d['begin']['@position']);
                self.RES_LAST   =   int(_res_d['end']['@position'])
            print(_ft_type)
        exit()
    def _get_sequence_details(self):
        '''
            
        '''
        _data           =   self._get_data(i=17);
        self.NUM_AA     =   int(_data['@length'])
        self.SEQUENCE   =   _data['$']
        self.MOLWT      =   int(_data['@mass'])
        self._logger.info('%-25s : %s'%('N. AA',self.NUM_AA))
        self._logger.info('%-25s : %s'%('Mol. Wt',self.MOLWT));
    def _get_ec_id(self):
        '''
            get the ec id
        '''
        _data       =   self._get_data(i=12);
        self.EC_ID  =   _data[0]['@id'];
        self._logger.info('%-25s : %s'%('EC ID',self.EC_ID))

    def _get_pdb_list(self):
        _db_data_list    =   self._get_data(i=12);
        for _pdb_data in _db_data_list:
            _d_type =   _pdb_data['@type']
            if(_d_type=='PDB'):
                self._save_pdb_properties(_pdb_data);
        self._logger.info('%-25s : %s'%('N. PDBs',len(self.PDB_LIST)))

    def _save_pdb_properties(self,_pdb_data):
        
        _id         =   _pdb_data['@id'];
        _method     =   "";
        _chains     =   [];
        _resolution =   0;
        _property_d =   _pdb_data['property'];
        self._dict_chains   =   {};
        for _d in _property_d:
            if(_d['@type']=='method'):
                _method =   _d['@value'];
            if(_d['@type']=='resolution'):
                _resolution =   float(_d['@value']);
            if(_d['@type']=='chains'):
                _c  = _d['@value'];
                self._get_chain_details(_c);
        
        _cl =   list(self._dict_chains.keys());
        _rl = [];
        for _c in _cl:
            _rl.append(self._dict_chains[_c])
        _d  =   self._pdbprop(pdbid=_id,method=_method,resolution=_resolution,numchains=len(_cl),chains=_cl,resrange=_rl);
        self.PDB_LIST.append(_d)
    def _get_chain_details(self,data):
        data   =   data.split(',')
        for _data in data:
            self._get_chain_properties(_data);
            
    def _get_chain_properties(self,_chain_d):
        '''
            A/B/C/D=2-161
            A/B=1-35, A/B=38-164
        '''
        _cl,_rr     =   _chain_d.split('=');
        _cl         =   _cl.split('/'); # even if only one still returns a list
        #_rr         =   _rr.split('-');
        #_rA         =   int(_rr[0]);    _rB         =   int(_rr[1]);
        mrr=[];
        for _chain in _cl:
            _i="";
            if(_chain in self._dict_chains):
                _i    =   self._dict_chains[_chain];
                mrr.append(_i);
            mrr.append(_rr);
            self._dict_chains[_chain]=mrr;
            mrr=[]; #reinitialize


    def _get_data(self,i=1):
        k=self._keys[i];
        data=self._content[k]
        return data
    def _read_txt(self):
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

        #if(len(self._sequence)  ==  self._num_aa):
        #    self._list_sequence=list(self._sequence)
        #print(self._list_features)

    def manager(self,dict_params):
        self._uniprot_id    =   dict_params['uniprotid']
        self._uni_record    =   self._uniprot_id+"."+self._fmt;  
        self.process_uniprot_record()
