#!/usr/bin/python3
'''
    @release_date  : $release_date
    @version       : $release_version
    @author        : Sarath Chandra Dantu
    

     Copyright (C) 2020 Sarath Chandra Dantu & Alessandro Pandini

     This file is part of: Dynamics based Network cOmparisons in Python (DyNoPy)

     DyNoPy is a free software: you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.

     DyNoPy is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with DyNoPy.  If not, see <http://www.gnu.org/licenses/>.

'''
import logging
import dynoio.fileio as fileIO
from dynolib.sequencelib import SeqTools
import dynoutil.dynoMath as dMath
import numpy as np

class Sequence(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
            Constructor
        '''
        self._logger=logging.getLogger('Dyno SEQ');
        self._logger.info('%s','-'*100)
        self._initialize();

    def _initialize(self):
        self._default_N         =   9;
        self.FASTA_SEQUENCE     =   ""; self.DICT_FASTA_SEQUENCE={};    self.LENGTH_SEQ=0;
        self.LIST_SEQUENCES     =   [];
        self.MATRIX_COEV        =   np.zeros((self._default_N,self._default_N));
        self.MATRIX_COEV_RHO    =   np.zeros((self._default_N,self._default_N));
        self._seq_utils         = SeqTools();
    def _get_fasta(self):
        _seq,_dict_seq              =   fileIO.read_fasta(self._dict_params['file_fasta']);
        self.FASTA_SEQUENCE         =   _seq;
        self.DICT_FASTA_SEQUENCE    =   _dict_seq;
        self.LENGTH_SEQ             =   len(_seq);
        self._logger.info('%-25s : %d','No. of aa',self.LENGTH_SEQ);
    
    def _check_for_reference(self):
        for count,seq in enumerate(self.LIST_SEQUENCES):
            if(seq==self.FASTA_SEQUENCE):
                self._logger.info('%-25s : %d','REF_SEQ_FOUND',count+1)

    def _info_alignment(self):
        self._get_fasta()
        self._get_alignment();
        self._check_for_reference()
        
    def _get_alignment(self):
        self.LIST_SEQUENCES         =   fileIO.read_aln(self._dict_params['file_aln']);

    def _coevolution_matrix_analysis(self):
        _num_of_residues    = self.MATRIX_COEV_RHO.shape[0];
        _out                =   "%6s%6s%12s%12s%12s\n"%("res_a","res_b","coe","scoe","nscoe");
        _avg                =   1.0;
        _max                =   np.max(self.MATRIX_COEV_RHO);
        _list_ana           =   [];
        for i in range(_num_of_residues):
            for j in range(_num_of_residues):
                if(i>j):
                    _coev=self.MATRIX_COEV[i,j];    _rho=self.MATRIX_COEV_RHO[i,j]; _nrho=0;
                    if(_rho>0):
                        _nrho=_rho/_max;
                    _list_ana.append((i+1,j+1,_coev,_rho,_nrho));
                    _out+="%6d%6d%12.3f%12.3f%12.3f\n"%(i+1,j+1,_coev,_rho,_nrho);
        _fname="PairsList-%s.txt"%(self._dict_params['file_lab'])
        fileIO.save_file(_fname,_out)
    def _process_matrix(self):
        self.MATRIX_COEV            =   fileIO.read_matrix(self._dict_params['file_mat']);
        self.MATRIX_COEV_RHO        =   dMath.normalize_matrix(self.MATRIX_COEV);
        self._coevolution_matrix_analysis();
        self._logger.info('%-25s : %s','MATRIX_SHAPE',self.MATRIX_COEV.shape)
        _fname="RHO-%s.mat"%(self._dict_params['file_lab'])
        fileIO.save_matrix(_fname,self.MATRIX_COEV_RHO)

    def _prccs(self):
        self._get_matrix();
        _prccs_data=dMath.calc_prccs(self.MATRIX_COEV_RHO,method=2);
        _fname="PRCCS-%s.txt"%(self._dict_params['file_lab'])
        fileIO.save_file(_fname,_prccs_data)

    def _calc_log_odd_matrix(self):
        fName="%s.h5"%(self._dict_params['file_lab']);
        self._info_alignment();
        self._get_matrix();
        self._seq_utils.calc_log_odd_matrix(self.LIST_SEQUENCES,fName);

    def analysis_manager(self,dict_params):
        self._dict_params   =   dict_params;
        self._type_analysis =   self._dict_params['type_ana'];
        if(self._type_analysis  ==   0):
            self._get_fasta();
        if(self._type_analysis  ==  1):
            self._info_alignment();
        if(self._type_analysis  ==  2):
            self._process_matrix();
        if(self._type_analysis  ==  3):
            self._prccs();
        if(self._type_analysis  == 4):
            self._calc_log_odd_matrix();
        if(self._type_analysis<0)and(self._type_analysis>4):
            self._logger.info('Choosen methods not implemented yet. Choose between 0-4')
            exit()

