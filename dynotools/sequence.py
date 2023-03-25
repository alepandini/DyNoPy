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
from Bio import pairwise2
from Bio import SeqUtils
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
        self._write_alignment = True

    def get_fasta(self,file_fasta="None"):
        _seq=""
        _dict_seq={}

        if(file_fasta=="None"):
            _seq,_dict_seq = fileIO.read_fasta(self._dict_params['file_fasta']);
        else:
            _seq,_dict_seq = fileIO.read_fasta(file_fasta);

        self.FASTA_SEQUENCE         =   _seq;
        self._logger.info('%s',_seq)
        self.DICT_FASTA_SEQUENCE    =   _dict_seq;
        self.LENGTH_SEQ             =   len(_seq);
        self._logger.info('%-25s : %d','No. of aa',self.LENGTH_SEQ);
    
    def _check_for_reference(self):
        for count,seq in enumerate(self.LIST_SEQUENCES):
            if(seq==self.FASTA_SEQUENCE):
                self._logger.info('%-25s : %d','REF_SEQ_FOUND',count+1)

    def _info_alignment(self):
        self.get_fasta()
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

    def compare_fasta_and_pdb_sequence(self,fasta_sequence="",pdb_sequence=""):
        #pass
        self._compare_pdb_fasta(fasta_sequence,pdb_sequence)
    def _compare_pdb_fasta(self,fasta_sequence,pdb_sequence,aln_file="test.aln"):
        '''
            Compare pdb to pasta
        '''
        _match=2;   _mismatch=-2;   _gap=-1;   _ext_gap=-.5;
        
        self._logger.info('####################################################################')
        self._logger.info('Performing FASTA vs. PDB sequence comparison....')
        self._logger.info('Sequence comparison parameters')
        self._logger.info('%-10s : %10d; %-10s : %10d; %-10s : %10d; %-10s : %10d','Match',_match,'Mismatch',_mismatch,'Gap',_gap,'Ext. gap',_ext_gap);
        self._logger.info('%-15s : %5d'%("Expected Max. Score",len(fasta_sequence)*_match))

        #generate alignments        
        for a in pairwise2.align.globalms(fasta_sequence,pdb_sequence,_match,_mismatch,_gap,_ext_gap):
            print(pairwise2.format_alignment(*a))
        #self._alignment_list=list(alignments[0]);
        #perfect score=perfect alignement
        #perfect_score=len(fasta_sequence)*_match
        '''
        if(self._write_alignment):
            out="";
            for a in alignments:
                out+="%s\n"%str(pairwise2.format_alignment(a))
                self._logger.debug('%s',pairwise2.format_alignment(a))
            print(pairwise2.format_alignment(alignments));
            
            fileIO.save_file(aln_file,out)
        '''
        '''
        self._logger.debug('Alignment analysis...')
        self._logger.debug('%-10s : %10d; %-10s : %10d;','Perfect score',perfect_score,'Align. score',self._alignment_list[2])
        
        pdb_residue_ids=list(self.pdbdata.RESIDUES.keys())
        
        list_of_fragments=[];   temp_fragment=""; code=99;    start_id=0;    last_id=0;
        _count_match=0; _count_gaps=0;  _count_mismatch=0;  count=0;

        for i in range(len(self.FASTA_SEQUENCE)):
            resid=i+1;
            fasta=self._alignment_list[0][i];   pdb=self._alignment_list[1][i];    
            if(pdb!="-"): #no gap i.e. either match or mis-match
                if(fasta==pdb):#match
                    code=pdb_residue_ids[count];
                    count+=1
                    temp_fragment+=pdb;
                    _count_match+=1;
                    #only store resid when temp_fragment length is 1
                    if(len(temp_fragment)==1):   
                        start_id=resid;
                else:#mis-match
                    code=1; last_id=resid;
                    list_of_fragments.append([temp_fragment,start_id,last_id]);
                    temp_fragment="";
                    _count_mismatch+=1;
            else:#gaps
                code=-1;
                _count_gaps+=1;
                if(len(temp_fragment)>0):
                    last_id=resid-1
                    
                    list_of_fragments.append([temp_fragment,start_id,last_id]);
                #re-initializing
                temp_fragment="";
            
            #if we have reached end of the list and there is still a match
            if(i==len(self.FASTA_SEQUENCE)-1):
                last_id=resid
                list_of_fragments.append([temp_fragment,start_id,last_id]);
            #
            self.RESIDUES_MASTER_LIST[resid]=[fasta,pdb,code]
        
        temp_list=[];
        for i in range(len(list_of_fragments)):
            fragment_length=len(list_of_fragments[i][0]);
            self._logger.debug('%-10s : %10d; %-10s : %10d;','Fragment ID',i,'Length',fragment_length);
            temp_list.append(fragment_length)
        
        max_frag_index=temp_list.index(numpy.max(temp_list));   max_fragment_length=temp_list[max_frag_index];
        
        self.START_ID=list_of_fragments[max_frag_index][1]; self.LAST_ID=list_of_fragments[max_frag_index][2];
        #_comparison_log+="\n\n %s \n\n"%(list_of_fragments[max_frag_index][0]);
        #self._logger.debug('PDB residues list: %s',self.pdbdata.RESIDUES)
        #self._logger.debug('Master residues list: %s',self.RESIDUES_MASTER_LIST);
        if(max_fragment_length!=len(self.pdbdata.PDB_SEQUENCE)):
            self._logger.warn('$$$ THIS PDB WILL NOT BE PROCESSED FURTHER $$$')
        if(max_fragment_length==len(self.pdbdata.PDB_SEQUENCE)):
            self._logger.info('FASTA and PDB sequence match has no issues. Proceeding ahead....')
        '''
        '''
        _aligned_seq=self._alignment_list[1];
        n_gaps=_aligned_seq.count('-')
        print(self.RESIDUES_MASTER_LIST)
        
        print("Number of Gaps %12d"%(n_gaps))
        print("Alignment score %12d (PS: %12d)"%(self._alignment_list[2],perfect_score))
        exit()
        '''
        

    def analysis_manager(self,dict_params):
        self._dict_params   =   dict_params;
        self._type_analysis =   self._dict_params['type_ana'];
        if(self._type_analysis == 0):
            self.get_fasta();
        if(self._type_analysis == 1):
            self._info_alignment();
        if(self._type_analysis == 2):
            self._process_matrix();
        if(self._type_analysis == 3):
            self._prccs();
        if(self._type_analysis == 4):
            self._calc_log_odd_matrix();
        if(self._type_analysis<0)and(self._type_analysis>4):
            self._logger.info('Choosen methods not implemented yet. Choose between 0-4')
            exit()


