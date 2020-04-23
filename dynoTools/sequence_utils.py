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
import dynoIO.fileIO as fileIO
import dynoutil.hash_maps as hash_maps
import logging,collections
import numpy as np
from itertools import combinations
from _operator import itemgetter

class SeqTools(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
            Constructor
        '''
        self._logger=logging.getLogger('Dyno SEQ');
        self._initialize()

    def _initialize(self):
        self.LIST_FIXED_AA              =   hash_maps.aalist();
        self.DICT_RESIDUE_FREQUENCIES   =   {};
        self.NUM_OF_SEQUENCES           =   0;
        self.NUM_OF_RESIDUES            =   0;

    def per_residue_aa_frequencies(self):
        '''
            From the list of aligned sequences
            iterate over each position and count freq of each aa at each position
        '''
        
        self.NUM_OF_SEQUENCES   =   len(self._aligned_sequences);
        fac=1.0/self.NUM_OF_SEQUENCES;                            

        self.DICT_RESIDUE_FREQUENCIES   =   {};
        self._list_sorted_freq          =   [];
        
        self.NUM_OF_RESIDUES            =   len(self._aligned_sequences[0]);
        '''
            iterate over each position
        '''
        for i in range(self.NUM_OF_RESIDUES):
            _aadict =   hash_maps.aadict();
            for sequence in list_sequences: 
                _aa_at_pos_i    =   sequence[i];
                if _aa_at_pos_i in _aadict:
                    _aadict[_aa_at_pos_i]=_aadict[_aa_at_pos_i]+fac;
            # re-order the dict of frequencies as per the FIXED_AA_LIST
            self._order_the_frequencies()
            self.DICT_RESIDUE_FREQUENCIES[i+1]=self._list_sorted_freq

    def per_pair_log_odds_score(self,res_A,res_B):

        self._dict_aa_pair      =   hash_maps.generate_combi_aa();
        self._freq_A       =   ;
        self._dict_freq_B       =   hash_maps.aadict();
        self._evolution_freq    = [];
        
        self._calculate_per_pair_scores();
        
    def _calculate_per_pair_scores(self):
        aa_pair_list=list(self._dict_aa_pair.keys());
        for aa_pair in aa_pair_list:
            aa_A=aa_pair[0];                aa_B=aa_pair[1];
            f_A=self._dict_freq_A[aa_A];    n_A=hash_maps.aafreq_from_literature(aa_A);
            f_B=self._dict_freq_B[aa_B];    n_B=hash_maps.aafreq_from_literature(aa_B);
            numerator   =   f_A*f_B
            denominator =   n_A*n_B
            score=0;
            if(numerator>0)and(denominator>0):
                    score=np.log(numerator/denominator)
            self._dict_aa_pair[aa_pair]=score;


    def calc_log_odd_matrix(self,aligned_sequence):
        #calculate the frequencies of aa for each residue
        self.per_residue_aa_frequencies();

        for res_i in range(self.NUM_OF_RESIDUES):
            for res_j in range(self.NUM_OF_RESIDUES):
                if(res_i!=res_j):
                    self.ppef(res_i,res_j)
                self._logger.info('%-15s : %5s %5s','Processing',res_pair[0],res_pair[1])
        fName="%s.scr"%(dict_params['fout'])
        fileIO.saveFile(fName, self._out_per_pair_scores)
    def _order_the_frequencies(self):
        '''
            since dictionary order is not fixed 
            make a list of frequencies based on fixed aa list
            
            does not sort but uses the list with fixed order
        '''
        self._list_sorted_freq=[];
        for aa in self.LIST_FIXED_AA:
            self._list_sorted_freq.append(self._aadict[aa])
    def _sort_frequencies(self):
        self._logger.info('Sorting frequencies by descending order..')
        temp=sorted(self._aadict.items(),key=itemgetter(1),reverse=True)
        self._aadict={};
        for k,v in temp:
            self._aadict[k]=v;
