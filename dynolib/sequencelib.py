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
import dynoio.fileio as fileIO
import dynoutil.hash_maps as hash_maps
import logging,collections,h5py,tqdm
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
        self._logger.info('%-15s : %s','Calculating','Per Residue AA frequency')
        self.NUM_OF_SEQUENCES   =   len(self._aligned_sequences);
        fac=1.0/self.NUM_OF_SEQUENCES;                            

        self.DICT_RESIDUE_FREQUENCIES   =   {};
        self._list_sorted_freq          =   [];
        
        self.NUM_OF_RESIDUES            =   len(self._aligned_sequences[0]);
        '''
            iterate over each position
        '''
        _pbar   = tqdm.tqdm(total=100,desc="Dyno SEQ - INFO - AA Freq @ Residue",position=0)
        for i in range(self.NUM_OF_RESIDUES):
            self._aadict =   hash_maps.aadict();
            for sequence in self._aligned_sequences: 
                _aa_at_pos_i    =   sequence[i];
                if _aa_at_pos_i in self._aadict:
                    self._aadict[_aa_at_pos_i]=self._aadict[_aa_at_pos_i]+fac;
            # re-order the dict of frequencies as per the FIXED_AA_LIST
            self._order_the_frequencies()
            self.DICT_RESIDUE_FREQUENCIES[i+1]=self._list_sorted_freq
            _pbar.update(100/452)

    def per_pair_log_odds_score(self,res_A,res_B):
        self._dict_joint_freq   =   hash_maps.generate_combi_aa();
        self._freq_A            =   self.DICT_RESIDUE_FREQUENCIES[res_A];
        self._freq_B            =   self.DICT_RESIDUE_FREQUENCIES[res_B];
        self._n_fixed_aa        =   len(self.LIST_FIXED_AA);
        self._matrix_log        =   np.zeros((self._n_fixed_aa,self._n_fixed_aa))
        self._matrix_frq        =   np.zeros((self._n_fixed_aa,self._n_fixed_aa))
        self._tup_aa_frq        =   [];
        self._calculate_per_pair_scores();
        _ds_name                =   "%d-%d"%(res_A,res_B);

        self._tup_aa_frq=sorted(self._tup_aa_frq,key=lambda tup: tup[1],reverse=True)

        self._h5_logOdds.create_dataset(_ds_name, data=self._matrix_log);       
        self._h5_frq.create_dataset(_ds_name,data=np.array(self._tup_aa_frq,dtype="S"));

    def _calculate_per_pair_scores(self):
        for i in range(self._n_fixed_aa):
            aa_A    =   self.LIST_FIXED_AA[i];
            f_A     =   self._freq_A[i];
            n_A     =   hash_maps.aafreq_from_literature(aa_A);

            for j in range(self._n_fixed_aa):
                aa_B        =   self.LIST_FIXED_AA[j];
                f_B         =   self._freq_B[j];    
                n_B         =   hash_maps.aafreq_from_literature(aa_B);
                numerator   =   f_A*f_B
                denominator =   n_A*n_B
                score       =   0;
                aa_pair     =   aa_A+aa_B;
                if(numerator>0)and(denominator>0):
                    score=np.log(numerator/denominator)
                
                self._dict_joint_freq[aa_pair]  =   numerator;
                self._matrix_frq[i][j]          =   numerator;
                self._matrix_log[i][j]          =   score;
                self._tup_aa_frq.append((aa_pair,score));
    def calc_log_odd_matrix(self,aligned_sequence,fName):
        #calculate the frequencies of aa for each residue
        self._aligned_sequences =   aligned_sequence;
        
        self._h5_logOdds        = h5py.File('LogODD-'+fName, 'w')
        self._h5_frq            = h5py.File('EVFRQ-'+fName, 'w')

        self.per_residue_aa_frequencies();
        _pbar   = tqdm.tqdm(total=100,desc="Dyno SEQ - INFO - Log ODD @ Residue",position=0)
        for res_i in range(self.NUM_OF_RESIDUES):
            res_i+=1;
            for res_j in range(self.NUM_OF_RESIDUES):
                res_j+=1;
                #self._logger.info('%-15s : %5s %5s','Processing',res_i,res_j)
                if(res_i!=res_j):
                    self.per_pair_log_odds_score(res_i,res_j)
            _pbar.update(100.0/self.NUM_OF_RESIDUES);
        self._h5_logOdds.close();
        self._h5_frq.close();
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
