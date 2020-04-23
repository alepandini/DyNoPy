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
import fileIO.FileIO as fileIO
import utilities.hash_maps as hash_maps
import logging,collections
import numpy as np
from itertools import combinations
from _operator import itemgetter

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
        
        
    def _frequency_manager(self):
        
        self._logger.info('Calculating frequencies per aa. position from aln file...%s',self._pdbid)
        self._aligned_sequences=[];
        #self._read_a3m();
        self._read_aln();
        self._aa_freq()
    def _aa_freq(self):
        '''
            From list of aligned sequences
            iterate over each position and count freq of each aa at each position
            
        '''
        _num_of_sequences   =   len(self._aligned_sequences); # number of sequences
        fac=1.0/_num_of_sequences;                            # normalization factor

        self.PER_RESIDUE_FREQUENCIES={};
        self._list_sorted_freq=[];
        
        n_pos=len(self._aligned_sequences[0]);                # number of positions i.e. number of residues
        '''
            iterate over each position
        '''
        
        
        for i in range(n_pos):                                # position 
            self._aadict=hash_maps.aadict();
            for sequence in self._aligned_sequences:          # sequence from alignment 
                #print(sequence)
                
                j_aa=sequence[i];
                #print(j_aa,i)
                if j_aa in self._aadict:
                    self._aadict[j_aa]=self._aadict[j_aa]+fac;

            self._order_the_frequencies()
            self.PER_RESIDUE_FREQUENCIES[i+1]=self._list_sorted_freq


    def _per_pair_evolution_frequencies(self,res_pair):
        '''
            since dictionary order is not fixed 
            make a list of frequencies based on fixed aa list
            
        '''
        self._dict_aa_pair      =   hash_maps.generate_combi_aa();
        self._dict_freq_A       =   hash_maps.aadict();
        self._dict_freq_B       =   hash_maps.aadict();
        self._evolution_freq    = [];
        
        _fac=1.0/(len(self._aligned_sequences))
        
        
        for sequence in self._aligned_sequences:
            seq="";
            for i in range(len(res_pair)):
                res=res_pair[i]
                i_aa=sequence[res-1]
                seq+=i_aa
                if(i==0):
                    self._dict_freq_A[i_aa]+=1
                if(i==1):
                    self._dict_freq_B[i_aa]+=1
                    
            if("-" not in seq):
                if(seq in self._dict_aa_pair):
                    self._dict_aa_pair[seq]=self._dict_aa_pair[seq]+1
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

    def _update_perpair_score_out(self,respair,top=10):
        temp=sorted(self._dict_aa_pair.items(),key=itemgetter(1),reverse=True)
        N=len(temp);
        print(N)
        for i in range(0,top,1):
            k=temp[i][0];  v=temp[i][1];
            self._out_per_pair_scores+="%5d%5d%5s%5s%5s%12.5f\n"%(respair[0],respair[1],k[0],k[1],k,v);
        for i in range(N-top,N,1):
            k=temp[i][0];  v=temp[i][1];
            self._out_per_pair_scores+="%5d%5d%5s%5s%5s%12.5f\n"%(respair[0],respair[1],k[0],k[1],k,v);

    def calculate_per_pair_score(self,dict_params):
        self._out_per_pair_scores="%5s%5s%5s%5s%5s%12s\n"%('r_A','r_B','aa_A','aa_B','r_AB','score');
        self._get_community_list(dict_params['fcom'])
        self._get_community_dict(dict_params['fcrd'])
        for cname,npairs in self._list_community:
            self._get_top_pair_from_a_community(cname,npairs)
            for pair,jscore in self._final_list_of_pairs:
                pair=pair.split('-')
                res_pair=[int(pair[0]),int(pair[1])];
                self._per_pair_evolution_frequencies(res_pair)
                self._update_perpair_score_out(res_pair,top=dict_params['npairs'])
                self._logger.info('%-15s : %5s %5s','Processing',res_pair[0],res_pair[1])
        fName="%s.scr"%(dict_params['fout'])
        fileIO.saveFile(fName, self._out_per_pair_scores)

             
    def _order_the_frequencies(self):
        '''
            since dictionary order is not fixed 
            make a list of frequencies based on fixed aa list
            
            does not sort but uses the list with fixed order
        '''
        aalist = hash_maps.aalist();
        self._list_sorted_freq=[];
        for aa in aalist:
            self._list_sorted_freq.append(self._aadict[aa])
        
    
    def _sort_frequencies(self):
        self._logger.info('Sorting frequencies by descending order..')
        temp=sorted(self._aadict.items(),key=itemgetter(1),reverse=True)
        self._aadict={};
        for k,v in temp:
            self._aadict[k]=v;
