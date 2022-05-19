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
import os,timeit,logging
import numpy as np
import dynolib.pwielib as pwielib
import dynolib.resmatrixlib as rmlib 
import dynoio.fileio as fileio

class ResMatrix(object):
    def __init__(self):
        self._initialize();
        
    def _initialize(self):
        self._logger    =   logging.getLogger('Dyno ReMa')
        '''
            matrices
        '''
        self._matrix_geometric      =   [];
        
    def _get_geometric_data(self):
        '''
            get the data in PCA/distance vector 
        '''
        self._logger.info('%-20s : %s'%('Processing', 'geometrical data...'))
        self._matrix_geometric  =   fileio.read_data_to_matrix(self._dict_params['file_gem']);
        #self._matrix_geometric  =   fileio.read_matrix(self._dict_params['file_gem']);
        self._logger.info('%-20s : %s'%('No. of GEM vectors',self._matrix_geometric.shape))
    
    def _correlation_manager(self):
        self._logger.info('%-20s : %s'%('Calculating...','Rho-Matrix'))
        
        '''
            generates a list of tuples with residue pairs (i,j). Each thread will get a residue tuple and a dictionary

        '''
        _res_first  =   self._dict_params['resi_fst'];
        _res_last   =   self._dict_params['resi_lst'];
        _list_of_pairs  =   pwielib.gen_pair_list(_res_first,_res_last,gap=1)
        #for i in range(_res_first,_res_last+1,1):
        #    for j in range(_res_first,_res_last+1,1):
        #        if(j>i):
        #            _list_of_pairs.append((i,j))

        _corr_params    =   {};
        _corr_params    =   {1 :   _list_of_pairs,
                             2 :   self._matrix_geometric,
                             3 :   self._dict_params['num_rep'],
                             4 :   self._dict_params['file_lab'],
                             5 :   self._dict_params['num_thr'],
                             6 :   self._dict_params['corr_met'],
                             7 :   self._dict_params['fold_iex'],
                             8 :   self._dict_params['corr_vec']
                             }

        self._list_of_correlations  =   rmlib.correlation_calculator(_corr_params);
        self._save_list_to_matrix();

    def _save_list_to_matrix(self):
        _matrix_correlations    =   np.zeros((self._dict_params['resi_lst'],self._dict_params['resi_lst']));
        _out_data="#%11s%12s"%('resA','resB');
        for i in range(len(self._list_of_correlations[0][2])):
            _s  =   "wvec_%d"%(i+1);
            _out_data+="%12s"%(_s)
        _out_data+="\n"
        for data in self._list_of_correlations:
            _out_data+="%12d%12d"%(data[0],data[1]);
            for i in data[2]:
                _out_data+="%12.5f"%(i);
            _out_data+="\n"
            #_matrix_correlations[_i][_j]  =   data[2];
            #_matrix_correlations[_j][_i]  =   data[2];
        rhotype="PEA";
        if(self._dict_params['corr_met']    ==  1):
            rhotype="SPM"
        elif(self._dict_params['corr_met']    ==  2):
            rhotype="NMI"
        
        fname="%s-%s.txt"%(rhotype,self._dict_params['file_lab'])
        #fileio.save_matrix(flabel,_matrix_correlations)
        fileio.save_file(fname,_out_data)

    def manager(self,dict_params):
        _start   =   timeit.default_timer();
        
        self._dict_params   =   dict_params;
        self._get_geometric_data();
        self._correlation_manager();
        _stop    =   timeit.default_timer();
        self._logger.info('%-20s : %.2f (s)'%('FINISHED_IN',_stop-_start));
