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

import dynolib.resmatrixlib as rmlib 
import dynoio.fileio as fileio

class ResMatrix(object):
    def __init__(self):
        self._initialize();
        
    def _initialize(self):
        
        self._serial=True;
        
        self._logger    =   logging.getLogger('Dyno ReMa')
        '''
            booleans
        '''
        self._use_coevolution       =   False
        '''
            matrices
        '''
        self._matrix_coevolution    =   [];  
        self._matrix_geometric      =   [];
        self._matrix_rho            =   [];
        self._matrix_ie             =   [];

        '''
            default values
        '''

        self._lambda_j            =   0.5;
        self._avg_c_score         =   0.0;    
        self._max_c_score         =   1.0;
        '''
            cut-offs
        '''
        self._scoe_cut            =   1.0;
        self._rho_cut             =   0.0;
        '''
            lists
        '''
        self._list_dist_data      =   []; 
        self._list_ie_data        =   [];
        
    def _get_coev_matrix(self):
        self._logger.info('%-20s : %s'%('STORING','coevolution matrix...'))
        if(self._dict_params['file_coe']==None):
            self._use_coevolution    =  False;
            self._logger.info('%-20s : %s'%('COE',self._dict_params['file_coe']))
            self._logger.info('%-20s : %s'%('MATRIX TYPE','RHO ONLY'))
        else:
            self._logger.info('%-20s : %s'%('COE',self._dict_params['file_coe']))
            self._logger.info('%-20s : %s'%('MATRIX TYPE','J'))

            self._matrix_coevolution  =   fileio.read_matrix(self._dict_params['file_coe'])
            if(self._matrix_coevolution.shape[0]>0):
                self._coe_matrix_stats();

    def _coe_matrix_stats(self):    
        self._logger.info('%-20s'%('COE PROPERTIES'))
        self._logger.info('%-20s : %d x %d'%('No. of residues',self._matrix_coevolution.shape[0],self._matrix_coevolution.shape[1]))
        self._avg_c_score         =   np.average(self._matrix_coevolution);
        self._max_c_score         =   np.max(self._matrix_coevolution);
        self._logger.info('%-20s : %.2f'%('Average',self._avg_c_score));
        self._logger.info('%-20s : %.2f'%('Maximum',self._max_c_score));

    def _calc_j_score(self,scaled_coe,rho):
        global scoe_cut,rho_cut,lambda_j
        jvalue=0.0; rho=np.abs(rho);
        
        if(scaled_coe>=scoe_cut)and(rho>=rho_cut):
            jvalue  =   lambda_j*(scaled_coe/max_c_score)+(1-lambda_j)*rho;
        return jvalue

    def _get_geometric_data(self):
        '''
            get the data in PCA/distance vector 
        '''
        self._logger.info('%-20s : %s'%('Processing', 'geometrical data...'))
        self._matrix_geometric  =   fileio.read_data_to_matrix(self._dict_params['file_gem']);
        #self._matrix_geometric  =   fileio.read_matrix(self._dict_params['file_gem']);
        self._logger.info('%-20s : %s'%('No. of GEM vectors',self._matrix_geometric.shape))
    
    def _calculate_scaled_c(c_value):
        global avg_c_score
        scaled_c=0;
        if(c_value>0):
            scaled_c=c_value/avg_c_score;
        return scaled_c

    def _correlation_manager(self):
        
        #out_data="%5s%5s%12s%12s%12s%12s\n"%('Res_a','Res_b','Cab','Sab','Rab','j');
        out_data="%8s%8s%12s%12s%12s\n"%('Res_a','Res_b','Cab','Sab','Rab');
        if(self._use_coevolution    ==  False):
            self._logger.info('%-20s : %s'%('Calculating...','Rho Matrix'))
        else:
            self._logger.info('%-20s : %s'%('Calculating...','J Matrix'))
        
        '''
            generate the list of tuples with residue pairs (i,j)
        '''

        _res_first  =   self._dict_params['resi_fst'];
        _res_last   =   self._dict_params['resi_lst'];
        _list_of_pairs  =   [];
        for i in range(_res_first,_res_last+1,1):
            for j in range(_res_first,_res_last+1,1):
                if(j>i):
                    _list_of_pairs.append((i,j))
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
        for data in self._list_of_correlations:
            _i  =   data[0]-1;  _j  =   data[1]-1
            _matrix_correlations[_i][_j]  =   data[2];
            _matrix_correlations[_j][_i]  =   data[2];
        
        rhotype="PEA";
        if(self._dict_params['corr_met']    ==  1):
            rhotype="SPM"
        elif(self._dict_params['corr_met']    ==  2):
            rhotype="NMI"
        
        flabel="%s-%s.mat"%(rhotype,self._dict_params['file_lab'])
        fileio.save_matrix(flabel,_matrix_correlations)

    def manager(self,dict_params):
        _start   =   timeit.default_timer();
        
        self._dict_params   =   dict_params;

        self._get_geometric_data();
        self._get_coev_matrix();
        self._correlation_manager();
        _stop    =   timeit.default_timer();
        self._logger.info('%-20s : %.2f (s)'%('FINISHED_IN',_stop-_start));
