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

import dynolib.networkslib as nwlib 
import dynoio.fileio as fileio

class Networks(object):
    def __init__(self):
        self._initialize();
        
    def _initialize(self):
        
        self._logger    =   logging.getLogger('Dyno ReMa')
        '''
            booleans
        '''
        self._use_coevolution       =   False
        '''
            matrices
        '''
        self._matrix_coevolution    =   [];  
        self._matrix_rho            =   [];

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
        if(self._dict_params['file_coe']!=None):
            self._matrix_coevolution  =   fileio.read_matrix(self._dict_params['file_coe'])
            if(self._matrix_coevolution.shape[0]>0):
                self._matrix_stats(self._matrix_coevolution);

    def _get_rho_matrix(self):
        self._logger.info('%-20s : %s'%('STORING','RHO matrix...'))
        if(self._dict_params['file_rho']!=None):
            self._matrix_rho  =   fileio.read_matrix(self._dict_params['file_rho'])
            if(self._matrix_rho.shape[0]>0):
                self._matrix_stats(self._matrix_rho);

    def _matrix_stats(self,_matrix):    
        self._logger.info('%-20s'%('MATRIX PROPERTIES'))
        self._logger.info('%-20s : %d x %d'%('No. of residues',_matrix.shape[0],_matrix.shape[1]))
        _avg    =   np.average(_matrix);
        _max    =   np.max(_matrix);
        self._logger.info('%-20s : %.2f'%('Average',_avg));
        self._logger.info('%-20s : %.2f'%('Maximum',_max));
    def _check_matrix_size(self):
        N_coe   =   self._matrix_coevolution.shape[0];
        N_rho   =   self._matrix_rho[-1][1];
        if(N_rho != N_coe):
            self._logger.error('%-20s : N (COE): %5d , N (RHO): %5d'%("RESIDUE_NUM_MISMATCH_ERROR",N_coe,N_rho))
            exit()
    def _network_manager(self):
        self._matrix_j  =   [];
        if(self._dict_params['file_coe']!=None)and(self._dict_params['file_rho']!=None):
            self._check_matrix_size();
            nvectors    =   self._matrix_rho.shape[1]-2;
            self._logger.info('CALCULATING NETWORK PROPERTIES')
            for i in range(nvectors):
                vec     =   i+2; 
                label   =   "%s-v%d"%(self._dict_params['file_lab'],vec);
                self._logger.info('........ for vector : %5d'%(vec))
                self._matrix_j  =   nwlib.calculate_jmatrix(self._matrix_rho,self._matrix_coevolution,rcutoff=0.5,coecutoff=1.0,vector_num=vec);
                nwlib.calculate_network_properties(self._matrix_j,label);
        
        if(self._dict_params['file_coe']==None)and(self._dict_params['file_rho']!=None):
            self._logger.info('CALCULATING NETWORK PROPERTIES')
            nwlib.calculate_network_properties(self._matrix_rho,self._dict_params['file_lab']);


    def manager(self,dict_params):
        _start   =   timeit.default_timer();
        self._dict_params   =   dict_params;
        self._get_coev_matrix();
        self._get_rho_matrix();
        self._network_manager();
        _stop    =   timeit.default_timer();
        self._logger.info('%-20s : %.2f (s)'%('FINISHED_IN',_stop-_start));
