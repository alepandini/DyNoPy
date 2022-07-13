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

class JMatrix(object):
    def __init__(self):
        self._initialize();
        
    def _initialize(self):
        
        self._serial=True;
        
        self._logger    =   logging.getLogger('Dyno JMAT')
        '''
            booleans
        '''
        self._use_coevolution       =   False
        '''
            matrices
        '''
        self._matrix_coevolution    =   [];  
        self._matrix_rho            =   [];
        self._dict_rho_data         =   {};
        self._dict_J_data           =   {};


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
        
    def _get_coev_matrix(self):
        self._logger.info('%-20s : %s'%('STORING','coevolution matrix...'))
        self._matrix_coevolution  =   fileio.read_matrix(self._dict_params['file_coe'])
        if(self._matrix_coevolution.shape[0]>0):
            self._coe_matrix_stats();
        
        self._calculate_scaled_matrix();

    def _get_rho_matrix(self):
        self._logger.info('%-20s : %s'%('STORING','dynamics matrix...'))
        self._max_scaled_coe=np.max(self._scaled_matrix)

        self._out_j="%12s%12s%12s%12s"%('Res_a','Res_b','Cab','Sab');
        #self._out_coe="%8s%8s%12s%12s%12s\n"%('Res_a','Res_b','Cab','Sab','Rab');

        #self._matrix_rho  =   fileio.read_matrix(self._dict_params['file_rho'])
        self._rho_data=open(self._dict_params['file_rho']).readlines()

    def _calculate_j_score(self):
        self._logger.info('%-20s'%('Calculating J-matrix....'))
        self._logger.info("%-15s : %s"%("Lambda",self._dict_params['lambda']))
        self._logger.info("%-15s : %s"%("Cut-off Rho",self._dict_params['rhocutoff']))
        self._logger.info("%-15s : %s"%("Scaling Coevolution",self._dict_params['scoe']))
        self._logger.info("%-15s : %s"%("Output File",self._dict_params['file_jmat']))
     
        count=1
        for i in (self._rho_data[0].split()[3:]):
            self._out_j+="%12s"%("J_"+"Vec_"+"%d"%(count))
            count+=1
        self._out_j+="\n"

        self._dict_rho_data={}
        for i in range(1,len(self._rho_data)):
            line_data=self._rho_data[i].strip().split()
            res1=int(line_data[0]); res2=int(line_data[1])

            i_key="%d-%d"%(res1,res2)
            i_scoe=self._scaled_matrix[res1-1,res2-1]
            i_coe=self._matrix_coevolution[res1-1,res2-1]
            self._out_j+="%12d%12d%12.2f%12.2f"%(res1,res2,i_coe,i_scoe)
            rho_values=[];
            for j in range(2,len(line_data)):
                j_rho=float(line_data[j])
                self._out_j+="%12.2f"%(self._calc_j_score(i_scoe,j_rho))
                rho_values.append(float(line_data[j]))
            self._out_j+="\n"
            self._dict_rho_data[i_key]=rho_values
        fileio.save_file(self._dict_params['file_jmat'],self._out_j)
    def _calc_j_score(self,scaled_coe,rho):
        jvalue=0.0; rho=np.abs(rho);
        if(scaled_coe>=1)and(rho>=self._dict_params["rhocutoff"]):
            jvalue  =   self._dict_params["lambda"]*(scaled_coe)+(1-self._dict_params["lambda"])*rho;
        return jvalue

    def _coe_matrix_stats(self):    
        self._logger.info('%-20s'%('COE PROPERTIES'))
        self._logger.info('%-20s : %d x %d'%('No. of residues',self._matrix_coevolution.shape[0],self._matrix_coevolution.shape[1]))
        self._avg_c_score         =   np.average(self._matrix_coevolution);
        self._max_c_score         =   np.max(self._matrix_coevolution);

        self._logger.info('%-20s : %.2f'%('Average',self._avg_c_score));
        self._logger.info('%-20s : %.2f'%('Maximum',self._max_c_score));

    def _calculate_scaled_matrix(self):
        self._scaled_matrix=self._matrix_coevolution;
        if(self._dict_params["scoe"]==True):
            self._scaled_matrix=self._matrix_coevolution/np.average(self._matrix_coevolution)


    def manager(self,dict_params):
        _start   =   timeit.default_timer();
        
        self._dict_params   =   dict_params;

        self._get_coev_matrix();
        self._get_rho_matrix();
        self._calculate_j_score()
        _stop    =   timeit.default_timer();
        self._logger.info('%-20s : %.2f (s)'%('FINISHED_IN',_stop-_start));
