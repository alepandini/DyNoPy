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
from numba import njit
from numba import jit
import os,timeit,sys,glob
import numpy as np
from scipy import stats
from sklearn.metrics import mutual_info_score, normalized_mutual_info_score
import h5py,ray,logging
class ResMatrix(object):
    def __init__(self):
        self._initialize();
        
    def _initialize(self):
        
        self._serial=True;
        
        self._logger    =   logging.getLogger('Dyno ReMa')

        self._matrix_coevolution  =   [];  
        self._matrix_geometric    =   [];
        self._matrix_rho          =   np.zeros((res_last,res_last));
        self._lambda_j            =   0.5;
        self._avg_c_score         =   0.0;    
        self._max_c_score         =   1.0;

        self._scoe_cut            =   1.0;
        self._rho_cut             =   0.0;

        self._list_dist_data      =   []; 
        
        self._list_ie_data        =   [];

    def _get_coev_matrix(self):
        self._logger.info('%-20s : %s'%('STORING','coevolution matrix...'))
        

        self._matrix_coevolution  =   np.loadtxt(self._dict_params['file_coe'])
        self._avg_c_score         =   np.average(matrix_coevolution);
        self._max_c_score         =   np.max(matrix_coevolution);

    def _calc_j_score(self,scaled_coe,rho):
        global scoe_cut,rho_cut,lambda_j
        jvalue=0.0; rho=np.abs(rho);
        
        if(scaled_coe>=scoe_cut)and(rho>=rho_cut):
            jvalue  =   lambda_j*(scaled_coe/max_c_score)+(1-lambda_j)*rho;
        return jvalue

    def _get_geometric_data(self):
        self._logger.info('%-20s : %s'%('STORING', 'distance data...'))
        exit()
        

    def _read_ie_h5(fName):
        global list_ie_data
        list_ie_data=[];
        hf = h5py.File(fName, 'r');
        n1 = np.array(hf.get('d1'))
        list_ie_data =  n1[:,0]+n1[:,1];
        #print(len(list_ie_data)
        #exit()
        #print("\n",np.sum(list_ie_data))
        
    def _calculate_scaled_c(c_value):
        global avg_c_score
        scaled_c=0;
        if(c_value>0):
            scaled_c=c_value/avg_c_score;
        return scaled_c

    def _rho_manager(self):
        global matrix_rho,matrix_coevolution,avg_c_score
        
        #out_data="%5s%5s%12s%12s%12s%12s\n"%('Res_a','Res_b','Cab','Sab','Rab','j');
        out_data="%8s%8s%12s%12s%12s\n"%('Res_a','Res_b','Cab','Sab','Rab');
        print('%-20s : %s'%('CALCULATING','rho...'))
        start=timeit.default_timer()
        for i in range(res_first,res_last+1,1):
            
            for j in range(res_first,res_last+1,1):
                if(j>i):
                    flabel                  =   '%d-%d'%(i,j);
                    
                    fname="%s/IE-%d-%d-PRO-%d-CG-2lcb-WT.h5"%(folder_iexvg,i,j,nreplica)
                    fname="%s/IE-%d-%d-%d-%s.h5"%(folder_iexvg,i,j,nreplica,pdbID)
                    fname="%s/IE-%d-%d-%s.h5"%(folder_iexvg,i,j,pdbID)
                    
                    #fName                   =   "%s/IE-%d-%d-PRO-1-CG-2lcb-WT.dat"%(folder_iexvg,i,j)
                    print_t='%-20s : %5d %5d'%('PROCESSING_PAIR',i,j);
                    sys.stdout.write("\r{0}".format(print_t));
                    sys.stdout.flush()

                    rho_value               =   calculate_r(fname)

                    c_value                 =   matrix_coevolution[i-1][j-1]
                    matrix_rho[i-1,j-1]     =   rho_value;
                    matrix_rho[j-1,i-1]     =   rho_value;
                    
                    scaled_c                =   calculate_scaled_c(c_value);
                    #out_data                +=  "%5d%5d%12.5f%12.5f%12.5f%12.5f\n"%(i,j,c_value,scaled_c,rho_value,calc_j_score(scaled_c, rho_value));
                    out_data                +=  "%8d%8d%12.5f%12.5f%12.5f\n"%(i,j,c_value,scaled_c,rho_value);
                    #print(out_data)
                    #exit()
            if(i==9999):
                stop=timeit.default_timer()
                print('\nFinished in : %12.8f (s)'%(stop-start))
                exit()
            #break
        rhotype="PC"
        if(corr_method==0):
            rhotype="PC"
        elif(corr_method==1):
            rhotype="SC"
        if(corr_method==2):
            rhotype="NMI"
        flabel="%s-%d-%s.mat"%(rhotype,nreplica,pdbID)
        fwrite=open('J-'+flabel,'w');    
        fwrite.write(out_data);
        fwrite.close()
        #np.savetxt(flabel,matrix_rho,fmt="%12.5f")
        #np.savetxt(rhotype+'-CoE-2lcb.mat',matrix_coevolution/avg_c_score,fmt="%12.5f")


    def calculate_r(fName):
        rho_value       =   0;

        if(os.path.isfile(fName)==True):
            read_ie_h5(fName)
        if(len(list_ie_data)==len(list_dist_data)):
            rho_value       =   correlation_analysis()
        else:
            rho_value=0;
        return rho_value

    def correlation_analysis():
        rho_value=0.0
        if(serial==True):
            #print('serial')
            rho_value=serial_calculate_correlation()
        else:
            #print('parallel')
            rho_value=ray.get(parallel_correlation_analysis.remote())
        #print(" ",rho_value)
        return rho_value
    def serial_calculate_correlation():
        global corr_method
        global list_dist_data,list_ie_data
        rho = 0.0;
        if(corr_method==0):
            rho=np.corrcoef(list_dist_data,list_ie_data)[0,1];
        elif(corr_method==1):
            rho,p=stats.spearmanr(list_dist_data,list_ie_data);
        elif(corr_method==2):
            score = normalized_mutual_info_score(list_dist_data,list_ie_data,average_method="arithmetic")
            rho = np.nan_to_num(score)
        return rho
    #@ray.remote
    def parallel_correlation_analysis():
        global corr_method
        global list_dist_data,list_ie_data
        rho = 0.0;
        if(corr_method==0):
            rho=np.corrcoef(list_dist_data,list_ie_data)[0,1];        
        elif(corr_method==1):
            rho,p=stats.spearmanr(list_dist_data,list_ie_data);
        elif(corr_method==2):
            score = normalized_mutual_info_score(list_dist_data,list_ie_data,average_method="arithmetic")   
            rho = np.nan_to_num(score)
        return rho

    def manager(self,dict_params):
        _start   =   timeit.default_timer();
        
        self._dict_params   =   dict_params;

        self._get_geometric_data();
        self._get_coev_matrix();
        self._rho_manager();
        _stop    =   timeit.default_timer();
        self._logger.info('\n%-20s : %8.2f (s)'%('FINISHED_IN',stop-start));
