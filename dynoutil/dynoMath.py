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
import numpy as np
import logging
logger=logging.getLogger('DyNo MAT')
def normalize_matrix(input_matrix,method=1):
    _factor=1;
    _shape_sum=np.sum(input_matrix.shape)

    if(_shape_sum==0):
        logger.info('%-25s : EXITING','ERROR_EMPTY_MATRIX')
    if(method==1):
        _factor=np.average(input_matrix)
    new_matrix=0;
    if(_factor>0):
        new_matrix=input_matrix/_factor;
    else:
        new_matrix=input_matrix;
    return new_matrix
def _per_res_sum(ires_data,cut_off=1.0,method=1):
    _sum=0; _num_gt=0;
    for i in ires_data:
        if(i>=cut_off):
            if(method==1):
                _sum+=i;
            if(method==2):
                _sum+=(i-1);
            _num_gt+=1;
    return _sum
def calc_prccs(matrix,method=1):
    _avg_cutoff=1.0;
    prccs="";
    N_rows  =   matrix.shape[0];
    for i in range(N_rows):
        prccs+="%12d%12.3f\n"%(i+1,_per_res_sum(matrix[i],cut_off=1.0,method=method))
    return prccs
