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
import cProfile
import os,timeit,sys,glob
import numpy as np
from scipy import stats
from sklearn.metrics import mutual_info_score, normalized_mutual_info_score
import h5py,ray
#ray.init(num_cpus=4,num_gpus=0)

folder_main     =   "/home/scdantu/Projects/alessandro/prog_out/1lym/engine_data_analysis"
folder_coev     =   "%s/coevolution"%(folder_main);
folder_iexvg    =   "%s/iexvg"%(folder_main);

folder_main     =   os.getcwd();
folder_iexvg    =   "/mnt/pandinilab/scdantu/Projects/AMY7/myosin_sirahv20918_amber/iexvg/"
pdbID           =   "AMY7";

folder_main     =   os.getcwd();
folder_iexvg    =   "/home/scdantu/Projects/alessandro/blactamase/analysis/residue_identification/iexvg/"
pdbID           =   "2FFY";

#os.chdir(folder_main);

if(len(sys.argv)<5):
    print('jmatrix.py <co_ev-matrix> <dist> <replica number> <corr method>')
    print('Corr method: 0 PC; 1 SC; 2 NMI')
    exit()

file_mat    =   sys.argv[1];
file_dist   =   sys.argv[2];
nreplica    =   int(sys.argv[3])
corr_method =   int(sys.argv[4])
#sys_label   =   sys.argv[3];

serial=True;

res_first   =   1;    
res_last    =   358;

matrix_coevolution  =   [];  
matrix_rho          =   np.zeros((res_last,res_last));
lambda_j            =   0.5;
avg_c_score         =   0.0;    
max_c_score         =   1.0;

scoe_cut            =   1.0;
rho_cut             =   0.0;

list_dist_data      =   []; 
list_ie_data        =   [];

def get_coev_matrix():
    global avg_c_score,matrix_coevolution,max_c_score
    print('%-20s : %s'%('STORING','coevolution matrix...'))
    matrix_coevolution  =   np.loadtxt(file_mat)
    avg_c_score         =   np.average(matrix_coevolution);
    max_c_score         =   np.max(matrix_coevolution);

def calc_j_score(scaled_coe,rho):
    global scoe_cut,rho_cut,lambda_j
    jvalue=0.0; rho=np.abs(rho);
    
    if(scaled_coe>=scoe_cut)and(rho>=rho_cut):
        jvalue  =   lambda_j*(scaled_coe/max_c_score)+(1-lambda_j)*rho;
        
    return jvalue
def get_dist_data():
    global list_dist_data,file_dist
    list_dist_data   =   [];
    print('%-20s : %s'%('STORING', 'distance data...'))
    #fName="%s/D1-PRO-1-CG-2lcb-WT.xvg"%(folder_dist);
    fOpen=open(file_dist,'r');
    for line in fOpen:
        line        =   line.strip()
        split_line  =   line.split()
        list_dist_data.append(float(split_line[1]));
    list_dist_data=np.asarray(list_dist_data,dtype=np.float64)
    

def read_ie_h5(fName):
    global list_ie_data
    list_ie_data=[];
    hf = h5py.File(fName, 'r');
    n1 = np.array(hf.get('d1'))
    list_ie_data =  n1[:,0]+n1[:,1];
    #print(len(list_ie_data)
    #exit()
    #print("\n",np.sum(list_ie_data))
    
def calculate_scaled_c(c_value):
    global avg_c_score
    scaled_c=0;
    if(c_value>0):
        scaled_c=c_value/avg_c_score;
    return scaled_c

def rho_manager():
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

def main():
    start   =   timeit.default_timer();
    get_dist_data();
    get_coev_matrix();
    rho_manager();
    stop    =   timeit.default_timer();
    print('\n%-20s : %8.2f (s)'%('FINISHED_IN',stop-start));
main()
