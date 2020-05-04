#!/usr/bin/python3
import os,timeit,sys
import numpy as np
import h5py,tqdm
from multiprocessing import Pool,Queue
from itertools import combinations 
dirN="";
labF="";
def gen_pair_list(fn,ln,gap=1):
    return list(combinations(np.arange(fn,ln+gap,gap), 2))
def run_command(cpp_command):
    #start	=	timeit.default_timer();
    os.system(cpp_command[1])
    gen_h5(cpp_command[2])
    #stop   =       timeit.default_timer();
    #print('%-10s%10s%8d : %8.2f (s)'%('FINISHED','pair_set',cpp_command[0],stop-start))

def gen_h5(pairs_list):
    global dirN,labF
    for res_a,res_b in pairs_list:
        datF="%s/IE-%d-%d-%s.dat"%(dirN,res_a,res_b,labF)
        h5pF="%s/IE-%d-%d-%s.h5"%(dirN,res_a,res_b,labF)
        if(os.path.isfile(datF)==True):
            npdata=returnData(datF);
            saveh5(h5pF,npdata);
            os.system('rm %s'%(datF))
def returnData(inF):
    n1=np.loadtxt(inF);
    ncols=n1.shape[1];
    if(ncols==3):
        n2 = np.empty([n1.shape[0],2],dtype=float)
        for index,i in enumerate(n1):
            n2[index][0]=i[1];
            n2[index][1]=i[2];
    if(ncols==2):
        n2=n1;
    return n2

def saveh5(outh5,npdata):
    hf = h5py.File(outh5, 'w')
    hf.create_dataset('d1',data=npdata,compression="gzip");
    hf.close()

def run_cpptraj_parallel_old(list_cpp,max_threads):
    _pbar   = tqdm.tqdm(total=100,desc="Dyno PIE - INFO - PIE @ ",position=0)
    pool=Pool(processes=max_threads)
    pool.map(run_command,[cpp_command for cpp_command in list_cpp]);
    pool.close()
    _pbar.update(100.0/self.NUM_OF_RESIDUES);

def run_cpptraj_parallel(list_cpp,max_threads,dirp,label):
    global dirN,labF
    labF=label; dirN=dirp;
    '''
        https://github.com/tqdm/tqdm/issues/484
    '''
    pool=Pool(processes=max_threads)
    for _ in tqdm.tqdm(pool.imap_unordered(run_command,list_cpp),desc="Dyno PIE - INFO - PIE Calculation           ",total=len(list_cpp)):
        pass
    pool.close()
    pool.join()
