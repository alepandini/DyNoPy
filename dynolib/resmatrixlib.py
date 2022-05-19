import numpy as np
import h5py,tqdm,logging,timeit
from multiprocessing import Pool,Queue
from itertools import combinations 
from scipy import stats
from sklearn.metrics import mutual_info_score, normalized_mutual_info_score

import dynoio.fileio as fileio
import dynoio.fileutils as futils
dirN="";
labF="";
corr_params =   {}; list_ie_data=[];
logger=logging.getLogger('DyNo ReMa')
def calculate_correlation(pair_tuple):
    '''
        corr_params [1] : list of pairs
        corr_params [2] : matrix of geometric data
        corr_params [4] : file label
        corr_params [3] : number of replicas
        corr_params [5] : number of threads
        corr_params [6] : correlation method 1,2
        corr_params [7] : folder containing the data
        corr_params [8] : number of vectors to calculate correlation with
    '''
    pair_a  =   pair_tuple[0];
    pair_b  =   pair_tuple[1];
    ie_data =   get_ie_data(pair_a,pair_b);
    ie_len  =   ie_data.shape[0];
    di_len  =   corr_params[2].shape[0];
    corr_m  =   corr_params[6];

    vec_coef  =   [];
    max_vecs    =   1;

    if(ie_len==di_len):
        for i in range(1,corr_params[8]+1):
            #cor.append(np.corrcoef(d[:,i][:200000],n1[:,1])[0,1])
            a   =   corr_params[2][:,i]
            vec_coef.append(get_correlation(a,ie_data));
    if(ie_len!=di_len):
        logger.info('%-20s : %4d %4d %8d %8d'%('LENGTH_MISMATCH_EXITING',pair_a,pair_b,ie_len,di_len))
        exit()
    return (pair_a,pair_b,vec_coef)


def get_correlation(a,b):
    r_coef  =   0.0;
    corr_m  =   corr_params[6];
    if(corr_m==0):
        r_coef      =   np.corrcoef(a,b)[0,1];
    elif(corr_m==1):
        r_coef,p    =   stats.spearmanr(a,b);
    elif(corr_m==2):
        score       =   normalized_mutual_info_score(a,b,average_method="arithmetic")
        r_coef      =   np.nan_to_num(score)
    return r_coef

def get_ie_data(pair_a,pair_b):
    list_ie_data    =   [];
    '''
        IE-1-32-R1-CG-2gdn.h5
    '''
    for i in range(1,corr_params[3]+1):
        fname_ie="%s/IE-%d-%d-%s.h5"%(corr_params[7],pair_a,pair_b,corr_params[4])
        if(futils.checkfile_with_return(fname_ie)==True):
            ie_matrix   =   fileio.read_h5_to_matrix(fname_ie);

            # calculating interaction energy by adding VdW+Elec
            ie          =   np.add(ie_matrix[:,0],ie_matrix[:,1])
            for j in range(len(ie)):
                list_ie_data.append(ie[j]);
    if(len(list_ie_data)>1):
        list_ie_data=np.asarray(list_ie_data)
    if(len(list_ie_data)<=1):
        list_ie_data=np.zeros((1,))
    return list_ie_data

def correlation_calculator(c_params):
    global corr_params,logger
    serial=False
    corr_params =   c_params
    '''
        https://github.com/tqdm/tqdm/issues/484
    '''
    start_t =   timeit.default_timer();

    list_results=[];
    if(serial==False):
        pool=Pool(processes=corr_params[5])
        #pool=Pool(processes=1);
        csize   = round(len(corr_params[1])*0.2)
        if(csize > 40):
            csize=40;
        logger.info('%-20s : %s'%('Chunk size',csize))
        for results in tqdm.tqdm(pool.imap(calculate_correlation,corr_params[1],chunksize=csize),desc="Dyno ReMa - INFO - R Calculation        ",total=len(corr_params[1])):
            list_results.append(results)
    else:
        for l in corr_params[1]:
            list_results.append(calculate_correlation(l))
    end_t   =   timeit.default_timer();
    logger.info('%-20s : %.2f (s)'%('Took',end_t-start_t))
    return list_results
