import numpy as np
import logging,collections
import dynoio.fileio as fileio
import dynoio.fileutils as futils
import dynoutil.rlib   as rlib
logger=logging.getLogger('DyNo ReMa')

def matrix_for_R(_matrix,_f_edges,_f_nodes):
    res_i   =   1;      res_j   =   1;
    out_edges   =   ""; out_nodes   =   "";
    cut_off     =   0.5;
    _dict_res   =   {};
    for row in _matrix:
        res_j=1;
        for data in row:
            if(res_j > res_i)and(np.abs(data)>=cut_off):
                _dict_res[res_i]=1;
                _dict_res[res_j]=1;
                out_edges   +=  "%5d%5d%12.5f\n"%(res_i,res_j,data)
                #out_nodes   +=  "%5d%5d\n"%(res_i,res_j)
            res_j+=1;
        res_i+=1;
    
    sorted_res_list =   collections.OrderedDict(sorted(_dict_res.items(),key=lambda t: t[0]))
    for res,res_count in sorted_res_list.items():
        out_nodes   +=  "%5d%5d\n"%(res,res_count)
    
    fileio.save_file(_f_edges,out_edges)
    fileio.save_file(_f_nodes,out_nodes)

def calculate_network_properties(_matrix,_label):
    _f_edges    =   "%s.edg"%(_label)
    _f_nodes    =   "%s.nod"%(_label)
    _f_r        =   "%s.R"%(_label)
    matrix_for_R(_matrix,_f_edges,_f_nodes)
    rlib.save_communities_rscript(_label,_f_r)
    rlib.run_rscript(_f_r)
