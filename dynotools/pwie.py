#!/usr/bin/python3
import numpy as np
from multiprocessing import Pool,Queue
from itertools import combinations 
import timeit,os,sys,logging
import dynoio.fileio as fileIO
import dynoio.fileutils as fUtils
import dynolib.pwielib as pwielib

class PWIE(object):

    def __init__(self):
        self._list_pairs    =   [];
        self._list_cpp      =   [];
        self._echo_out      =   "";
        self._dir_name      =   "iedata"
        self._logger        =   logging.getLogger('Dyno PIE')
    def _check_files(self):
        fUtils.check_file(self._file_top,'')
        fUtils.check_file(self._file_trj,'')
        fUtils.check_dir(self._dir_name);

    def _get_pairs(self):
        self._list_pairs    =   pwielib.gen_pair_list(self._resi_fst,self._resi_lst,gap=1)
        self._logger.info('%-25s : %8d'%('First Residue',self._resi_fst))
        self._logger.info('%-25s : %8d'%('Last  Residue',self._resi_lst))
        self._logger.info('%-25s : %8d'%('NPairs',len(self._list_pairs)))

 
    def _ech_head(self):	
        out_head="parm %s\n"%(self._file_top);
        out_head+="trajin %s 1 last\n"%(self._file_trj)
        return out_head

    def _ech_pair(self,fnum,lnum):
        out="pairwise :%d,%d out %s/IE-%d-%d-%s.dat cutevdw 0 cuteelec 0\ngo\n"%(fnum,lnum,self._dir_name,fnum,lnum,self._file_lab)
        return out

    def _gen_command(self):
        comm="cpptraj >/dev/null 2>&1 <<-EOF\n%s\nEOF"%(self._ech_out)
        return comm
    def _run_cpptraj(self):
        self._logger.info('%-25s : %8s'%('No. of pair sets',len(self._list_cpp)))
        pwielib.run_cpptraj_parallel(self._list_cpp,self._thrd_max,self._dir_name,self._file_lab);

    def _gen_cpp_command_list(self):
        _pair_count     =   0;	
        _pair_set_num   =   1;
        self._npairs    =   len(self._list_pairs);
        self._ech_out   =   self._ech_head();
        _temp_pairs     =   [];
        for count,pair in enumerate(self._list_pairs):
            self._ech_out+="%s"%(self._ech_pair(pair[0],pair[1]));
            _temp_pairs.append((pair[0],pair[1]))
            if(_pair_count==self._pair_max):
                self._list_cpp.append([_pair_set_num,self._gen_command(),_temp_pairs]);
                self._ech_out     =   self._ech_head();	
                _temp_pairs  =   [];    
                _pair_count  =   0;
                _pair_set_num+=1;

            if(_pair_count<self._pair_max)and(count+1==self._npairs):
                self._list_cpp.append([_pair_set_num,self._gen_command(),_temp_pairs])
            _pair_count+=1;	
 
    def analysis_manager(self,dict_params):
        self._pair_max  =   dict_params['pair_max'];
        self._resi_fst  =   dict_params['resi_fst'];
        self._resi_lst  =   dict_params['resi_lst'];
        self._file_lab  =   dict_params['file_lab'];
        self._thrd_max  =   dict_params['thrd_max'];
        self._file_trj  =   dict_params['file_trj'];
        self._file_top  =   dict_params['file_top'];
        self._check_files();
        self._get_pairs();
        self._gen_cpp_command_list();
        self._run_cpptraj();
        #check_files();
