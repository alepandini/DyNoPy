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
import logging
import dynoutil.options as argParser
from dynotools.resmatrix import ResMatrix

logger="";  dict_params={};
def initiate_logging():
    global logger
    '''
        input: geometrical variable*,folder with interaction energy data, coevolution matrix (optional), first res, last res, number of threads, correlation method, number of replicas
            -- calculate Rho matrix
            -- calculate J-Matrix if Coevolution matrix is provided
    '''
    logging.basicConfig(format='%(name)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.INFO);
    logger=logging.getLogger('Dyno ReMa')

def main():
    global dict_params,logger
    initiate_logging();
    args    =   argParser.resma();
    #dict_params['file_coe']     =   args.coe;
    dict_params['file_gem']     =   args.gem;
    dict_params['fold_iex']     =   args.fiex;
    dict_params['resi_fst']     =   int(args.fst)
    dict_params['resi_lst']     =   int(args.lst)
    dict_params['num_rep']      =   int(args.nrep)
    dict_params['file_lab']     =   args.label
    dict_params['num_thr']      =   int(args.tmax)
    dict_params['corr_met']     =   int(args.corr)
    dict_params['corr_vec']     =   int(args.nvec);
    object_rema                 =   ResMatrix();
    object_rema.manager(dict_params);
main()
