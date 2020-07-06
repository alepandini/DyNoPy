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
from dynotools.pwie import PWIE

logger="";
dict_params={};
def initiate_logging():
    global logger
    '''
        input: trj,top,
            -- calculate PairWise Interaction Energies 
            -- generates pairs
            -- assign weights to each pair 
    '''
    logging.basicConfig(format='%(name)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.INFO);
    logger=logging.getLogger('Dyno PWIE')
def main():
    global dict_params,logger
    initiate_logging();
    args    =   argParser.pwie();
    dict_params['file_trj']     =   args.trj
    dict_params['file_top']     =   args.top
    dict_params['resi_fst']     =   int(args.fst)
    dict_params['resi_lst']     =   int(args.lst)
    dict_params['pair_max']     =   int(args.pmax)
    dict_params['file_lab']     =   args.label
    dict_params['thrd_max']     =   int(args.tmax)
    object_pwie                 =   PWIE();
    object_pwie.analysis_manager(dict_params);
main()
