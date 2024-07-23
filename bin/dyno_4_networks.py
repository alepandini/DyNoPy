#!/usr/bin/env python

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
from dynotools.networks import Networks

logger = ""
dict_params = {}


def initiate_logging():
    global logger
    '''
        input: geometrical variable*,folder with interaction energy data, coevolution matrix (optional), first res, last res, number of threads, correlation method, number of replicas
            -- calculate Rho matrix
            -- calculate J-Matrix if Coevolution matrix is provided
    '''
    logging.basicConfig(
        format='%(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.INFO)
    logger = logging.getLogger('Dyno Netw')


def main():
    global dict_params, logger
    initiate_logging()
    args = argParser.networks()
    dict_params['file_jmat'] = args.fjmat
    dict_params['cut_mod'] = float(args.mod)
    dict_params['file_gml'] = args.fgml
    dict_params['nsteps'] = int(args.nsteps)
    dict_params['vec_num'] = int(args.nvec)
    dict_params['file_evc'] = args.fevc
    dict_params['file_stats'] = args.fstats
    dict_params['scan_q'] = False
    dict_params['vec_cutoff'] =  0.0
    dict_params['npairs'] = 10
    object_netw = Networks()
    object_netw.manager(dict_params)

main()
