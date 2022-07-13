#!/usr/bin/python3.9
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
from dynotools.jmatrix import JMatrix

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
    args    =   argParser.jmatrix();
    dict_params['file_coe']     =   args.coe;
    dict_params['file_rho']     =   args.rho;
    dict_params['file_jmat']    =   args.out;
    dict_params['lambda']       =   float(args.dlambda)
    dict_params['scoe']         =   args.scalescore
    dict_params['rhocutoff']    =   float(args.rhocutoff)
    object_jma                 =   JMatrix();
    object_jma.manager(dict_params);

main()
