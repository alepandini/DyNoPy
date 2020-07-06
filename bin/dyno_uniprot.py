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
from dynotools.uniprot import Uniprot

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
    logger=logging.getLogger('ResRank')
def main():
    global dict_params,logger
    initiate_logging();
    args    =   argParser.resrank();
    dict_params['uniprotid']    =   args.uniprot
    object_uniprot              =   Uniprot();
    object_uniprot.manager(dict_params);
main()
