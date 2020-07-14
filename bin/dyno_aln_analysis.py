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
import timeit,os,sys,logging
import dynoio.fileio as fileIO
import dynoutil.options as argParser
from dynotools.sequence import Sequence

logger="";
dict_params={};
def initiate_logging():
    global logger
    '''
        input: coevolution matrix, sequence alignment, fasta sequence
            -- calculate frequencies of A and B
            -- generates pairs
            -- assign weights to each pair 
                @log frequencies
                @
    '''

    logging.basicConfig(format='%(name)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.INFO);
    logger=logging.getLogger('Dyno CANA')
def main():
    global dict_params,logger
    initiate_logging();
    args    =   argParser.opts_aln_analysis();
    dict_params['file_fasta']   =   args.fasta
    dict_params['file_aln']     =   args.aln
    dict_params['file_mat']     =   args.matrix
    dict_params['file_lab']     =   args.label
    dict_params['type_ana']     =   args.type
    object_seq_ana              =   Sequence();
    object_seq_ana.analysis_manager(dict_params);
main()
