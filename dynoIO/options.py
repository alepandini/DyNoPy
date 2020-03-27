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

     DPD is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with DyNoPy.  If not, see <http://www.gnu.org/licenses/>.

'''
import argparse as argp
def opts_coevolution():
    parser = argp.ArgumentParser()
    parser.add_argument("-i", "--pdbid",help="pdbID; NOTE: expects pdbID.fasta file in the directory")
    parser.add_argument("-d", "--database",default="UniRef30_2020_02",help="uniprot database ; Supports only UniRef30_2020_02")
    args = parser.parse_args()
    if(args.pdbid==None):
        print('Provide pdbid. Exiting')
        exit()
        
    return args
