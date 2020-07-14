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
import timeit,os,sys
import multiprocessing as mp
import dynoutil.dependencies as dependency
import dynoio.fileio as fileIO
import dynoutil.options as argParser

def main():
    global file_aln,file_ccm,pdbID,nthreads,dict_hhv
    args    =   argParser.opts_coevolution();
    pdbID   =   args.pdbid
    hhdb    =   args.database
    nthreads=   args.numthreads
main()
