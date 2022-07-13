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
import numpy as np
import multiprocessing as mp
import h5py,tqdm
from multiprocessing import Pool,Queue
import dynoutil.dependencies as dependency
import dynoio.fileio as fileIO
import dynoutil.options as argParser
import logging 
import dynolib.pwielib as pwielib

folname="iedata"
logger=""
numdata=1000; pdbID="1TST";
nthreads=4; res_first=1; res_last=100;
stddev=20;
def initiate_logging():
    global logger
    '''
        input: trj,top,
            -- calculate PairWise Interaction Energies 
            -- generates pairs
            -- assign weights to each pair 
    '''
    logging.basicConfig(format='%(name)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.INFO);
    logger=logging.getLogger('Dyno TEST')
def generate_matrix(n_rows,n_columns,file_name="Matrix.out"):
    global logger
    rand_matrix=np.random.rand(n_rows,n_columns)
    np.savetxt(file_name,rand_matrix,fmt="%12.2f")
    logger.info("%-15s : %s"%("Saved file",file_name))

def generate_pwie():
    global pdbID,nthreads,numdata
    global logger
    logger.info("Generating PWIE files...")

    if(os.path.isdir(folname)==False):
        os.mkdir(folname)
    list_pairs    =   pwielib.gen_pair_list(res_first,res_last,gap=1)
    logger.info("%-15s : %s"%("No. of pairs",len(list_pairs)))

    #gen_pwie_serial(list_pairs)
    gen_pwie_parallel(list_pairs)
def gen_pwie_pair(pair):
    global numdata,pdbID,logger
    out="";
    for i in range(numdata):
        a=np.random.normal(scale=stddev)
        b=np.random.normal(scale=stddev)
        out+="%12.2f %12.2f %12.2f\n"%(a,b,a+b)
    datF="%s/IE-%d-%d-%s.dat"%(folname,pair[0],pair[1],pdbID)
    outh5="%s/IE-%d-%d-%s.h5"%(folname,pair[0],pair[1],pdbID)

    f=open(datF,'w')
    f.write(out)

    npdata=pwielib.returnData(datF)
    pwielib.saveh5(outh5,npdata)
    os.system('rm %s'%(datF))

def gen_pwie_serial(list_pairs):
    global pdbID,nthreads,numdata
    for pair in list_pairs:
        gen_pwie_pair(pair)

def gen_pwie_parallel(list_pairs):
    global nthreads
    '''
        https://github.com/tqdm/tqdm/issues/484
    '''
    pool=Pool(processes=nthreads)
    for _ in tqdm.tqdm(pool.imap_unordered(gen_pwie_pair,list_pairs),desc="Dyno PIE - INFO - PIE Calculation           ",total=len(list_pairs)):
        pass
    pool.close()
    pool.join()
def generate_numvec():
    global pdbID,nthreads,numdata,numvectors
    out=""
    for i in range(numdata):
        out+="%12d"%(i+1)
        for j in range(numvectors):
            out+="%12.2f"%(np.random.normal(scale=stddev))
        out+="\n"
    f=open("Vectors-%s.dat"%(pdbID),'w')
    f.write(out)

def main():
    global pdbID,nthreads,numdata,numvectors
    global res_first,res_last
    initiate_logging()
    args    =   argParser.opts_test();
    pdbID   =   args.pdbid
    res_first   =  int(args.firstres)
    res_last    =   int(args.lastres)
    numdata     =   int(args.numdata)
    numvectors  =   int(args.numvectors)
    nthreads=   int(args.numthreads)
    generate_numvec()
    generate_pwie()
    generate_matrix(res_last,res_last,file_name="Coevolution-"+pdbID+".mat")
    generate_matrix(res_last,res_last,file_name="Dynamics-"+pdbID+".mat")

main()
