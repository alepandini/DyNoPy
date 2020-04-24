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
import numpy as np
import os,logging,h5py
import dynoIO.fileUtils as fUtils
logger=logging.getLogger('DyNo IO ')
def save_file(fName,data):
    f=open(fName,'w')
    f.write(data)
    f.close();

def read_aln(fName):
    global logger
    '''
        Make list of the aligned sequences
        each item is one sequence
    '''
    list_sequences=[];
    logger.info('Reading alignment...')
    fUtils.check_file(fName,' ');
    fileOpen    =   open(fName,'r');
    for line in fileOpen:
        a=(line.split('\n')[0])
        list_sequences.append(a)
    N=len(list_sequences);
    logger.info('%-25s : %s','Sequences (N)',N);
    if(N==0):
        logger.info('%-25s : EXITING.','ERROR_ALN_BLANK')
        exit()
    if(N==1):
        logger.info('%-25s : EXITING.','ERROR_ALN_SINGLE');
        exit()
    return list_sequences

def read_fasta(fastaF):
    global logger
    fUtils.check_file(fastaF,'')
    dict_fasta_sequence={};
    str_fasta_sequence="";
    logger.info('%-25s : %s','FASTA file',fastaF)
    fileObject=open(fastaF,'r');
    count=1;
    for line in fileObject:
        if(line[0]!=">"):
            line=line.strip()
            for i in line:
                dict_fasta_sequence[count]=i;
                str_fasta_sequence+=i
                count+=1;
    return str_fasta_sequence,dict_fasta_sequence

def read_matrix(matF):
    fUtils.check_file(matF,'')
    matrix=np.loadtxt(matF)
    return matrix
def save_matrix(fName,matrix):
    np.savetxt(fName,matrix,fmt='%10.5f')
def convert_h5_to_ascii(in_h5,pair,out_txt):
    hf  =   h5py.File(in_h5,'r')
    out_txt="%s-%s"%(pair,out_txt)
    if(pair in hf):
        data    =   hf.get(pair);
        data    =   np.array(data);
        x,y     =   data.shape;
        if(x==y):
            save_matrix(out_txt,data);
        if(x>y):
            out="";
            for d in data:
                out+="%6s %8.3f\n"%(d[0].decode('UTF-8'),float(d[1].decode('UTF-8')))
            save_file(out_txt,out)

