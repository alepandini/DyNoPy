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
import dynoio.fileutils as fUtils
import dynolib.pwielib as pwielib
logger=logging.getLogger('DyNo IO ')

def read_h5_to_matrix(fName,keyw="d1"):
    fUtils.check_file(fName)
    hf      =   h5py.File(fName, 'r')
    mat     =   np.array(hf.get(keyw))
    return mat
def save_file(fName,data):
    f=open(fName,'w')
    f.write(data)
    f.close();

def read_file(filename):
    fUtils.check_file(filename)
    f=open(filename,'r');
    data="";
    for line in f:
        data+=line
    return data

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
    fUtils.check_file(in_h5);

    hf  =   h5py.File(in_h5,'r')
    if(pair in hf):
        data    =   hf.get(pair);
        data    =   np.array(data);
        x,y     =   data.shape;
        if(x==y):
            save_matrix(out_txt,data);
        if(x>y):
            out="";
            for count,d in enumerate(data):
                out+="%12d%12.5f%12.5f\n"%(count,float(d[0]),float(d[1]))
            save_file(out_txt,out)
def compress_h5(in_h5,pair,out_txt,dt=1):
    fUtils.check_file(in_h5);
    out=[]
    hf  =   h5py.File(in_h5,'r')
    #print(hf.keys())
    if(pair in hf):
        data    =   hf.get(pair);
        data    =   np.array(data);
        x,y     =   data.shape
        out=np.array(data[:,0]+data[:,1])
        out=out[0::dt]
        saveh5(out_txt,out)
def saveh5(outh5,npdata):
    hf = h5py.File(outh5, 'w')
    hf.create_dataset('d1',data=npdata,compression="gzip");
    hf.close()

def read_fasta(seq_file):
    global logger
    fUtils.check_file(seq_file)
    FASTA_SEQUENCE="";
    fileObject=open(seq_file,'r');
    for line in fileObject:
        if(line[0]!=">"):
            line=line.strip()
            for i in line:
                FASTA_SEQUENCE+="%s"%(i);
    return FASTA_SEQUENCE

def remove_grammar(word):
    if(word.find(",")>-1):
        word=word.split(",")[0]
    if(word.find(";")>-1):
        word=word.split(";")[0]
    return word

def save_record(fname,record):
    global logger
    logger.info('%-15s : %s','Saving',fname)
    with open(fname,'wb') as file:
        file.write(record.content)

def read_data_to_matrix(fName,limit=999999999):
    '''
        this method on a 600000,11 data file takes 2.2s while np.loadtxt takes 3.845s
    '''
    fUtils.check_file(fName)
    fOpen=open(fName,'r');
    fLine=fOpen.readlines();

    ncol=len(fLine[0].split())
    N=len(fLine);   count=0;
    if(N<=limit):
        Array_Data=np.empty([N,ncol],dtype=float)
    elif(N>limit):
        Array_Data=np.empty([limit,ncol],dtype=float)
    for i in range(N):
        if(i<limit):
            line=fLine[i].strip().split();
            for j in range(len(line)):
                Array_Data[i][j]=float(line[j]);
            count+=1;
    return Array_Data

