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
import os,logging
logger=logging.getLogger('DyNo IO  ')
def check_file(fName,cue_message=''):
    global logger
    if(os.path.isfile(fName)==False):
        logger.info('%-25s : %s. %s. EXITING','FILE_NOT_FOUND',fName,cue_message)
        exit()
    else:
        logger.debug('%-20s : %s','FILE_FOUND',fName)

#def check_exe(eName):
#    global logger
#    if(os.path.isfile(eName)==False):
#        logger.info('%-25s : %s. EXITING','EXE_NOT_FOUND',eName)
#        exit()
#    else:
#        logger.info('%-25s : %s','EXE_FOUND',eName)
def check_dir(dirN):
    if(os.path.isdir(dirN)==False):
        os.system('mkdir %s'%(dirN))
def checkfile_with_return(fileName):
    return os.path.isfile(fileName)

def checkfile_with_message(fileName,message):
    global logger
    if(os.path.isfile(fileName))==False:
        logger.info('%s does not exist. %s',fileName,message)
        exit();
def get_file_type(file_name):
    extension=os.path.splitext(file_name)
    return extension

def read_fasta(file_fasta):
    '''
        reads a fasta file
        returns a fasta_sequence dictionary 
        assumption the sequence numbering starts with 1
    '''
    global logger
    dict_fasta_sequence={}
    logger.debug("%-15s : %s"%("FASTA file",file_fasta))

    fileObject=open(file_fasta,'r')

    count=1
    for line in fileObject:
        if(line[0]!=">"):
            line=line.strip()
            for i in line:
                dict_fasta_sequence[count]=i;
                count+=1;
    return dict_fasta_sequence
