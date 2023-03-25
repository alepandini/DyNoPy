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
     but WITHOUT ANY WARRANTY without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with DyNoPy.  If not, see <http://www.gnu.org/licenses/>.

'''
import os,gzip,numpy,logging
import dynoutil.hash_maps as hash_maps
import dynoio.fileutils as futils
class PDB(object):
    '''
    classdocs
    '''
    def __init__(self):
        
        '''
            Constructor
        '''
        self._initiate_logger()

        
        
        self.pdb_sequence=""
        
        self._file_name=""
        
        self.RESIDUES={}
        self.COORDINATES={}
        self.HETATM={}
 
    def _initiate_logger(self):
        self._logger=logging.getLogger('DyNo PDB')

    def get_pdb_sequence(self):
        '''
            generate sequence from the PDB file
        '''
        self._logger.debug('Extracting sequence from pdb...')
        i_keys=list(self.RESIDUES.values())
        seq=""
        for i in range(len(i_keys)):
            seq+="%s"%(i_keys[i])
        self.pdb_sequence=seq
        self._logger.info('%s',self.pdb_sequence)
        self._logger.info('%-20s : %s','Sequence length',len(self.pdb_sequence))
        if(len(self.pdb_sequence)==0):
            self._logger.info('PDB length is zero. File might be empty. Exiting')
            exit()
            
    def read_pdb(self,file_pdb=""):

        self._logger.debug('%-20s : %s','PDB file',file_pdb)
        futils.check_file(file_pdb)
        
        #pdb_lines=gzip.open(file_pdb,'rb')
        pdb_file_data=open(file_pdb,'r')

        for line in pdb_file_data:
            #decoded_line=line.decode('utf-8')
            line_split=line.split()
            if(line_split[0]=="ATOM"):
                self._process_pdb_line(line)
        
    def _process_pdb_line(self,line):
 
        residue     =       line[17:20]
        chain       =       line[21]
        resnum      =       int(line[22:26].strip())
        x           =       float(line[30:38].strip())
        y           =       float(line[38:46].strip())
        z           =       float(line[46:54].strip())
        
        res_code=hash_maps.res2aa(residue)
        '''
            ONLY STORES RESIDUES in 'A' format to the dict RESIDUES
        '''
        if(res_code=="UNK"):
            self.HETATM[resnum]=residue
        else:
            self.RESIDUES[resnum]=res_code
