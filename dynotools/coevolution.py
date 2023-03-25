#!/usr/bin/env python
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
import timeit,os,sys,logging
import multiprocessing as mp
import dynoutil.dependencies as dependency
import dynoio.fileutils as futils
import dynoutil.options as argParser
from dynoio.pdb import PDB
from dynotools.sequence import Sequence

class Coevolution(object):
    def __init__(self):

        self._perf_out="Performance stats...\n"
        self._nthreads=1

        self._file_aln=""	
        self._path_hhdb="" 
        self._file_label=""

        self._dict_hhv={}
        self._logger=""
        self._initiate_logging()
        
        self._object_pdb=PDB()
        self._object_sequence=Sequence()

    def _initiate_logging(self):
        # initiate self._logger
        logging.basicConfig(format='%(name)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.INFO)
        self._logger=logging.getLogger('Dyno CoEv')
    
    def run_hhblits(self):
        # find similar sequences using the uniprot database
        # generate sequence alignment
        # filter sequence alignment using hhfilter
        start = timeit.default_timer()
        hhr="%s-%d.hhr"%(self._file_label,self._nthreads)
        a3m="%s-%d.a3m"%(self._file_label,self._nthreads)
        oa3m="%s.a3m"%(self._file_label)
        ### check if fasta file is present
        
        futils.check_file(self._file_fasta)
        #dbname="uniclust30_2017_10"
        #hh_database="%s/%s/%s"%(hhpath,dbname,dbname)
        ### run hhblits
        self._logger.info('Running hhblits on : %s'%(self._file_fasta))
        hh_com="%s -B 100000 -v 2 -n 4 -cpu %d -neffmax 20 -nodiff -maxfilt 100000 -d %s -i %s -o %s -oa3m %s"%(self._dict_hhv['hhblits'],self._nthreads,self._dict_hhv['hhdb'],self._file_fasta,hhr,a3m)
        os.system(hh_com)
        futils.check_file(a3m,'hhblits run might have failed...')

        ### run hhfilter
        self._logger.info('Running hhfilter....')
        hh_filter="%s -id 90 -cov 75 -v 2 -i %s -o %s"%(self._dict_hhv['hhfilter'],a3m,oa3m)
        os.system(hh_filter)
        futils.check_file(oa3m,'hhfilter run might have failed...')
            
        ### check for unique sequences
        a3mtoaln="egrep -v \"^>\" %s | sed 's/[a-z]//g' | sort -u > %s"%(oa3m,self._file_aln)
        os.system(a3mtoaln)
        futils.check_file(self._file_aln)
        stop = timeit.default_timer()
        self._logger.info("%-15s : %12.2f(s)  (N_THREADS=%4d)\n"%("HHBLITS",stop-start,self._nthreads))

    def run_ccmpred_gpu():
        # run ccmpred on GPU using the alignment file from HHBLITS
        # output is a coevolution matrix file (.mat)
        start = timeit.default_timer()
        self._logger.info("Starting CCMPRED")
        ccmpred_com="ccmpred %s %s -n 75"%(self._file_aln,self._file_ccm)
        os.system(ccmpred_com)
        stop = timeit.default_timer()
        self._logger.info("%-15s : %12.2f(s)  (GPU)\n"%("CCMPRED",stop-start))
        self._logger.info('%-15s : %12s\n'%('Coevolution Matrix file',self._file_ccm))

    def _checkmaxthreads(self):
        # check user input on max threads
        # if max threads>available threads, set them to 80% of the CPU power
        self._nthreads=int(self._nthreads)
        _threads_max =   mp.cpu_count()
        if(self._nthreads>_threads_max):
            self._nthreads=round(_threads_max*0.8)
            self._logger.info("Setting # Threads : %12d (80%% of max available)"%(self._nthreads))
        else:
            self._logger.info("Setting # Threads : %12d"%(self._nthreads))

    def run_coevolution(self,args):

        self._initiate_logging()

        self._logger.info("Starting Coevolution analysis...")
        self._logger.info("Will check if hhblits and ccmpred are available...")
        self._logger.info("Also will check for the chosen uniprot database")

        self._args = args
        self._file_fasta = self._args.fasta
        self._file_label = self._args.label
        self._hhdb = self._args.database
        self._nthreads = self._args.numthreads
        self._file_aln = "%s.aln"%(self._file_label)
        self._file_ccm = "%s.mat"%(self._file_label)

        # check if  nthreads set is more than available

        self._checkmaxthreads()
        
        #check if hhblits and database are installed
        self._dict_hhv = dependency.check_hhblits(self._hhdb)

        # to control output file names

        #run hhblits+tools
        self._run_hhblits()
        #run CCMpred
        self._run_ccmpred_gpu()

    def generate_residue_dictionary(self,args):
        self._file_pdb = args.pdb
        self._file_fasta = args.fasta

        self._object_sequence.get_fasta(self._file_fasta)
        
        self._object_pdb.read_pdb(self._file_pdb)
        self._object_pdb.get_pdb_sequence()
        self.pdb_sequence=self._object_pdb.pdb_sequence
        if(self._object_pdb.pdb_sequence!=self._object_sequence.FASTA_SEQUENCE):
            self._logger.warn("Sequence from the PDB file and the fasta file are not the same")
            self._logger.warn("Please trim the fasta sequence to match the pdb file")
            self._logger.warn("Future implementations of DyNoPy will be able to handle this mis-match. Sorry!!!")
            self._object_sequence.compare_fasta_and_pdb_sequence(self._object_sequence.FASTA_SEQUENCE,self.pdb_sequence)
        else:
            self._logger.info("Sequence from PDB and fasta file are a match. No issues. Phew!!!")
            # only for debugging purposes
            #self._object_sequence.compare_fasta_and_pdb_sequence(self._object_sequence.FASTA_SEQUENCE,self.pdb_sequence)






