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
import timeit,os,sys
import dynoutil.dependencies as dependency
import dynoIO.fileIO as fileIO
import dynoIO.options as argParser

perf_out="Performance stats...\n";

file_aln="";	path_hhdb=""; pdbID="";

def run_hhblits(num_threads=4):
    global perf_out
    start = timeit.default_timer()
    fasta="%s.fasta"%(pdbID);	hhr="%s-%d.hhr"%(pdbID,num_threads);	a3m="%s-%d.a3m"%(pdbID,num_threads)
    oa3m="%s.a3m"%(pdbID);
    ### check if fasta file is present
    fileIO.check_file(fasta);
    
    #dbname="uniclust30_2017_10"
    #hh_database="%s/%s/%s"%(hhpath,dbname,dbname)
    ### run hhblits
    print('Running hhblits on : %s'%(fasta))
    hh_com="hhblits -B 100000 -v 2 -n 4 -cpu %d -neffmax 20 -nodiff -maxfilt 100000 -d %s -i %s -o %s -oa3m %s"%(num_threads,path_hhdb,fasta,hhr,a3m)
    os.system(hh_com)
    fileIO.check_file(a3m,'hhblits run might have failed...')

    ### run hhfilter
    print('Running hhfilter....')
    hh_filter="hhfilter -id 90 -cov 75 -v 2 -i %s -o %s"%(a3m,oa3m)
    os.system(hh_filter)
    fileIO.check_file(oa3m,'hhfilter run might have failed...')
	
    ### check for unique sequences
    a3mtoaln="egrep -v \"^>\" %s | sed 's/[a-z]//g' | sort -u > %s"%(oa3m,file_aln)
    os.system(a3mtoaln)
    fileIO.check_file(file_aln)
    stop = timeit.default_timer()
    perf_out+="%-15s : %12.2f(s) ; (N_THREADS=%4d)\n"%("HHBLITS",stop-start,num_threads);

def run_ccmpred_gpu(ccm_file,num_threads=0):
    global perf_out
    start = timeit.default_timer()
    print ("Running CCMPRED on %s .........."%(pdbID))

    ccmpred_com="ccmpred %s %s -n 75"%(file_aln,file_ccm)
    os.system(ccmpred_com)
    stop = timeit.default_timer()
    perf_out+="%-15s : %12.2f(s) ; (GPU)\n"%("CCMPRED",stop-start);
    perf_out+='%-15s : %12s\n'%('Co-evo Matrix',ccm_file)

def main():
    global file_aln,file_ccm,pdbID
    args    =   argParser.opts_coevolution();
    pdbID   =   args.pdbid
    hhdb    =   args.database
    
    #check if hhblits and database are installed
    dependency.check_hhblits(hhdb);

    # to control output file names
    file_aln="%s.aln"%(pdbID);
    file_ccm="%s.mat"%(pdbID);
    
    #run hhblits+tools
    run_hhblits(num_threads=24)
    #run CCMpred
    run_ccmpred_gpu(num_threads=1)
    
    print (perf_out)

main()
