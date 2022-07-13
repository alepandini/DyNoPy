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
import dynoio.fileutils as fUtils
import dynoutil.options as argParser

perf_out="Performance stats...\n";

nthreads=1;

file_aln="";	path_hhdb=""; pdbID="";

dict_hhv={};

def run_hhblits():
    global perf_out
    start = timeit.default_timer()
    fasta="%s.fasta"%(pdbID);	hhr="%s-%d.hhr"%(pdbID,nthreads);	a3m="%s-%d.a3m"%(pdbID,nthreads)
    oa3m="%s.a3m"%(pdbID);
    ### check if fasta file is present
    fUtils.check_file(fasta);
    #dbname="uniclust30_2017_10"
    #hh_database="%s/%s/%s"%(hhpath,dbname,dbname)
    ### run hhblits
    print('Running hhblits on : %s'%(fasta))
    hh_com="%s -B 100000 -v 2 -n 4 -cpu %d -neffmax 20 -nodiff -maxfilt 100000 -d %s -i %s -o %s -oa3m %s"%(dict_hhv['hhblits'],nthreads,dict_hhv['hhdb'],fasta,hhr,a3m)
    os.system(hh_com)
    fUtils.check_file(a3m,'hhblits run might have failed...')

    ### run hhfilter
    print('Running hhfilter....')
    hh_filter="%s -id 90 -cov 75 -v 2 -i %s -o %s"%(dict_hhv['hhfilter'],a3m,oa3m)
    os.system(hh_filter)
    fUtils.check_file(oa3m,'hhfilter run might have failed...')
	
    ### check for unique sequences
    a3mtoaln="egrep -v \"^>\" %s | sed 's/[a-z]//g' | sort -u > %s"%(oa3m,file_aln)
    os.system(a3mtoaln)
    fUtils.check_file(file_aln)
    stop = timeit.default_timer()
    perf_out+="%-15s : %12.2f(s) ; (N_THREADS=%4d)\n"%("HHBLITS",stop-start,nthreads);

def run_ccmpred_gpu():
    global perf_out
    start = timeit.default_timer()
    print ("Running CCMPRED on %s .........."%(pdbID))
    ccmpred_com="ccmpred %s %s -n 75"%(file_aln,file_ccm)
    os.system(ccmpred_com)
    stop = timeit.default_timer()
    perf_out+="%-15s : %12.2f(s) ; (GPU)\n"%("CCMPRED",stop-start);
    perf_out+='%-15s : %12s\n'%('Co-evo Matrix',file_ccm)
def checkmaxthreads():
    global nthreads
    nthreads=int(nthreads)
    threads_max =   mp.cpu_count();
    if(nthreads>threads_max):
        nthreads=round(threads_max*0.8)
def main():
    global file_aln,file_ccm,pdbID,nthreads,dict_hhv
    args    =   argParser.opts_coevolution();
    pdbID   =   args.pdbid
    hhdb    =   args.database
    nthreads=   args.numthreads
    # check if  nthreads set is more than available
    checkmaxthreads()
    #check if hhblits and database are installed
    dict_hhv=dependency.check_hhblits(hhdb);
    # to control output file names
    file_aln="%s.aln"%(pdbID);
    file_ccm="%s.mat"%(pdbID);
    #run hhblits+tools
    run_hhblits()
    #run CCMpred
    run_ccmpred_gpu()
    print (perf_out)

main()
