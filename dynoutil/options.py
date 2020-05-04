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
import argparse as argp
def opts_coevolution():
    parser = argp.ArgumentParser()
    parser.add_argument("-i", "--pdbid",help="pdbID; NOTE: expects pdbID.fasta file in the directory")
    parser.add_argument("-d", "--database",default="UniRef30_2020_02",help="uniprot database ; Supports only UniRef30_2020_02")
    parser.add_argument("-n", "--numthreads",default=1,help="number of threads to use")
    args = parser.parse_args()
    if(args.pdbid==None):
        print('Provide pdbid. Exiting')
        exit()
        
    return args
def opts_aln_analysis():
    _usage_info="Sequence analysis tool helps you analyse the following: \n FASTA \n coevolution matrix."
    parser = argp.ArgumentParser(prog                   =   "dyno_aln_analysis.py",
                                 usage="",description   =   _usage_info,
                                 formatter_class        =   argp.RawTextHelpFormatter
                                 );
    parser.add_argument("-f","--fasta",
                        default=None,
                        required=True,
                        help="fasta file with the reference protein sequence"
                        );
    parser.add_argument("-m","--matrix",
                        default=None,
                        help="coevolution matrix file");
    parser.add_argument("-a", "--aln",
                        default=None,
                        help="number of threads to use");
    parser.add_argument("-l", "--label",
                        default="SeqAna",
                        help="the label of the output file");
    parser.add_argument("-t", "--type",
                        default =   0,
                        type    =   int,
                        choices =   [0,1,2,3,4],
                        help    =   " 0   : Print FASTA info\n 1   : Print alignment info\n 2   : Normalized coevolution matrix\n 3   : Per-residue coevolution cummulative score (PRCCS)\n 4   : log-odds matrix (Not Implemented)\n 5   : Not implemented\n99   : Not Implemented");
    args = parser.parse_args()
    return args
def conv_h5_to_ascii():
    _usage_info="Extract the data stored in h5 file with the pair key:\n dyno_h5_to_ascii.py -i Data.h5 -o 1-2-Log.txt -p 1-2"
    parser = argp.ArgumentParser(prog                   =   "dyno_h5_to_ascii.py",
                                 usage="",description   =   _usage_info,
                                 formatter_class        =   argp.RawTextHelpFormatter
                                 );
    parser.add_argument("-i","--h5",
                        default=None,
                        required=True,
                        help="fasta file with the reference protein sequence"
                        );
    parser.add_argument("-o","--out",
                        default="PAIR-ASCII.txt",
                        help="Output ASCII file"
                        );

    parser.add_argument("-p", "--pair",
                        default="1-2",
                        help="Pairs to extract. MULTIPLE PAIRS NOT SUPPORTED YET"
                        );
    args = parser.parse_args()
    return args

def pwie():
    _usage_info="Calculate Pairwise Interaction Energies for selected residue range:\n dyno_pwie.py -f xtc/binpos -l label -t topology  "
    parser = argp.ArgumentParser(prog                   =   "dyno_pwie.py",
                                 usage="",description   =   _usage_info,
                                 formatter_class        =   argp.RawTextHelpFormatter
                                 );
    parser.add_argument("-x","--trj",
                        default=None,
                        required=True,
                        help="trajectory file"
                        );

    parser.add_argument("-t","--top",
                        default=None,
                        required=True,
                        help="topology file"
                        );

    parser.add_argument("-f","--fst",
                        default=1,
                        type=int,
                        help="first residue"
                        );

    parser.add_argument("-l","--lst",
                        default=10,
                        type=int,
                        help="last residue"
                        );
    parser.add_argument("-p","--pmax",
                        default=10,
                        type=int,
                        help="maximum number of pairs"
                        );
    parser.add_argument("-n","--tmax",
                        default=10,
                        type=int,
                        help="number of threads to use"
                        );
    parser.add_argument("-o","--label",
                        default="PWIE",
                        help="Label to add to all the files"
                        );
    args = parser.parse_args()
    return args
