#!/usr/bin/python3.9
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
def opts_test():
    parser = argp.ArgumentParser()
    parser.add_argument("-i", "--pdbid",help="A pdbID to label all the files.",default="1TST")
    parser.add_argument("-f", "--firstres",default=1,help="First residue ID")
    parser.add_argument("-l", "--lastres",default=100,help="Last residue ID")
    parser.add_argument("-n", "--numdata",default=1000,help="number of data points for interaction energy/geometry vectors")
    parser.add_argument("-v", "--numvectors",default=1,help="number of geometry vectors")
    parser.add_argument("-t", "--numthreads",default=1,help="number of threads to use")

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
    parser = argp.ArgumentParser()
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
    parser.add_argument("-n", "--dt",
                        default=1,
                        help="Number of points"
                        );

    args = parser.parse_args()
    return args

def pwie():
    parser = argp.ArgumentParser()
    _usage_info="Calculate Pairwise Interaction Energies for selected residue range:\n"
    _usage_info+="dyno_pwie.py -x xtc/binpos -l label -t topology -f first_residue_id -l last_residue_id -p num_pairs -n num_threads\n"
    _usage_info+="e.g.: dyno_pwie.py -x KPC2.binpos -t KPC2_nowat.prmtop -f 1 -l 10 -p 10 -n 8 -o KPC2 \n";

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
def resrank():
    _usage_info="In a protein rank the residues by their importance:\n"
    _usage_info+="dyno_resrank.py -h\n"
    _usage_info+="e.g.: dyno_resrank.py -u P10191 \n";

    parser = argp.ArgumentParser(prog                   =   "resrank.py",
                                 usage="",description   =   _usage_info,
                                 formatter_class        =   argp.RawTextHelpFormatter
                                 );
    parser.add_argument("-u","--uniprot",
                        default=None,
                        help="uniprot ID"
                        );

    parser.add_argument("-p","--pdbid",
                        default=None,
                        help="pdb ID"
                        );
    args = parser.parse_args()
    return args
def resma():
    _usage_info="Extract the correlation matrix\n"
    _usage_info+="dyno_matrix.py -h\n"
    _usage_info+="e.g.: dyno_matrix.py -i iexvg -f 1 -l 10 -n 1\n";
    #parser = argp.ArgumentParser(prog                   =   "dyno_matrix.py",usage="",description   =   _usage_info,formatter_class        =   argp.RawTextHelpFormatter);
    parser = argp.ArgumentParser()
    parser.add_argument("-g","--gem",help="File with geometrical variables");
    parser.add_argument("-i","--fiex",help="Folder with interaction energy files iedata/");
    parser.add_argument("-f","--fst",help="First residue",default=1);
    parser.add_argument("-l","--lst",help="last residue",default=10);
    parser.add_argument("-n","--nrep",help="Number of replicas [NOT FUNCTIONAL]",default=1);
    parser.add_argument("-o","--label",help="label for output files");
    parser.add_argument("-t","--tmax",help="Number of threads. By default will use 80%% of the available threads");
    parser.add_argument("-m","--corr",help="Type of correlation method 0: Pearson; 1: Spearmen; 2: NMI",default=0);
    parser.add_argument("-v","--nvec",help="Number of Geom vectors to use 1...N default 1",default=1);
    args = parser.parse_args()
    return args

def jmatrix():
    _usage_info="Calculate the J-matrix:\n"
    _usage_info+="dyno_jmatrix.py -h\n"
    _usage_info+="e.g.: dyno_matrix.py -c coevolution.mat -r rho.mat -l 0.5 -o label\n";

    parser = argp.ArgumentParser(prog                   =   "dyno_jmatrix.py",
                                 usage="",description   =   _usage_info,
                                 formatter_class        =   argp.RawTextHelpFormatter
                                 );
    parser.add_argument("-c","--coe",
                        default=None,
                        help="coevolution matrix file"
                        );

    parser.add_argument("-r","--rho",
                        default=None,
                        help="File with RHO matrix"
                        );
    parser.add_argument("-o","--out",
                        default="J-Matrix.out",
                        help="Output File with J-matrix"
                        );

    parser.add_argument("-l","--dlambda",
                        default=0.5,
                        help="Set Lambda value, default is 0.5. [FULL SCAN NOT IMPLEMENTED]"
                        );
    parser.add_argument("-s",
                        "--scalescore",
                        default=True,
                        action=argp.BooleanOptionalAction,
                        help="Scale coevolution scores with average and set scores less then average to zero"
                        );
    parser.add_argument("--rhocutoff",
                        default=0.5,
                        help="Cut-off for absolute rho values "
                        );

    args = parser.parse_args()
    return args
