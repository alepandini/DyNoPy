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
from optparse import OptionParser
print('IMPORTING_BAD_CODE options_optparser.py. EXITING')
exit()
def _list_archer():
    global logger
    out_usage ="dpd_list_simulations runs equilibration and generates files for archer.\n\n\n"
    out_usage+="DPD_Test.py -s list_of_mutations -l Run-1 -c cut_off -v 1 \n"
    
    parser = OptionParser(usage=out_usage)
    parser.add_option("-i", "--pdb",        dest="mpdb",   default="1lym",          help="pdbid (2abk)");
    parser.add_option("-s", "--list",       dest="mlist",   default="1lym.txt",     help="pdbid (2abk)");
    parser.add_option("-l", "--run",        dest="mrun",    default="Run-X",        help="Run label");
    parser.add_option("-n", "--nod",        dest="mnodes",    default="2",          help="num nodes");
    parser.add_option("-v", "--verbose",    dest="mverbose",default="-1",           help="-1\t:\tSILENT\n0\t:\tINFO\n1\t:\tDEBUG");
    (options, args) = parser.parse_args()
    if (int(options.mverbose)<-1)and (int(options.mverbose)>1):
        print('VERBOSE \t\t -1')
        options.mverbose=-1;
    
    return options

def opts_aln_analysis():
    usage_info="Ton of tools to analyze a sequence alignment and the coevolution matrix"
    parser = OptionParser(usage=out_usage)

    parser.add_option("-f","--fasta",
                      dest="fastaF",
                      default=None,
                      help="fasta file with the reference protein sequence"
                      );
    parser.add_option("-m","--matrix",
                      dest="matrixF",
                      default=None,
                      help="coevolution matrix file"
                      );
    parser.add_option("-a", "--aln",
                      dest="alnF"
                      default=None,
                      help="alignment file"
                      );
    (options, args) = parser.parse_args()
    return options
