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
import dynoIO.fileIO as fileIO
import shutil,os,subprocess
def check_exe(exe_for_search):
    #shutil.which(exe_for_search)
    #rc=subprocess.getoutput("which %s"%(exe_for_search))
    return os.path.isfile()
def check_hhblits(dbname):
    hhlib   =   os.getenv("HHLIB")
    dict_hhv={};
    dict_hhv['hhlib']   =   hhlib;
    dict_hhv['hhblits'] =   hhlib+"/build/bin/hhblits";
    dict_hhv['hhfilter']=   hhlib+"/build/bin/hhfilter";
    if(hhlib==None):
        print('HHLIB variable not set. HHLIB should point to hh-suite directory')
        print('DyNoPy expects:\n\thhblits and hhfilter @ $HHLIB/build/bin/')
        print('\t uniprot files @ $HHLIB\database')
        print('Exiting.')
        exit()
    fileIO.check_exe(dict_hhv['hhblits'])
    fileIO.check_exe(dict_hhv['hhfilter'])
    fileIO.check_exe("ccmpred")
    '''
        add a more thorough check of the database folder
    '''
    file_hhdb="%s/database/%s_a3m.ffdata"%(hhlib,dbname)
    fileIO.check_file(file_hhdb,cue_message=dbname+" files not found. Download the files again.")
    dict_hhv['hhdb']="%s/database/%s"%(hhlib,dbname)
    return dict_hhv
