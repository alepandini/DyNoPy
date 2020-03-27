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
import shutil
import subprocess
def check_exe(exe_for_search):
    #shutil.which(exe_for_search)
    rc=subprocess.getoutput("which %s"%(exe_for_search))
    return rc
def check_hhblits(dbname):
    hhb_exe="hhblits";    
    hhf_exe="hhfilter";
    rc=check_exe(hhb_exe)
    if(check_exe(hhb_exe)==None):
        print('hhblits NOT_FOUND. INSTALL hh-suite. exiting...')
        exit()
    else:
        print('FOUND \t : hhblits @ %s'%(check_exe(hhb_exe)))
    if(check_exe(hhf_exe)==None):
        print('hhbfilter NOT_FOUND. INSTALL hh-suite. exiting...')
        exit()
    else:
        print('FOUND \t : hhfilter @ %s'%(check_exe(hhf_exe)))

    exit()
    '''
        WARNING: Set the name of uniclust database. Should be present in $hhpath
    '''
    hh_database="%s/%s/%s"%(hhpath,dbname,dbname)
    if(os.path.isdir(hhpath+'/'+dbname)==False):
        print('HHBLITS uniclust database not found in path: %s'%(hh_database))
        exit()
    else:
        print('HHBLITS DATABASE directory found...')

