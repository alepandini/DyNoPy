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
import dynoio.fileutils as fUtils
import shutil,os,subprocess
import logging
logger=""
logger=logging.getLogger('Dyno Chks')

def check_exe(exe_for_search):
    global logger
    exloc=shutil.which(exe_for_search)
    if(exloc==None):
        logger.error('EXE_NOT_FOUND %s'%(exe_for_search))
        exit()
    else:
        logger.info('EXE_FOUND     @ %s'%(exloc))
    #rc=subprocess.getoutput("which %s"%(exe_for_search))
def check_hhblits(dbname):
    global logger
    hhlib   =   os.getenv("HHLIB")
    if(hhlib==None):
        logger.error("HHLIB : %s"%(hhlib))
        logger.error('HHLIB variable not set. HHLIB should point to hh-suite directory')
        logger.error('DyNoPy expects:\n\thhblits and hhfilter @ $HHLIB/build/bin/')
        logger.error('\t uniprot files @ $HHLIB\database')
        logger.error('Exiting.')
        exit()

    dict_hhv={};
    dict_hhv['hhlib']       =   hhlib;
    #dict_hhv['hhblits']     =   hhlib+"/build/bin/hhblits";
    #dict_hhv['hhfilter']    =   hhlib+"/build/bin/hhfilter";
    dict_hhv['hhblits']     =   shutil.which("hhblits");
    dict_hhv['hhfilter']    =   shutil.which("hhfilter")

    check_exe(dict_hhv['hhblits'])
    check_exe(dict_hhv['hhfilter'])
    check_exe("ccmpred")
    '''
        add a more thorough check of the database folder
    '''
    file_hhdb="%s/database/%s/%s_a3m.ffdata"%(hhlib,dbname,dbname)
    fol_hhdb="%s/database/%s"%(hhlib,dbname)
    if(os.path.isdir(fol_hhdb)==False):
        logger.error("%s does not contain the %s database in %s."%(hhlib,dbname,fol_hhdb))
        logger.error("Please check if the database has been downloaded/named properly. Exiting")
        exit()
    fUtils.check_file(file_hhdb,cue_message=dbname+" files not found. Download the files again.")
    dict_hhv['hhdb']="%s/database/%s/%s"%(hhlib,dbname,dbname)
    return dict_hhv
