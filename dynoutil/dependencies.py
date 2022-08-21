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
import logging,glob
logger=""
logger=logging.getLogger('Dyno Chks')

def check_exe(exe_for_search):
    global logger
    exloc=shutil.which(exe_for_search)
    if(exloc==None):
        logger.error("%-25s @ %s"%('EXE_NOT_FOUND',exe_for_search))
        exit()
    else:
        logger.info("%-25s @ %s"%("EXE_FOUND",exloc))
    #rc=subprocess.getoutput("which %s"%(exe_for_search))
def check_hhblits(dbname):
    hhlib   =   os.getenv("HHLIB")
    check_executables(hhlib)
    check_hhdb(hhlib,dbname)

def check_executables(hhlib):
    global logger
    logger.info("%-25s : %s"%("HHLIB",hhlib))
    logger.info('DyNoPy expects hhblits and hhfilter @ $HHLIB/build/bin/')
    logger.info('uniprot files @ $HHLIB\database')

    if(hhlib==None):
        logger.error('HHLIB variable not set. HHLIB should point to hh-suite directory')
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

def check_folder(fol_name):
    if(os.path.isdir(fol_name)==False):
        logger.error("%-25s : %s"%("FOLDER_NOT_FOUND",fol_name))
        logger.error("EXITING")
        exit()
    else:
        logger.info("%-25s : %s"%("FOLDER_FOUND",fol_name))

def check_hhdb(hhlib,dbname):
    global logger
    '''
        add a more thorough check of the database folder
    '''
    logger.info("%-25s : %s"%("Searching for database",dbname))
    fol_db      =   "%s/database"%(hhlib)
    fol_hhdb    =   "%s/%s"%(fol_db,dbname)
    file_hhdb   =   "%s/%s_a3m.ffdata"%(fol_hhdb,dbname)
    
    check_folder(fol_db)
    check_folder(fol_hhdb)
    
    list_of_database=glob.glob(fol_hhdb+"/")
    print(list_of_database)
    #if(os.path.isdir(fol_hhdb)==False):
    #    logger.error("%s does not contain the %s database in %s."%(hhlib,dbname,fol_hhdb))
    #    logger.error("Please check if the database has been downloaded/named properly. Exiting")
    #    exit()

    fUtils.check_file(file_hhdb,cue_message=dbname+" files not found. Download the files again.")
    dict_hhv['hhdb']="%s/database/%s/%s"%(hhlib,dbname,dbname)
    return dict_hhv
