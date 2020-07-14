#!/bin/python3
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

     ResPy is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with DyNoPy.  If not, see <http://www.gnu.org/licenses/>.

'''

import requests as req
import logging
import dynoio.fileio as fileIO

url_uniprot="https://www.uniprot.org/uniprot"

logger=logging.getLogger('Dyno IO ')

def remove_grammar(word):
    if(word.find(",")>-1):
        word=word.split(",")[0]
    if(word.find(";")>-1):
        word=word.split(";")[0]
    return word

def get_chains(word):
    first=word.split("=");
    n_chains=len(first[0].split("/"))
    return n_chains

def get_uniprot_record(uniprot_id,fmt="txt"):
    global logger
    logger.info('%-15s : %s',"Source URL",url_uniprot);
    logger.info('%-15s : %s',"Acquiring ID",uniprot_id);
    
    url_uni =   "%s/%s.%s"%(url_uniprot,uniprot_id,fmt)
    uni_rec =   get_record(url_uni)
    fname   =   "%s.%s"%(uniprot_id,fmt)

    fileIO.save_record(fname,uni_rec)
    
def get_record(myurl):
    uni_rec=req.get(myurl)
    return uni_rec

