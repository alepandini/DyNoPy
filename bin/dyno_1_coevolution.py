#!/usr/bin/env python
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

import dynoutil.options as argParser
from dynotools.coevolution import Coevolution
import logging
logger = ""


def initiate_logging():
    # initiate logger
    global logger
    logging.basicConfig(
        format='%(name)s - %(levelname)s - %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.INFO)
    logger = logging.getLogger('Dyno CoEv')

def main():
    initiate_logging()
    args = argParser.opts_coevolution()
    object_coevolution = Coevolution()
    #object_coevolution.generate_residue_dictionary(args)
    object_coevolution.run_coevolution(args)

main()
