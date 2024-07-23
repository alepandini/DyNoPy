#!/bin/python3

import multiprocessing as mp
import threading
def get_max_threads():
    max_threads =   mp.cpu_count()
    return max_threads

def set_num_threads(nice_level=1):
    n_threads   =   1;

    max_threads =   get_max_threads();

    if(nice_level>0):
        '''
            10
            nl=1
                10*1/10     = 1
            nl=10
                10*10/10    = 10
        '''
        if(max_threads>1):
            n_threads   =   max_threads*(nice_level/10)
    return n_threads

