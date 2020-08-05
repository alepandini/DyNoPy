#!/bin/python3
import os
import dynoio.fileio as fileio
import dynoio.fileutils as futils
def save_communities_rscript(label,fileName):
    rout="library(igraph)\n"
    rout+="source(\"/home/scdantu/Utilities/apps/mylibs/tempR/dpdnet.R\")\n"
    rout+="\n\n"
    rout+="#load the files\n"
    rout+="ssrn <- read.table('%s.nod',header=F)\n"%(label)
    rout+="ssre <- read.table('%s.edg',header=F)\n"%(label)
    rout+="ssrnet <- graph_from_data_frame(d=ssre,vertices=ssrn,directed=F)\n"
    rout+="ssrcle <- cluster_leading_eigen(ssrnet)\n"
    rout+="#save data\n"
    rout+="save_modularity(ssrcle)\n"
    rout+="save_comm_list(ssrnet,\"CL-%s.txt\")\n"%(label)
    rout+="save_graph_prop(ssrnet,\"GP-%s.txt\")\n"%(label)
    fileio.save_file(fileName,rout)
    
def run_rscript(fileName):
    futils.check_file(fileName)
    os.system("Rscript %s >/dev/null 2>&1"%(fileName))
    #os.system("Rscript %s "%(fileName))
