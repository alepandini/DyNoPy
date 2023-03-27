#!/usr/bin/python3
from gmlparser import *
f="../../GML-rmsd_j_wt.gml"
gml=GMLParser(f)
gml.create_edges_dict("Cab")
#a=gml.__get_res_number("R_0001")
#print(a)
