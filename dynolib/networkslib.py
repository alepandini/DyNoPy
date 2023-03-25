import numpy as np
import logging
import collections
import dynoio.fileio as fileio
import dynoio.fileutils as futils
import dynoutil.rlib as rlib
import pandas as pd
logger = logging.getLogger('DyNo ReMa')


def calculate_jmatrix(m_rho, m_coe, rcutoff=0.5, coecutoff=1.0,
                      jcutoff=0.0, vector_num=1, jlambda=0.5):
    nrows = m_rho.shape[0]
    ncols = m_rho.shape[1]
    # adjusting for first two columns containing residue numbers A, B
    # a,b,rho_1,rho_2,rho_3,....,rho_n
    vector_num = 1 + vector_num
    '''
        bug proof?
    '''
    avg_coe = np.average(m_coe)
    n_res = m_coe.shape[0]
    matrix_j = np.zeros([n_res, n_res])
    for i in range(nrows):
        i_row = m_rho[i]
        res_a = int(i_row[0])
        res_b = int(i_row[1])
        i_rho = np.abs(i_row[vector_num])
        i_coe = (m_coe[res_a - 1][res_b - 1]) / avg_coe

        if(i_coe < coecutoff):
            i_coe = 0
        if(i_rho < rcutoff):
            i_rho = 0
        i_j = (jlambda * i_coe) + ((1 - jlambda) * i_rho)
        print(res_a, res_b, i_rho, i_coe, i_j)
        matrix_j[res_a - 1][res_b - 1] = i_j
    return matrix_j


def matrix_for_R(_matrix, _f_edges, _f_nodes):
    '''
        obsolete code
        not using R-igraph anymore
    '''
    res_i = 1
    res_j = 1
    out_edges = ""
    out_nodes = ""
    cut_off = 0.5
    _dict_res = {}
    for row in _matrix:
        res_j = 1
        for data in row:
            if(res_j > res_i) and (np.abs(data) >= cut_off):
                _dict_res[res_i] = 1
                _dict_res[res_j] = 1
                out_edges += "%5d%5d%12.5f\n" % (res_i, res_j, data)
                #out_nodes   +=  "%5d%5d\n"%(res_i,res_j)
            res_j += 1
        res_i += 1
    sorted_res_list = collections.OrderedDict(
        sorted(_dict_res.items(), key=lambda t: t[0]))
    for res, res_count in sorted_res_list.items():
        out_nodes += "%5d%5d\n" % (res, res_count)

    fileio.save_file(_f_edges, out_edges)
    fileio.save_file(_f_nodes, out_nodes)


def calculate_network_properties(_matrix, _label):
    '''
        obsolete code
        not using R-igraph anymore
    '''
    _f_edges = "%s.edg" % (_label)
    _f_nodes = "%s.nod" % (_label)
    _f_r = "%s.R" % (_label)
    matrix_for_R(_matrix, _f_edges, _f_nodes)
    rlib.save_communities_rscript(_label, _f_r)
    rlib.run_rscript(_f_r)


def read_jmatrix_to_dataframe(file_jmatrix):
    futils.check_file(file_jmatrix, cue_message="")
    data_frame = pd.read_fwf(file_jmatrix)
    return data_frame


'''
    df.loc[:,["Res_a","Res_b","Cab"]]
    df.iloc[:,[0,1,2,3]]

    def _save_community_graph(self, community_graph, main_graph):
        visual_style = {}
        visual_style["vertex_size"] = 25
        visual_style["vertex_label_size"] = 10
        visual_style["edge_width"] = 1
        visual_style["vertex_frame_width"] = 0
        visual_style["vertex_label_color"] = "black"
        num_communities = np.max(community_graph.membership)
        palette = igh.RainbowPalette(n=num_communities + 1)
        for i, community in enumerate(community_graph):
            s1 = community_graph.subgraph(i)
            vertex_labels = [int(x) for x in s1.vs["name"]]
            main_graph.vs[community]["color"] = i
        visual_style["vertex_label_size"] = 10
        visual_style["edge_width"] = 1
        visual_style["vertex_frame_width"] = 0
        visual_style["vertex_label_color"] = "black"
        num_communities = np.max(community_graph.membership)
        palette = igh.RainbowPalette(n=num_communities + 1)
        for i, community in enumerate(community_graph):
            s1 = community_graph.subgraph(i)
            vertex_labels = [int(x) for x in s1.vs["name"]]
            main_graph.vs[community]["color"] = i
            # mygraph.vs[community]["name"] = vertex_labels
            main_graph.vs[community]["vertex_label"] = vertex_labels
            main_graph.vs[community]["vertex_color"] = palette[i]
        visual_style["vertex_label"] = main_graph.vs["vertex_label"]
'''
