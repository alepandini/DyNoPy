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
     but WITHOUT ANY WARRANTY without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with DyNoPy.  If not, see <http://www.gnu.org/licenses/>.

'''
import os
import re
import timeit
import logging
import numpy as np
import pandas as pd
import igraph as igh
import dynolib.networkslib as nwlib
import dynoio.fileio as fileio


class Networks(object):
    def __init__(self):
        self._initialize()

    def _initialize(self):

        self._logger = logging.getLogger('Dyno ReMa')
        '''
            booleans
        '''
        self._use_coevolution = False
        '''
            matrices
        '''
        self._matrix_coevolution = []
        self._matrix_rho = []

        '''
            default values
        '''

        self._lambda_j = 0.5
        self._avg_c_score = 0.0
        self._max_c_score = 1.0
        '''
            cut-offs
        '''
        self._scoe_cut = 1.0
        self._rho_cut = 0.0
        '''
            lists
        '''
        self._list_dist_data = []
        self._list_ie_data = []

    def _matrix_stats(self, _matrix):
        self._logger.info('%-20s' % ('MATRIX PROPERTIES'))
        self._logger.info(
            '%-20s : %d x %d' %
            ('No. of residues',
             _matrix.shape[0],
             _matrix.shape[1]))
        _avg = np.average(_matrix)
        _max = np.max(_matrix)
        self._logger.info('%-20s : %.2f' % ('Average', _avg))
        self._logger.info('%-20s : %.2f' % ('Maximum', _max))

    def _network_manager(self):
        self._matrix_j = []
        if (self._dict_params['file_coe'] is not None) and (
                self._dict_params['file_rho'] is not None):
            self._check_matrix_size()
            nvectors = self._matrix_rho.shape[1] - 2
            self._logger.info('CALCULATING NETWORK PROPERTIES')
            for i in range(nvectors):
                vec = i + 2
                label = "%s-v%d" % (self._dict_params['file_lab'], vec)
                self._logger.info('........ for vector : %5d' % (vec))
                self._matrix_j = nwlib.calculate_jmatrix(
                    self._matrix_rho,
                    self._matrix_coevolution,
                    rcutoff=0.5,
                    coecutoff=1.0,
                    vector_num=vec)
                nwlib.calculate_network_properties(self._matrix_j, label)

        if (self._dict_params['file_coe'] is None) and (
                self._dict_params['file_rho'] is not None):
            self._logger.info('CALCULATING NETWORK PROPERTIES')
            nwlib.calculate_network_properties(
                self._matrix_rho, self._dict_params['file_lab'])

    def _get_jmatrix_data(self):
        self._jmatrix_df = fileio.read_jmatrix(self._dict_params["file_jmat"])
        self._rename_residues()

    def _rename_residues(self):

        '''
            clean this in future
            is a quick fix to changing residue int to string
            pandas methods are not functional
        '''
        temp = []
        for x in self._jmatrix_df['Res_a']:
            temp.append(self._change_node_name(x))
        self._jmatrix_df['Res_a'] = temp

        temp = []
        for x in self._jmatrix_df['Res_b']:
            temp.append(self._change_node_name(x))
        self._jmatrix_df['Res_b'] = temp

    def _process_to_graph(self):
        self.full_graph = igh.Graph.DataFrame(
            self._jmatrix_df, directed=False, use_vids=False)
        self._num_nodes = len(self.full_graph.vs["name"])
        # assign default community id to 0
        self.full_graph.vs["community_id"] = self._num_nodes * [-1]

        # extract indices for each node
        self._dict_node_index = {}
        for index in self.full_graph.vs.indices:
            self._dict_node_index[self.full_graph.vs[index]["name"]] = index

    def _process_jmatrix(self):
        '''
            https://www.opentechguides.com/how-to/article/pandas/193/index-slice-subset.html
            https://igraph.org/python/tutorial/develop/tutorials/visualize_communities/visualize_communities.html
        '''

        self._vec_Q = []

        _vec_i = self._dict_params['vec_num']
        _choosen_vector = self._jmatrix_df.iloc[:, _vec_i]
        
        _vec_max = np.max(_choosen_vector)
        _vec_min = np.min(_choosen_vector)

        _vec_delta = (_vec_max - _vec_min) / self._dict_params['nsteps']
        _list_of_steps = np.arange(_vec_min, _vec_max + _vec_delta, _vec_delta)
        _cle = ""
        self._logger.info(
            "Starting calculation of Communities from the network")
        self._logger.info(
            "%-20s : %3d" %
            ("Analysis vector", self._dict_params['vec_num']))
        self._logger.info(
            "%-20s : %4.2f" %
            ("Modularity (Q) cut-off",
             self._dict_params['cut_mod']))
        self._logger.info("%-20s : %4.2f" % ("Vector Min.", _vec_min))
        self._logger.info("%-20s : %4.2f" % ("Vector Max.", _vec_max))
        self._logger.info(
            "%-20s : %4.2f" %
            ("No. of Steps for Search",
             self._dict_params['nsteps']))
        self.Q = 0 
        _list_Q = []
        #rename i
        for i in _list_of_steps:
            _filtered_df = self._jmatrix_df[_choosen_vector > i]
            _sliced_df = _filtered_df.iloc[:, [0, 1, _vec_i]]
            _mygraph = igh.Graph.DataFrame(
                _sliced_df, directed=False, vertices=None, use_vids=False)
            # igh.Graph.community_leading_eigenvector(_mygraph)
            _cle = _mygraph.community_leading_eigenvector()
            self.Q = _cle.modularity
            self._logger.info(
                "%-20s (Q) : %4.2f (%4.2f)" %
                ("Vector cut-off..", i, self.Q))

            if self.Q > self._dict_params['cut_mod']:
                self._logger.info("Desired Q, reached. Saving community graph")
                self._save_gml(_cle, _mygraph)
                self._save_residue_stats()
                break
        if self.Q < self._dict_params['cut_mod']:
            self._logger.info("Desired Q could not be reached. Either lower your Q cut-off or change your vector.")

    def _save_gml(self, community_graph, _sub_graph):
        self.get_evc_scores(_sub_graph)
        self._dict_resid_evc = {}

        num_communities = np.max(community_graph.membership)
        _community_id = 1
        _evc_out = "%15s%15s%15s\n" % ("residue_id", "community_id", "evc")

        for _subgraph in community_graph.subgraphs():
            list_of_node_names = _subgraph.vs["name"]
            for _node_name in list_of_node_names:
                _node_index = self._dict_node_index[_node_name]
                self.full_graph.vs[_node_index]["community_id"] = _community_id
                _evc_out += "%15d%15d%15.5f\n" % (self._get_int_residue_id(
                    _node_name), _community_id, self._evc_dictionary[_node_name])
                self._dict_resid_evc[_node_name] = _community_id

            _community_id += 1
        fileio.save_file(self._dict_params['file_evc'], _evc_out)
        igh.write(self.full_graph, self._dict_params['file_gml'])
    def _save_residue_stats(self):
        #get list of unique residues
        list_of_residue_ids=list(np.unique(self._jmatrix_df["Res_a"]))

        list_of_residue_ids.append(list(self._jmatrix_df["Res_b"])[-1])

        _j_col=list(self._jmatrix_df.columns)[4:]
        #_new_j_col=[]
        _out=""
        _out="Resid"+",CommunityID"+",EVC"+",C_CS"+",C_SCS"

        for i in _j_col:
            #_col_names.append("C_"+i)
            _out+=",%s"%(i)
        _out+="\n"

        for resid in list_of_residue_ids:
            int_resid=self._get_int_residue_id(resid)

            _df_a=self._jmatrix_df[self._jmatrix_df["Res_a"]==resid]
            _df_b=self._jmatrix_df[self._jmatrix_df["Res_b"]==resid]
            _community_id=-1

            # setting default evc value to -1
            _evc=-1

            if(resid in self._dict_resid_evc.keys()):
                _community_id=self._dict_resid_evc[resid]
                _evc=self._evc_dictionary[resid]
            _out+="%d,%d,%.5f"%(int_resid,_community_id,_evc)

            _res_stats=[]
            '''
                per residue calculate cummulative:
                                    coevolution score
                                    Scaled coevolution score
                                    J-value (J can be > 1)
                                    
            '''
            #print(int_resid,len(_df_a),len(_df_b),_df_a.empty,_df_b.empty)
            if(_df_a.empty):
                _res_stats=_df_b.sum(numeric_only=True)
            if(_df_b.empty):
                _res_stats=_df_a.sum(numeric_only=True)
            if(_df_a.empty==False and _df_b.empty==False):
                _df_a.add(_df_b)
                _res_stats=_df_a.sum(numeric_only=True)
            for _stat in _res_stats:
                _out+=",%.5f"%(_stat)

            _out+="\n"
        fileio.save_file(self._dict_params['file_stats'],_out)    
    def _change_node_name(self, node_name="NaN"):
        new_node_name = ""
        _node_name_int = np.rint(node_name).astype(int)
        node_name = str(_node_name_int)
        if (_node_name_int < 10):
            new_node_name = "R_000" + node_name
        if (_node_name_int >= 10 and _node_name_int <= 99):
            new_node_name = "R_00" + node_name
        if (_node_name_int >= 100 and _node_name_int <= 999):
            new_node_name = "R_0" + node_name
        if (_node_name_int >= 1000 and _node_name_int <= 9999):
            new_node_name = "R_" + node_name
        return new_node_name

    def _get_int_residue_id(self, node_name):
        return int(re.findall(r'\d+', node_name)[0])

    def get_evc_scores(self, graph):
        self._evc_list = igh.Graph.eigenvector_centrality(graph)
        self._evc_dictionary = {}
        count = 0
        for node_name in graph.vs['name']:
            self._evc_dictionary[node_name] = self._evc_list[count]
            count += 1

    def manager(self, dict_params):
        _start = timeit.default_timer()
        self._dict_params = dict_params
        self._get_jmatrix_data()
        self._process_to_graph()
        self._process_jmatrix()
        #self.save_residue_stats()
        _stop = timeit.default_timer()
        self._logger.info('%-20s : %.2f (s)' % ('FINISHED_IN', _stop - _start))
