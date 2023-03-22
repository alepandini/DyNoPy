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
import os,timeit,logging
import numpy as np
import pandas as pd
import igraph as igh
import dynolib.networkslib as nwlib 
import dynoio.fileio as fileio

class Networks(object):
    def __init__(self):
        self._initialize();
        
    def _initialize(self):
        
        self._logger    =   logging.getLogger('Dyno ReMa')
        '''
            booleans
        '''
        self._use_coevolution       =   False
        '''
            matrices
        '''
        self._matrix_coevolution    =   [];  
        self._matrix_rho            =   [];

        '''
            default values
        '''

        self._lambda_j            =   0.5;
        self._avg_c_score         =   0.0;    
        self._max_c_score         =   1.0;
        '''
            cut-offs
        '''
        self._scoe_cut            =   1.0;
        self._rho_cut             =   0.0;
        '''
            lists
        '''
        self._list_dist_data      =   []; 
        self._list_ie_data        =   [];
        
    def _matrix_stats(self,_matrix):    
        self._logger.info('%-20s'%('MATRIX PROPERTIES'))
        self._logger.info('%-20s : %d x %d'%('No. of residues',_matrix.shape[0],_matrix.shape[1]))
        _avg    =   np.average(_matrix);
        _max    =   np.max(_matrix);
        self._logger.info('%-20s : %.2f'%('Average',_avg));
        self._logger.info('%-20s : %.2f'%('Maximum',_max));

    def _network_manager(self):
        self._matrix_j  =   [];
        if(self._dict_params['file_coe']!=None)and(self._dict_params['file_rho']!=None):
            self._check_matrix_size();
            nvectors    =   self._matrix_rho.shape[1]-2;
            self._logger.info('CALCULATING NETWORK PROPERTIES')
            for i in range(nvectors):
                vec     =   i+2; 
                label   =   "%s-v%d"%(self._dict_params['file_lab'],vec);
                self._logger.info('........ for vector : %5d'%(vec))
                self._matrix_j  =   nwlib.calculate_jmatrix(self._matrix_rho,self._matrix_coevolution,rcutoff=0.5,coecutoff=1.0,vector_num=vec);
                nwlib.calculate_network_properties(self._matrix_j,label);
        
        if(self._dict_params['file_coe']==None)and(self._dict_params['file_rho']!=None):
            self._logger.info('CALCULATING NETWORK PROPERTIES')
            nwlib.calculate_network_properties(self._matrix_rho,self._dict_params['file_lab']);

    def _get_jmatrix_data(self):
        self._jmatrix_df    =   fileio.read_jmatrix(self._dict_params["file_jmat"])
        self._jmatrix_df.loc[:,"Res_a"]    =   fileio.convert_df_col_to_str(self._jmatrix_df,"Res_a")
        self._jmatrix_df.loc[:,"Res_b"]    =   fileio.convert_df_col_to_str(self._jmatrix_df,"Res_b")

    def _process_to_graph(self):
        self.full_graph    =   igh.Graph.DataFrame(self._jmatrix_df,directed=False,use_vids=False)
        self._num_nodes = len(self.full_graph.vs["name"])
        #assign default community id to 0
        self.full_graph.vs["community_id"]=self._num_nodes*[-1]

        #extract indices for each node
        self._dict_node_index={}
        for index in self.full_graph.vs.indices:
            self._dict_node_index[self.full_graph.vs[index]["name"]]=index

    def _process_jmatrix(self):
        '''
            https://www.opentechguides.com/how-to/article/pandas/193/index-slice-subset.html
            https://igraph.org/python/tutorial/develop/tutorials/visualize_communities/visualize_communities.html
        '''
       
        self._vec_Q =   []

        _vec_i      =   self._dict_params['vec_num']
        _choosen_vector = self._jmatrix_df.iloc[:,_vec_i]

        _vec_max    =   int(np.max(_choosen_vector))
        _vec_min    =   int(np.min(_choosen_vector))

        _vec_delta  =   (_vec_max-_vec_min)/self._dict_params['nsteps']
        _list_of_steps  = np.arange(_vec_min,_vec_max+_vec_delta,_vec_delta)
        _cle=""
        self._logger.info("Starting calculation of Communities from the network")
        self._logger.info("%-20s : %3d"%("Analysis vector",self._dict_params['vec_num']))
        self._logger.info("%-20s : %4.2f"%("Modularity (Q) cut-off",self._dict_params['cut_mod']))
        self._logger.info("%-20s : %4.2f"%("Vector Min.",_vec_min))
        self._logger.info("%-20s : %4.2f"%("Vector Max.",_vec_max))
        self._logger.info("%-20s : %4.2f"%("No. of Steps for Search",self._dict_params['nsteps']))


        Q=0; _list_Q=[]
        for i in _list_of_steps:
            _filtered_df    =   self._jmatrix_df[_choosen_vector>i]
            _sliced_df  =   _filtered_df.iloc[:,[0,1,_vec_i]] 
            _mygraph    =   igh.Graph.DataFrame(_sliced_df,directed=False,vertices=None,use_vids=False)
            _cle        =   _mygraph.community_leading_eigenvector()#igh.Graph.community_leading_eigenvector(_mygraph)
            Q           =   _cle.modularity
            self._logger.info("%-20s (Q) : %4.2f (%4.2f)"%("Vector cut-off..",i,Q))

            if(_cle.modularity>self._dict_params['cut_mod']):
                self._logger.info("Desired Q, reached. Saving community graph")
                #self._save_community_graph(_cle,_mygraph);
                self._save_gml(_cle,_mygraph)
                break
        self._logger.info("Desired Q could not be reached. Either lower your Q cut-off or change your vector.")
    def _save_gml(self,community_graph,_sub_graph):

        num_communities=np.max(community_graph.membership)
        _community_id=1
        for _subgraph in community_graph.subgraphs():
            list_of_node_names = _subgraph.vs["name"];
            for _node_name in list_of_node_names:
                _node_index=self._dict_node_index[_node_name]
                self.full_graph.vs[_node_index]["community_id"]=_community_id

            _community_id+=1
        igh.write(self.full_graph,self._dict_params['file_out'])
    def _change_node_name(self,node_name="NaN"):
        new_node_name=""
        _node_name_int=int(node_name)
        if(_node_name_int<10):
            return _node_name="R_000"+node_name
        if(_node_name_int>=10 and _node_name_int<=99):
            return _node_name="R_00"+node_name
        if(_node_name_int>=100 and _node_name_int<=999):
            return _node_name="R_"+node_name

    def _save_community_graph(self,community_graph,main_graph):
        visual_style = {}
        visual_style["vertex_size"] = 25
        visual_style["vertex_label_size"] = 10
        visual_style["edge_width"]=1
        visual_style["vertex_frame_width"]=0
        visual_style["vertex_label_color"] = "black"
        num_communities=np.max(community_graph.membership)
        palette = igh.RainbowPalette(n=num_communities+1)
        for i, community in enumerate(community_graph):
            s1=community_graph.subgraph(i)
            vertex_labels=[int(x) for x in s1.vs["name"]]
            main_graph.vs[community]["color"] = i
            #mygraph.vs[community]["name"] = vertex_labels
            main_graph.vs[community]["vertex_label"] = vertex_labels
            main_graph.vs[community]["vertex_color"] = palette[i]
        visual_style["vertex_label"] =main_graph.vs["vertex_label"]

        igh.plot(main_graph,self._dict_params['file_out'],palette=palette,**visual_style,)
        self._logger.info("Community saved in %s"%(self._dict_params['file_out']))
    def manager(self,dict_params):
        _start   =   timeit.default_timer();
        self._dict_params   =   dict_params;
        self._get_jmatrix_data()
        self._process_to_graph()
        self._process_jmatrix()
        _stop    =   timeit.default_timer();
        self._logger.info('%-20s : %.2f (s)'%('FINISHED_IN',_stop-_start));
