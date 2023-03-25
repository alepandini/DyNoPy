import igraph as ig
import re


class GMLParser:
    def __init__(self, gml_filename):
        self.gml_filename = gml_filename
        self.graph = self.__read_gml_file()

    def __read_gml_file(self):
        graph = ig.read(self.gml_filename)
        return graph

    def __get_res_number(self, resid):
        res_number = None
        resid_match = re.search("^(R_)?(\d+)", resid)
        if resid_match is not None:
            res_number = int(resid_match.groups()[1])
            return res_number
        else:
            return None

    def create_nodes_dict(self, property_name, integer_cast = False):
        nodes_dict = dict()
        for vertex in self.graph.vs:
            id_key = self.__get_res_number(vertex.attributes()['name'])
            if integer_cast:
                nodes_dict[id_key] = int(vertex.attributes()[property_name])
            else:
                nodes_dict[id_key] = vertex.attributes()[property_name]
        return nodes_dict

    def create_edges_dict(self, property_name, integer_cast = False):
        edges_dict = dict()
        for edge in self.graph.es:
            s_id = edge.source
            t_id = edge.target
            s_id_key = self.__get_res_number(self.graph.vs[s_id].attributes()['name'])
            t_id_key = self.__get_res_number(self.graph.vs[t_id].attributes()['name'])
            if integer_cast:
                edges_dict[(s_id_key, t_id_key)] = int(edge.attributes()[property_name])
            else:
                edges_dict[(s_id_key, t_id_key)] = edge.attributes()[property_name]
        return edges_dict

    def generate_groups_dict(self, property_name, integer_cast = True):
        group_dict = {}
        nodes_dict = self.create_nodes_dict(property_name, integer_cast)
        value_list = list(nodes_dict.values())
        key_list = list(nodes_dict.keys())
        for v in set(value_list):
            k_indices = [i for i, x in enumerate(value_list) if x == v]
            v_key_list = [key_list[i] for i in k_indices]
            group_dict[v] = v_key_list
        return group_dict
