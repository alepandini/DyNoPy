import igraph as ig
import numpy as np


class GMLStats:
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

    def get_weighted_adjacency(self, property_name):
        adjacency_igraph_matrix = self.graph.get_adjacency(attribute = property_name)
        adjacency_matrix = np.array(adjacency_igraph_matrix.data)
        return adjacency_matrix

    def get_row_sum(self, property_name):
        adjacency_matrix = self.get_weighted_adjacency(property_name)
        row_sum = adjacency_matrix.sum(axis = 1)
        return row_sum

    def get_col_sum(self, property_name):
        adjacency_matrix = self.get_weighted_adjacency(property_name)
        col_sum = adjacency_matrix.sum(axis = 0)
        return col_sum
