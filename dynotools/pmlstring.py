from dynotools.gmlparser import *
import igraph as ig
import matplotlib as mpl
import numpy as np

class PMLString:
    def __init__(self, gml_parser):
        self.gmlparser = gml_parser

    def __colour_gradient(self, start_colour, end_colour, step, n_steps = 100):
        s_rgb_c = mpl.colors.to_rgb(start_colour)
        e_rgb_c = mpl.colors.to_rgb(end_colour)
        f = (step / n_steps)
        o_hex_c = mpl.colors.to_hex([(1-f)*s_rgb_c[i] + f*e_rgb_c[i] for i in range(len(s_rgb_c))])
        return o_hex_c

    def generate_group_selection(self, property_name):
        group_selection_strings = []
        group_dict = self.gmlparser.generate_groups_dict(property_name)
        key_list = list(group_dict.keys())
        key_list.sort()
        selection_template = "select {0}_{1}, resi {2}"
        for k in key_list:
            res_string = '+'.join([str(res_id) for res_id in group_dict[k]])
            group_selection_strings.append(selection_template.format(property_name, k, res_string))
        group_selection_strings.append("deselect")
        return group_selection_strings

    def alter_bfactor_by_value(self, property_name):
        bfactor_alteration_strings = []
        nodes_dict = self.gmlparser.create_nodes_dict(property_name)
        min_value = min(list(nodes_dict.values()))
        max_value = max(list(nodes_dict.values()))
        key_list = list(nodes_dict.keys())
        key_list.sort()
        bfactor_alt_template = "alter resi {0}, b={1:4.2f}"
        for k in key_list:
            b_value = nodes_dict[k]/max_value
            bfactor_alteration_strings.append(bfactor_alt_template.format(k, b_value))
        return bfactor_alteration_strings

    def add_network_edges(self, property_name, color_scale = True, scaling_factor = 10):
        start_colour = "blue"
        end_colour = "orange"
        network_edges_strings = []
        edge_dict = self.gmlparser.create_edges_dict(property_name)
        min_value = min(edge_dict.values())
        max_value = max(edge_dict.values())
        n_values = len(set(edge_dict.values()))
        distance_template = "distance e{0}_{1}, name CA and resi {0}, name CA and resi {1}"
        label_template = "hide label, e{0}_{1}"
        gap_template = "set dash_gap, 0, e{0}_{1}"
        width_template = "set dash_width, {0}, e{1}_{2}"
        color_template = "color 0x{0}, e{1}_{2}"
        for resids, value in edge_dict.items():
            if value > 0:
                weight = value / max_value * scaling_factor
                network_edges_strings.append(distance_template.format(resids[0], resids[1]))
                network_edges_strings.append(label_template.format(resids[0], resids[1]))
                network_edges_strings.append(gap_template.format(resids[0], resids[1]))
                network_edges_strings.append(width_template.format(weight, resids[0], resids[1]))
                if color_scale:
                    i = int(((value - min_value)/(max_value - min_value))*n_values)
                    edge_color = self.__colour_gradient(start_colour, end_colour, i, n_values)
                    network_edges_strings.append(color_template.format(edge_color[1:], resids[0], resids[1]))
        return network_edges_strings

    def save_pml_script(self, pml_cmd_list, outfilename):
        fout = open(outfilename, 'w')
        for pml_cmd in pml_cmd_list:
            fout.write("{0}\n".format(pml_cmd))
        fout.close()
