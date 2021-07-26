from contiguity_contrained_ward import *
import sys
import os
import pandas as pd
import networkx as nx
sys.path.append(os.path.abspath('../contiguity_constrained_ward/contiguity_constrained_ward.py'))
from contiguity_contrained_ward.contiguity_constrained_ward_pop_density_tiebreaker import get_ward_tree, get_SSD_one_region, get_ssd_current_regions, choose_most_dissimilar_pop_density

# Calculates hg*, the homogeneity gain, of a tree
# Input: An nx graph
# Output: hg_star (a float), e_star (tuple of 2 ints), Ta and Tb (nx graphs resulting when the best edge e* is removed)
def get_hg(tree, pol_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict, dens_data_dict, min_region_size):
    hg_star = -1 * float('inf')
    e_star = None
    D = list(tree.edges(data=True))
    D_a = None
    D_b = None

    r_u_ties = [] # list of tuples, each of which is a single ID (int)
    r_v_ties = []
    best_edges_trees = {} # key: a potential best edge (tuple of 2 ints), value: trees resulting when this edge is removed (list of nx graphs)

    for_loop_iterations = 0
    for edge in D:
        for_loop_iterations += 1
        if edge[0] == 74 and edge[1] == 399:
            x=0
        tree_temp = tree.copy()
        tree_temp.remove_edge(edge[0], edge[1])
        D_ia, D_ib = get_2_new_trees(tree_temp)
        #if len(D_ia.nodes()) < 8 or len(D_ib.nodes()) < 8: # subtrees must have a minimum of 8 nodes
        #    continue
        #if len(D_ib.nodes()) == 0:
        #    x=0
        h_D = get_SSD_one_region(list(tree_temp.nodes()), pol_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict)
        h_D_ia = get_SSD_one_region(list(D_ia.nodes()), pol_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict)
        h_D_ib = get_SSD_one_region(list(D_ib.nodes()), pol_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict)
        hg = h_D - h_D_ia - h_D_ib
        if hg > hg_star: # potential best edge found
            if len(D_ia.nodes()) >= 8 and len(D_ib.nodes()) >= 8:  # subtrees must have a minimum of 8 nodes
                hg_star = hg
                e_star = (edge[0], edge[1])
                D_a = D_ia
                D_b = D_ib
                r_u_ties = [tuple([edge[0]])] # remove previous potential r_us that have smaller hg
                r_v_ties = [tuple([edge[1]])]
                best_edges_trees = {e_star: [D_ia, D_ib]}
        elif hg == hg_star: # tie found
            if len(D_ia.nodes()) >= 8 or len(D_ib.nodes()) >= 8:  # subtrees must have a minimum of 8 nodes
                r_u_ties.append(tuple([edge[0]]))
                r_v_ties.append(tuple([edge[1]]))
                best_edges_trees[(edge[0], edge[1])] = [D_ia, D_ib]

    if len(r_u_ties) > 1 and len(r_v_ties) > 1: # Need to break tie
        e1, e2 = choose_most_dissimilar_pop_density(r_u_ties, r_v_ties, dens_data_dict)
        e_star = (e1[0], e2[0])
        D_a = best_edges_trees[e_star][0]
        D_b = best_edges_trees[e_star][1]

    #if D_a == None or D_b == None:
    #    x=0
    return hg_star, e_star, D_a, D_b


# Input: tree (nx graph), with the edge removed
# Output: new_tree1, new_tree2 -- nxgraphs
def get_2_new_trees(nx_graph):
    nx_graph_adj = dict(nx_graph.adjacency()) # key: ID, value: {ID1: {'weight': int}, ID2: {'weight': int}, ...}
    #new_tree1 = []
    new_tree1 = nx.Graph()
    appended = {node: False for node in nx_graph.nodes()} # dict with key: ID, value: boolean

    # Add first node to new_tree1, then add nodes adjacent to it, and the nodes adjacent to those, and so on
    first_node = list(nx_graph.nodes())[0]
    #new_tree1.append(first_node)
    new_tree1.add_node(first_node)
    appended[first_node] = True
    curr_ind = 0
    while curr_ind < len(new_tree1): # while there are still adjacent edges to add, add adjacent nodes using breadth first search
        #for node in nx_graph_adj[new_tree1[curr_ind]]:
        for node in nx_graph_adj[list(new_tree1.nodes())[curr_ind]]:
            if not appended[node]:
                #new_tree1.append(node)
                new_tree1.add_edge(list(new_tree1.nodes())[curr_ind], node, weight=nx_graph_adj[list(new_tree1.nodes())[curr_ind]][node]['weight'])
                appended[node] = True
        curr_ind += 1

    # The other tree consists of the remaining IDs
    #new_tree2 = [node for node in nx_graph.nodes() if node not in new_tree1]
    new_tree2 = nx.Graph()
    for node in nx_graph.nodes():
        if node not in new_tree1.nodes():
            if len(nx_graph_adj[node]) == 0:
                new_tree2.add_node(node)
            else:
                for adj_node in nx_graph_adj[node]:
                    new_tree2.add_edge(node, adj_node, weight=nx_graph_adj[node][adj_node]['weight'])

    return new_tree1, new_tree2

def tree_partitioning_regionalization(attrs=['POL IDEN', 'AVG_EXTROVERT', 'AVG_AGREEABLE', 'AVG_CONSCIENTIOUS', 'AVG_NEUROTIC', 'AVG_OPEN'], min_region_size=None):
    ward_tree = get_ward_tree(vars_used=attrs) # build tree using output of first algorithm

    # Load data from .csv files
    variables = pd.read_csv('all6variables_regionalization_final.xlsx - ALL6VARIABLES.csv', delimiter=',')
    contiguity_data = pd.read_csv('all6variables_regionalization_final.xlsx - CONTIGUITY.csv', delimiter=',')

    R = pd.DataFrame(variables, columns=['ID#']).values.tolist()  # regions (set of 500 data points)

    # Load data into dictionaries
    # political data
    pol_data_dict = {}
    if 'POL IDEN' in attrs:
        pol_data_list = pd.DataFrame(variables, columns=['ID#', 'POL IDEN']).values.tolist()
        for p in pol_data_list:
            pol_data_dict[p[0]] = p[1]

    # EXT data
    ext_data_dict = {}
    if 'AVG_EXTROVERT' in attrs:
        ext_data_list = pd.DataFrame(variables, columns=['ID#', 'AVG_EXTROVERT']).values.tolist()
        for e in ext_data_list:
            ext_data_dict[e[0]] = e[1]

    # AGR data
    agr_data_dict = {}
    if 'AVG_AGREEABLE' in attrs:
        agr_data_list = pd.DataFrame(variables, columns=['ID#', 'AVG_AGREEABLE']).values.tolist()
        for a in agr_data_list:
            agr_data_dict[a[0]] = a[1]

    # CON data
    con_data_dict = {}
    if 'AVG_CONSCIENTIOUS' in attrs:
        con_data_list = pd.DataFrame(variables, columns=['ID#', 'AVG_CONSCIENTIOUS']).values.tolist()
        for c in con_data_list:
            con_data_dict[c[0]] = c[1]

    # NEU data
    neu_data_dict = {}
    if 'AVG_NEUROTIC' in attrs:
        neu_data_list = pd.DataFrame(variables, columns=['ID#', 'AVG_NEUROTIC']).values.tolist()
        for n in neu_data_list:
            neu_data_dict[n[0]] = n[1]

    # OPE
    ope_data_dict = {}
    if 'AVG_OPEN' in attrs:
        ope_data_list = pd.DataFrame(variables, columns=['ID#', 'AVG_OPEN']).values.tolist()
        for o in ope_data_list:
            ope_data_dict[o[0]] = o[1]

    # Population Density
    dens_data_dict = {}
    dens_data_list = pd.DataFrame(variables, columns=['ID#', 'DENSITY']).values.tolist()
    for d in dens_data_list:
        dens_data_dict[d[0]] = float(d[1].replace(',', ''))

    # Output file
    tree_partitions_filename = 'tree_partitions'
    vars_abbreviated = {'POL IDEN': 'POL', 'AVG_EXTROVERT': 'EXT', 'AVG_AGREEABLE': 'AGR', 'AVG_CONSCIENTIOUS': 'CON',
                        'AVG_NEUROTIC': 'NEU', 'AVG_OPEN': 'OPE'}
    if len(attrs) < 6:
        for variable in attrs:
            tree_partitions_filename += '_' + vars_abbreviated[variable]
    if min_region_size != None:
        tree_partitions_filename += '_min_region_size_' + str(min_region_size)
    tree_partitions_filename += '.txt'
    tree_partitions_file = open(tree_partitions_filename, "w")

    # Step 1: Calculate homogeneity gain of initial tree (ward_tree), and the best tree with the largest hg
    hg_star, e_star, D_a, D_b = get_hg(ward_tree, pol_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict, dens_data_dict, min_region_size)
    trees = {ward_tree: hg_star} # key: tree (nxgraph), value: float. Denoted R in the article
    for k in range(2, 31):
        while len(list(trees.keys())) < k:
            # Step 2: Find the best tree with the largest hg*
            best_tree = [t for t in list(trees.keys()) if trees[t] == max(list(trees.values()))][0]
            # Step 3: Cut best_tree into two trees by removing the best edge
            hg_star, e_star, T_a, T_b = get_hg(best_tree, pol_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict, dens_data_dict, min_region_size)
            if T_a == None or T_b == None: # enters this statement if no cut found
                raise Exception('No cut found')
            # Step 4: Calculate hg* for T_a and T_b
            hg_star_T_a, e_star_T_a, T_a_a, T_a_b = get_hg(T_a, pol_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict, dens_data_dict, min_region_size)
            hg_star_T_b, e_star_T_b, T_b_a, T_b_b = get_hg(T_b, pol_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict, dens_data_dict, min_region_size)
            del trees[best_tree]
            trees[T_a] = hg_star_T_a
            trees[T_b] = hg_star_T_b

        tree_partitions_file.write("k = " + str(k) + "\n")
        for tree in trees.keys():
            tree_partitions_file.write(str(tree.nodes()) + "\n")
        current_regions = [list(tree.nodes()) for tree in trees.keys()]
        within_region_ssd = get_ssd_current_regions(current_regions, pol_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict)
        tree_partitions_file.write("Within-Region SSD: " + str(within_region_ssd) + "\n")
        tree_partitions_file.write("\n")


if __name__ == "__main__":
    tree_partitioning_regionalization(min_region_size=8)
    #tree_partitioning_regionalization()
