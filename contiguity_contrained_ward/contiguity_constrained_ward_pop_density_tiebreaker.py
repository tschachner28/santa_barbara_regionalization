# Conducts regionalization using the Contiguity Constrained Ward algorithm described in Guo 2009.
# In tiebreaker cases where there are multiple shortest edges, the edge connecting the regions that are most similar in population density is removed.

import networkx as nx
import pandas as pd

# Helper function to calculate within-region SSD (for one region)
# Input: r (list of IDs in the region)
# Output: SSD of r
def get_SSD_one_region(r, political_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict):
    ssd = 0

    r_pol_sum = sum(political_data_dict[point] for point in r) if len(r) > 1 else political_data_dict[r[0]]
    r_pol_mean = r_pol_sum / len(r)
    r_ext_sum = sum(ext_data_dict[point] for point in r) if len(r) > 1 else ext_data_dict[r[0]]
    r_ext_mean = r_ext_sum / len(r)
    r_agr_sum = sum(agr_data_dict[point] for point in r) if len(r) > 1 else agr_data_dict[r[0]]
    r_agr_mean = r_agr_sum / len(r)
    r_con_sum = sum(con_data_dict[point] for point in r) if len(r) > 1 else con_data_dict[r[0]]
    r_con_mean = r_con_sum / len(r)
    r_neu_sum = sum(neu_data_dict[point] for point in r) if len(r) > 1 else neu_data_dict[r[0]]
    r_neu_mean = r_neu_sum / len(r)
    r_ope_sum = sum(ope_data_dict[point] for point in r) if len(r) > 1 else ope_data_dict[r[0]]
    r_ope_mean = r_ope_sum / len(r)

    for point in r:
        ssd += pow(political_data_dict[point] - r_pol_mean, 2)
        ssd += pow(ext_data_dict[point] - r_ext_mean, 2)
        ssd += pow(agr_data_dict[point] - r_agr_mean, 2)
        ssd += pow(con_data_dict[point] - r_con_mean, 2)
        ssd += pow(neu_data_dict[point] - r_neu_mean, 2)
        ssd += pow(ope_data_dict[point] - r_ope_mean, 2)

    return ssd

# Helper function to calculate across-region SSD (for two regions)
# Inputs: r1 (list of IDs in first region), r2 (list of IDs in second region)
# Outputs: SSD of r1 and r2, SSD of r1, SSD of r2 (all floats)
def get_SSD_two_regions(r1, r2, political_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict):
    regions = [r1, r2]
    regions_str = ['r1', 'r2']
    ssd_list = [0, 0] # [r1_ssd, r2_ssd]
    sums_dict = {"r1": [], "r2": []} # stores sums of values for each region. List for each is [rx_ext_sum, rx_agr_sum, rx_con_sum, rx_neu_sum, rx_con_sum]

    # Calculate SSD of each region
    for i, r in enumerate(regions):
        r_pol_sum = 0
        r_ext_sum = 0
        r_agr_sum = 0
        r_con_sum = 0
        r_neu_sum = 0
        r_ope_sum = 0
        len_r = len(r)
        if len_r > 1:
            for point in r:
                r_pol_sum += political_data_dict[point]
                r_ext_sum += ext_data_dict[point]
                r_agr_sum += agr_data_dict[point]
                r_con_sum += con_data_dict[point]
                r_neu_sum += neu_data_dict[point]
                r_ope_sum += ope_data_dict[point]
        else:
            r_pol_sum += political_data_dict[r[0]]
            r_ext_sum += ext_data_dict[r[0]]
            r_agr_sum += agr_data_dict[r[0]]
            r_con_sum += con_data_dict[r[0]]
            r_neu_sum += neu_data_dict[r[0]]
            r_ope_sum += ope_data_dict[r[0]]

        sums_dict[regions_str[i]].append(r_pol_sum)
        r_pol_mean = r_pol_sum / len_r

        sums_dict[regions_str[i]].append(r_ext_sum)
        r_ext_mean = r_ext_sum / len_r

        sums_dict[regions_str[i]].append(r_agr_sum)
        r_agr_mean = r_agr_sum / len_r

        sums_dict[regions_str[i]].append(r_con_sum)
        r_con_mean = r_con_sum / len_r

        sums_dict[regions_str[i]].append(r_neu_sum)
        r_neu_mean = r_neu_sum / len_r

        sums_dict[regions_str[i]].append(r_ope_sum)
        r_ope_mean = r_ope_sum / len_r

        for point in r:
            ssd_list[i] += pow(political_data_dict[point] - r_pol_mean, 2)
            ssd_list[i] += pow(ext_data_dict[point] - r_ext_mean, 2)
            ssd_list[i] += pow(agr_data_dict[point] - r_agr_mean, 2)
            ssd_list[i] += pow(con_data_dict[point] - r_con_mean, 2)
            ssd_list[i] += pow(neu_data_dict[point] - r_neu_mean, 2)
            ssd_list[i] += pow(ope_data_dict[point] - r_ope_mean, 2)

    # Calculate SSD of the union of the regions
    ssd_r1_r2 = 0
    r1_r2_pol_mean = (sums_dict["r1"][0] + sums_dict["r2"][0]) / (len(r1) + len(r2)) # mean for pol iden
    r1_r2_ext_mean = (sums_dict["r1"][1] + sums_dict["r2"][1]) / (len(r1) + len(r2))  # mean for ext
    r1_r2_agr_mean = (sums_dict["r1"][2] + sums_dict["r2"][2]) / (len(r1) + len(r2))  # mean for agr
    r1_r2_con_mean = (sums_dict["r1"][3] + sums_dict["r2"][3]) / (len(r1) + len(r2))  # mean for con
    r1_r2_neu_mean = (sums_dict["r1"][4] + sums_dict["r2"][4]) / (len(r1) + len(r2))  # mean for neu
    r1_r2_ope_mean = (sums_dict["r1"][5] + sums_dict["r2"][5]) / (len(r1) + len(r2))  # mean for ope

    for point in r1:
        ssd_r1_r2 += pow(political_data_dict[point] - r1_r2_pol_mean, 2)
        ssd_r1_r2 += pow(ext_data_dict[point] - r1_r2_ext_mean, 2)
        ssd_r1_r2 += pow(agr_data_dict[point] - r1_r2_agr_mean, 2)
        ssd_r1_r2 += pow(con_data_dict[point] - r1_r2_con_mean, 2)
        ssd_r1_r2 += pow(neu_data_dict[point] - r1_r2_neu_mean, 2)
        ssd_r1_r2 += pow(ope_data_dict[point] - r1_r2_ope_mean, 2)

    for point in r2:
        ssd_r1_r2 += pow(political_data_dict[point] - r1_r2_pol_mean, 2)
        ssd_r1_r2 += pow(ext_data_dict[point] - r1_r2_ext_mean, 2)
        ssd_r1_r2 += pow(agr_data_dict[point] - r1_r2_agr_mean, 2)
        ssd_r1_r2 += pow(con_data_dict[point] - r1_r2_con_mean, 2)
        ssd_r1_r2 += pow(neu_data_dict[point] - r1_r2_neu_mean, 2)
        ssd_r1_r2 += pow(ope_data_dict[point] - r1_r2_ope_mean, 2)

    return ssd_r1_r2 - ssd_list[0] - ssd_list[1], ssd_list[0], ssd_list[1]

# Input: new_r (tuple of IDs), connecting_id (tuple), r_u_region and r_v_region (each is list of lists)
# Output: the SSD of the smallest first order edge (float), ID that this edge connects to in r1, ID that this edge connects to in r2
def get_min_ssd_and_first_order_edge(new_r, connecting_id, C, R, political_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict):
    min_ssd = float('inf')
    id1 = None
    id2 = None
    connecting_id_region = [r for r in R if connecting_id[0] in r][0]
    for id_r1 in new_r:
        for id_r2 in connecting_id_region:
            # For each contiguous pair of IDs, one from each region, choose the pair with the min SSD
            if id_r1 in C[id_r2] or id_r2 in C[id_r1]:
                ssd, ssd_r1, ssd_r2 = get_SSD_two_regions([id_r1], [id_r2], political_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict)
                if ssd < min_ssd:
                    min_ssd = ssd
                    id1 = id_r1
                    id2 = id_r2
    return min_ssd, id1, id2

# Output: SSD (float)
def get_ssd_current_regions(regions, political_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict):
    ssd = 0

    # Calculate SSD of each region
    for i, r in enumerate(regions):
        r_pol_sum = sum(political_data_dict[point] for point in r) if len(r) > 1 else political_data_dict[r[0]]
        r_pol_mean = r_pol_sum / len(r)

        r_ext_sum = sum(ext_data_dict[point] for point in r) if len(r) > 1 else ext_data_dict[r[0]]
        r_ext_mean = r_ext_sum / len(r)

        r_agr_sum = sum(agr_data_dict[point] for point in r) if len(r) > 1 else agr_data_dict[r[0]]
        r_agr_mean = r_agr_sum / len(r)

        r_con_sum = sum(con_data_dict[point] for point in r) if len(r) > 1 else con_data_dict[r[0]]
        r_con_mean = r_con_sum / len(r)

        r_neu_sum = sum(neu_data_dict[point] for point in r) if len(r) > 1 else neu_data_dict[r[0]]
        r_neu_mean = r_neu_sum / len(r)

        r_ope_sum = sum(ope_data_dict[point] for point in r) if len(r) > 1 else ope_data_dict[r[0]]
        r_ope_mean = r_ope_sum / len(r)

        for point in r:
            ssd += pow(political_data_dict[point] - r_pol_mean, 2)
            ssd += pow(ext_data_dict[point] - r_ext_mean, 2)
            ssd += pow(agr_data_dict[point] - r_agr_mean, 2)
            ssd += pow(con_data_dict[point] - r_con_mean, 2)
            ssd += pow(neu_data_dict[point] - r_neu_mean, 2)
            ssd += pow(ope_data_dict[point] - r_ope_mean, 2)

    return ssd


# Find the two regions with the min average population density difference (i.e. min(avg_pop_dens_ru - avg_pop_dens_rv)
# Inputs: e_star_ru_list: list of potential r_u's (list of tuples), e_star_rv_list: list of potential r_v's (list of tuples)
# Outputs: r_u and r_v, each of which is a list of tuples
def choose_most_similar_pop_density(e_star_ru_list, e_star_rv_list, dens_data_dict):
    print("Breaking tie between: " + str(e_star_ru_list) + " and " + str(e_star_rv_list))
    r_u = None
    r_v = None
    min_diff = float('inf')
    for i, r in enumerate(e_star_ru_list): # iterate through each pair of potential r_u and r_v's
        avg_pop_dens_ru = sum([dens_data_dict[point] for point in r]) / len(r) if len(r) > 1 else dens_data_dict[r[0]]
        avg_pop_dens_rv = sum([dens_data_dict[e_star_rv_list[i][j]] for j in range(0,len(e_star_rv_list[i]))]) / len(e_star_rv_list[i]) if len(r) > 1 else dens_data_dict[e_star_rv_list[i][0]]
        if abs(avg_pop_dens_ru - avg_pop_dens_rv) < min_diff:
            min_diff = abs(avg_pop_dens_ru - avg_pop_dens_rv)
            r_u = r
            r_v = e_star_rv_list[i]
    print("Chose: " + str(r_u) + " and " + str(r_v))
    return r_u, r_v


# Executes Ward's algorithm
# Input: spatially_contiguous, a boolean
def wards_algorithm(first_order_edges=False):
    # Load data from .csv files
    variables = pd.read_csv('all6variables_regionalization_final.xlsx - ALL6VARIABLES.csv', delimiter=',')
    contiguity_data = pd.read_csv('all6variables_regionalization_final.xlsx - CONTIGUITY.csv', delimiter=',')

    # Output files
    data_filename = ''
    final_regions_file = ''
    if not first_order_edges:
        data_filename = 'regionalization_pop_density_tiebreaker.txt'
        final_regions_file = 'final_regions_pop_density_tiebreaker.txt'
    else:
        data_filename = 'regionalization_pop_density_tiebreaker_first_order_edges.txt'
        final_regions_file = 'final_regions_pop_density_tiebreaker_first_order_edges.txt'
    headers = 'r_u, r_v, r_u SSD, r_v SSD, sum of within-region SSD\n'  # Column names
    file1 = open(data_filename, "w")
    file1.writelines(headers)


    # Log file
    log_filename = 'log_wards_first_order_edges.txt'
    log_file = open(log_filename, "w")

    # Step 1
    R = pd.DataFrame(variables, columns=['ID#']).values.tolist()  # regions (set of 500 data points)

    # Load data into dictionaries
    # political data
    political_data_dict = {}
    political_data_list = pd.DataFrame(variables, columns=['ID#', 'POL IDEN']).values.tolist()
    for p in political_data_list:
        political_data_dict[p[0]] = p[1]

    # EXT data
    ext_data_dict = {}
    ext_data_list = pd.DataFrame(variables, columns=['ID#', 'AVG_EXTROVERT']).values.tolist()
    for e in ext_data_list:
        ext_data_dict[e[0]] = e[1]

    # AGR data
    agr_data_dict = {}
    agr_data_list = pd.DataFrame(variables, columns=['ID#', 'AVG_AGREEABLE']).values.tolist()
    for a in agr_data_list:
        agr_data_dict[a[0]] = a[1]

    # CON data
    con_data_dict = {}
    con_data_list = pd.DataFrame(variables, columns=['ID#', 'AVG_CONSCIENTIOUS']).values.tolist()
    for c in con_data_list:
        con_data_dict[c[0]] = c[1]

    # NEU data
    neu_data_dict = {}
    neu_data_list = pd.DataFrame(variables, columns=['ID#', 'AVG_NEUROTIC']).values.tolist()
    for n in neu_data_list:
        neu_data_dict[n[0]] = n[1]

    # OPE
    ope_data_dict = {}
    ope_data_list = pd.DataFrame(variables, columns=['ID#', 'AVG_OPEN']).values.tolist()
    for o in ope_data_list:
        ope_data_dict[o[0]] = o[1]

    # Population Density
    dens_data_dict = {}
    dens_data_list = pd.DataFrame(variables, columns=['ID#', 'DENSITY']).values.tolist()
    for d in dens_data_list:
        dens_data_dict[d[0]] = float(d[1].replace(',', ''))

    # Step 2: Initialize graph's edge set to empty set
    G = nx.Graph()
    C_list = pd.DataFrame(contiguity_data, columns=['source_ID', 'nbr_ID']).values.tolist()
    C = {}  # contiguity dictionary
    for c in C_list:
        if c[0] not in C.keys():
            C[c[0]] = [c[1]]
        else:
            C[c[0]].append(c[1])

    # Step 3: Add edges for contiguous regions
    for r_u in R:
        for r_v in R:
            if (r_u[0] in C.keys() and r_v[0] in C[r_u[0]]) or (r_v[0] in C.keys() and r_u[0] in C[r_v[0]]):
                edge_weight, r_u_ssd, r_v_ssd = get_SSD_two_regions(r_u, r_v, political_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict)
                G.add_edge(tuple(r_u), tuple(r_v),
                           weight=edge_weight)  # need to convert r_u and r_v from lists to tuples first

    # Step 4: iteratively remove the minimum length edge to merge IDs into two major regions
    removed_edge_lengths = []
    G_min_edges = nx.Graph()  # This tree is used in Algorithm 2 (tree_partitioning_algorithm/tree_partitioning_algorithm.py). It contains the edges that are removed, as well as the final edge that remains.
    while_loop_repeats = 0
    while len(R) > 2 and len(list(G.edges())) > 1:
        e_star = min([e[2]['weight'] for i, e in enumerate(G.edges(data=True))])  # length of shortest edge
        removed_edge_lengths.append(e_star)
        e_star_r1_list = [val[0] for i, val in enumerate(G.edges(data=True)) if
                          val[2]['weight'] == e_star]  # regions that this shortest edge connects
        e_star_r2_list = [val[1] for i, val in enumerate(G.edges(data=True)) if val[2]['weight'] == e_star]
        r_u = ()
        r_v = ()

        # if multiple shortest edges, choose the one connecting the regions with the most similar population density
        if len(e_star_r1_list) > 1:
            r_u, r_v = choose_most_similar_pop_density(e_star_r1_list, e_star_r2_list, dens_data_dict)
        else:
            r_u = e_star_r1_list[0]
            r_v = e_star_r2_list[0]

        r_u_r_v_ssd = e_star
        r_u_ssd = get_SSD_one_region(r_u, political_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict)
        r_v_ssd = get_SSD_one_region(r_v, political_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict)

        G_min_edges.add_edge(r_u, r_v, weight=G.get_edge_data(r_u, r_v)[
            'weight'])  # add edge to graph that will be used in Algorithm 2
        G.remove_edge(r_u, r_v)  # remove e_star from edge set
        log_file.write("removed " + str(r_u) + ", " + str(r_v) + "\n")

        new_r = None
        if not first_order_edges:
            # Remove the regions the edge connects to from R, and replace it with the union of these regions
            new_r = tuple(r_u) + tuple(r_v)
            R.remove(list(r_u))  # remove r_u from R
            R.remove(list(r_v))  # remove r_v from R
            R.append(list(new_r))  # add r_u = r_u U r_v (union of both) to R
            # Update the id to region mappings
        else:
            # Remove the regions where the IDs the edge connects to are located, and replace it with the union of these regions
            r_u_region = [r for r in R if r_u[0] in r][0]
            r_v_region = [r for r in R if r_v[0] in r][0]
            new_r = r_u_region + r_v_region
            R.remove(r_u_region)
            R.remove(r_v_region)
            R.append(new_r)
            log_file.write("new R: " + str(R) + "\n")
            # Remove any other edges that connected r_u and r_v
            for id1 in r_u_region:
                for id2 in r_v_region:
                    if (tuple([id1]), tuple([id2])) in G.edges():
                        G.remove_edge(tuple([id1]), tuple([id2]))
                        log_file.write("removed " + str(tuple([id1])) + ", " + str(tuple([id2])) + "\n")


        edges_to_append = []  # list of lists, each of which contains [v1, v2, weight] representing the new edge

        # redirect edges incident on r_u to new_r, and edges incident on r_v to new_r, with the new weights
        if not first_order_edges: # Original Ward's algorithm
            for e in list(G.edges()):
                if e[0] == r_u or e[0] == r_v:
                    new_weight, ssd_ru, ssd_rv = get_SSD_two_regions(new_r, e[1], political_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict)
                    G.remove_edge(e[0], e[1])
                    edges_to_append.append([new_r, e[1], new_weight])
                if e[1] == r_u or e[1] == r_v:
                    new_weight, ssd_ru, ssd_rv = get_SSD_two_regions(e[0], new_r, political_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict)
                    G.remove_edge(e[0], e[1])
                    edges_to_append.append([e[0], new_r, new_weight])

        else: # When redirecting edges, add a first order edge between the ID in r_u and the ID in r_v with the smallest SSD between them
            for e in list(G.edges()):
                if e[0] == r_u or e[0] == r_v:
                    if while_loop_repeats == 100:
                        x=0
                    new_weight, id1, id2 = get_min_ssd_and_first_order_edge(new_r, e[1], C, R, political_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict)
                    G.remove_edge(e[0], e[1])
                    log_file.write("removed " + str(e[0]) + ", " + str(e[1]) + "\n")
                    edges_to_append.append([tuple([id1]), tuple([id2]), new_weight])
                    log_file.write("redirected edge " + str(tuple([id1])) + ", " + str(tuple([id2])) + ", weight=" + str(new_weight) + "\n")
                if e[1] == r_u or e[1] == r_v:
                    new_weight, id1, id2 = get_min_ssd_and_first_order_edge(new_r, e[0], C, R, political_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict)
                    G.remove_edge(e[0], e[1])
                    log_file.write("removed " + str(e[0]) + ", " + str(e[1]) + "\n")
                    edges_to_append.append([tuple([id1]), tuple([id2]), new_weight])
                    log_file.write("redirected edge " + str(tuple([id1])) + ", " + str(tuple([id2])) + ", weight=" + str(new_weight) + "\n")

        # Add the new edges to G
        for edge in edges_to_append:
            G.add_edge(edge[0], edge[1], weight=edge[2])

        # output_data = r_u, r_v, r_u_ssd, r_v_ssd, r_u_r_v_ssd
        # file1.write(str(output_data[0]) + ', ' + str(output_data[1]) + ', ' + str(output_data[2]) + ', ' + str(output_data[3]) + ', ' + str(output_data[4]) + '\n')
        regionalization_result_ssd = get_ssd_current_regions(R, political_data_dict, ext_data_dict, agr_data_dict, con_data_dict, neu_data_dict, ope_data_dict)
        output_data = r_u, r_v, r_u_ssd, r_v_ssd, regionalization_result_ssd
        file1.write(str(output_data[0]) + ', ' + str(output_data[1]) + ', ' + str(output_data[2]) + ', ' + str(
            output_data[3]) + ', ' + str(output_data[4]) + '\n')

        while_loop_repeats += 1

    file2 = open(final_regions_file, "w")
    print("Final R: " + str(R))
    file2.write("Final R: " + str(R) + "\n")
    for i, r in enumerate(R):
        print("IDs in r" + str(i) + ": " + str(r))
        file2.write("IDs in r" + str(i) + ": " + str(r) + "\n")
        print("r" + str(i) + " Length: " + str(len(r)))
        file2.write("r" + str(i) + " Length: " + str(len(r)) + "\n")
    if first_order_edges:
        file2.write("All edges removed: " + str(G_min_edges.edges(data=True)))

    return G, G_min_edges



# Returns a tree containing the set of edges removed, as well as the single edge remaining in the final cluster
def get_ward_tree():
    G, G_min_edges = wards_algorithm(first_order_edges=True)
    last_edge = list(G.edges())[0]
    last_edge_weight = G.get_edge_data(last_edge[0], last_edge[1])['weight']
    G_min_edges.add_edge(last_edge[0], last_edge[1], weight=last_edge_weight) # add final edge remaining in G
    print("Edges for tree in algorithm 2: " + str(list(G_min_edges.edges())))
    return G_min_edges


if __name__ == '__main__':
    #wards_algorithm()
    get_ward_tree()