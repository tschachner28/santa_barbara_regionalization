import networkx as nx
import pandas as pd

# Helper function to calculate within-region SSD (for one region)
# Input: r (list of IDs in the region)
# Output: SSD of r
def get_SSD_one_region(r):
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
def get_SSD_two_regions(r1, r2):
    regions = [r1, r2]
    regions_str = ['r1', 'r2']
    ssd_list = [0, 0] # [r1_ssd, r2_ssd]
    sums_dict = {"r1": [], "r2": []} # stores sums of values for each region. List for each is [rx_ext_sum, rx_agr_sum, rx_con_sum, rx_neu_sum, rx_con_sum]

    # Calculate SSD of each region
    for i, r in enumerate(regions):
        r_pol_sum = sum(political_data_dict[point] for point in r) if len(r) > 1 else political_data_dict[r[0]]
        sums_dict[regions_str[i]].append(r_pol_sum)
        r_pol_mean = r_pol_sum / len(r)
        for point in r:
            ssd_list[i] += pow(political_data_dict[point] - r_pol_mean, 2)

        r_ext_sum = sum(ext_data_dict[point] for point in r) if len(r) > 1 else ext_data_dict[r[0]]
        sums_dict[regions_str[i]].append(r_ext_sum)
        r_ext_mean = r_ext_sum / len(r)
        for point in r:
            ssd_list[i] += pow(ext_data_dict[point] - r_ext_mean, 2)

        r_agr_sum = sum(agr_data_dict[point] for point in r) if len(r) > 1 else agr_data_dict[r[0]]
        sums_dict[regions_str[i]].append(r_agr_sum)
        r_agr_mean = r_agr_sum / len(r)
        for point in r:
            ssd_list[i] += pow(agr_data_dict[point] - r_agr_mean, 2)

        r_con_sum = sum(con_data_dict[point] for point in r) if len(r) > 1 else con_data_dict[r[0]]
        sums_dict[regions_str[i]].append(r_con_sum)
        r_con_mean = r_con_sum / len(r)
        for point in r:
            ssd_list[i] += pow(con_data_dict[point] - r_con_mean, 2)

        r_neu_sum = sum(neu_data_dict[point] for point in r) if len(r) > 1 else neu_data_dict[r[0]]
        sums_dict[regions_str[i]].append(r_neu_sum)
        r_neu_mean = r_neu_sum / len(r)
        for point in r:
            ssd_list[i] += pow(neu_data_dict[point] - r_neu_mean, 2)

        r_ope_sum = sum(ope_data_dict[point] for point in r) if len(r) > 1 else ope_data_dict[r[0]]
        sums_dict[regions_str[i]].append(r_ope_sum)
        r_ope_mean = r_ope_sum / len(r)
        for point in r:
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


# Find the two regions with the min average population density difference (i.e. min(avg_pop_dens_ru - avg_pop_dens_rv)
# Inputs: e_star_ru_list: list of potential r_u's (list of tuples), e_star_rv_list: list of potential r_v's (list of tuples)
# Outputs: r_u and r_v, each of which is a list of tuples
def choose_most_similar_pop_density(e_star_ru_list, e_star_rv_list):
    print("Breaking tie between: " + str(e_star_ru_list) + " and " + str(e_star_rv_list))
    r_u = None
    r_v = None
    min_diff = float('inf')
    for i, r in enumerate(e_star_ru_list): # iterate through each pair of potential r_u and r_v's
        avg_pop_dens_ru = sum([dens_data_dict[id] for id in r]) / len(r)
        avg_pop_dens_rv = sum([dens_data_dict[e_star_rv_list[i][j]] for j in range(0,len(e_star_rv_list[i]))]) / len(r)
        if abs(avg_pop_dens_ru - avg_pop_dens_rv) < min_diff:
            min_diff = abs(avg_pop_dens_ru - avg_pop_dens_rv)
            r_u = r
            r_v = e_star_rv_list[i]
    print("Chose: " + str(r_u) + " and " + str(r_v))
    return r_u, r_v




# Load data from .csv files
variables = pd.read_csv('all6variables_regionalization_final.xlsx - ALL6VARIABLES.csv', delimiter=',')
contiguity_data = pd.read_csv('all6variables_regionalization_final.xlsx - CONTIGUITY.csv', delimiter=',')

# Output files
data_filename = 'regionalization_pop_density_tiebreaker.txt'
headers = 'r_u, r_v, r_u SSD, r_v SSD, r_u and r_v SSD\n'  # Column names
file1 = open(data_filename, "w")
file1.writelines(headers)
final_regions_file = 'final_regions_pop_density_tiebreaker.txt'

#widths = [5,5,5,5,5]
#format_spec = '{:{widths[0]}}  {:>{widths[1]}}  {:>{widths[2]}}  {:>{widths[3]}}  {:>{widths[4]}}'
#print(format_spec.format(*headers, widths=widths))

# Step 1
R = pd.DataFrame(variables, columns= ['ID#']).values.tolist() # regions (set of 500 data points)

# Load data into dictionaries
# political data
political_data_dict = {}
political_data_list = pd.DataFrame(variables, columns= ['ID#', 'POL IDEN']).values.tolist()
for p in political_data_list:
    political_data_dict[p[0]] = p[1]

# EXT data
ext_data_dict = {}
ext_data_list = pd.DataFrame(variables, columns= ['ID#', 'AVG_EXTROVERT']).values.tolist()
for e in ext_data_list:
    ext_data_dict[e[0]] = e[1]

# AGR data
agr_data_dict = {}
agr_data_list = pd.DataFrame(variables, columns= ['ID#', 'AVG_AGREEABLE']).values.tolist()
for a in agr_data_list:
    agr_data_dict[a[0]] = a[1]

# CON data
con_data_dict = {}
con_data_list = pd.DataFrame(variables, columns= ['ID#', 'AVG_CONSCIENTIOUS']).values.tolist()
for c in con_data_list:
    con_data_dict[c[0]] = c[1]

# NEU data
neu_data_dict = {}
neu_data_list = pd.DataFrame(variables, columns= ['ID#', 'AVG_NEUROTIC']).values.tolist()
for n in neu_data_list:
    neu_data_dict[n[0]] = n[1]

# OPE
ope_data_dict = {}
ope_data_list = pd.DataFrame(variables, columns= ['ID#', 'AVG_OPEN']).values.tolist()
for o in ope_data_list:
    ope_data_dict[o[0]] = o[1]

# Population Density
# OPE
dens_data_dict = {}
dens_data_list = pd.DataFrame(variables, columns= ['ID#', 'AVG_OPEN']).values.tolist()
for d in dens_data_list:
    dens_data_dict[d[0]] = d[1]

# Step 2: Initialize graph's edge set to empty set
G = nx.Graph()
C_list = pd.DataFrame(contiguity_data, columns= ['source_ID', 'nbr_ID']).values.tolist()
C = {} # contiguity dictionary
for c in C_list:
    if c[0] not in C.keys():
        C[c[0]] = [c[1]]
    else:
        C[c[0]].append(c[1])


# Step 3: Add edges for contiguous regions
for r_u in R:
    for r_v in R:
        if (r_u[0] in C.keys() and r_v[0] in C[r_u[0]]) or (r_v[0] in C.keys() and r_u[0] in C[r_v[0]]):
            edge_weight, r_u_ssd, r_v_ssd = get_SSD_two_regions(r_u, r_v)
            G.add_edge(tuple(r_u), tuple(r_v), weight=edge_weight) # need to convert r_u and r_v from lists to tuples first

# Step 4: iteratively remove the minimum length edge to merge IDs into two major regions
removed_edge_lengths = []
while_loop_repeats = 0
while len(R) > 2 and len(list(G.edges())) > 1:
    e_star = min([e[2]['weight'] for i, e in enumerate(G.edges(data=True))]) # length of shortest edge
    removed_edge_lengths.append(e_star)
    e_star_r1_list = [val[0] for i,val in enumerate(G.edges(data=True)) if val[2]['weight'] == e_star] # regions that this shortest edge connects
    e_star_r2_list = [val[1] for i, val in enumerate(G.edges(data=True)) if val[2]['weight'] == e_star]
    r_u = ()
    r_v = ()



    # if multiple shortest edges, choose the one connecting the regions with the smallest combined ID
    if len(e_star_r1_list) > 1:
        r_u, r_v = choose_most_similar_pop_density(e_star_r1_list, e_star_r2_list)


    else:
        r_u = e_star_r1_list[0]
        r_v = e_star_r2_list[0]

    r_u_r_v_ssd = e_star
    r_u_ssd = get_SSD_one_region(r_u)
    r_v_ssd = get_SSD_one_region(r_v)


    if while_loop_repeats == 5:
        y = 1

    G.remove_edge(r_u, r_v) # remove e_star from edge set
    new_r = tuple(r_u) + tuple(r_v)
    R.remove(list(r_u)) # remove r_u from R
    R.remove(list(r_v)) # remove r_v from R
    R.append(list(new_r)) # add r_u = r_u U r_v (union of both) to R
    edges_to_append = [] # list of lists, each of which contains [v1, v2, weight] representing the new edge


    # redirect edges incident on r_u to new_r, and edges incident on r_v to new_r, with the new weights
    for e in list(G.edges()):
        if e[0] == r_u or e[0] == r_v:
            new_weight, ssd_ru, ssd_rv = get_SSD_two_regions(new_r, e[1])
            G.remove_edge(e[0], e[1])
            edges_to_append.append([new_r, e[1], new_weight])
        if e[1] == r_u or e[1] == r_v:
            if while_loop_repeats == 432:
                z = 2
            new_weight, ssd_ru, ssd_rv = get_SSD_two_regions(e[0], new_r)
            G.remove_edge(e[0], e[1])
            edges_to_append.append([e[0], new_r, new_weight])

    # Add the new edges to G
    for edge in edges_to_append:
        if while_loop_repeats == 432:
            z = 2
        G.add_edge(edge[0], edge[1], weight=edge[2])

    output_data = r_u, r_v, r_u_ssd, r_v_ssd, r_u_r_v_ssd
    #print(format_spec.format(*output_data, widths=widths))
    file1.write(str(output_data[0]) + ', ' + str(output_data[1]) + ', ' + str(output_data[2]) + ', ' + str(output_data[3]) + ', ' + str(output_data[4]) + '\n')

    while_loop_repeats += 1




x=0

file2 = open(final_regions_file, "w")
print("Final R: " + str(R))
file2.write("Final R: " + str(R) + "\n")
for i, r in enumerate(R):
    print("IDs in r" + str(i) + ": " + str(r))
    file2.write("IDs in r" + str(i) + ": " + str(r) + "\n")
    print("r" + str(i) + " Length: " + str(len(r)))
    file2.write("r" + str(i) + " Length: " + str(len(r)) + "\n")



