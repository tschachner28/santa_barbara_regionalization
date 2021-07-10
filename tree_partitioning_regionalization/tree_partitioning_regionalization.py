from contiguity_contrained_ward import *
import sys
import os
import networkx as nx

# Calculates hg*, the homogeneity gain, of a tree
# Input: An nx graph
# Output: hg_star (a float)
def get_hg(tree):
    hg_star = 0
    e_star = None
    for e in list(tree.edges(data=True)):



sys.path.append(os.path.abspath('../contiguity_constrained_ward/contiguity_constrained_ward.py'))
from contiguity_contrained_ward.contiguity_constrained_ward_pop_density_tiebreaker import get_ward_tree
ward_tree = get_ward_tree() # build tree using output of first algorithm

# Step 1: Calculate homogeneity gain of initial tree (ward_tree)


x=0

