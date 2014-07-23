#!/usr/bin/env python

import sys, newick3, phylo3

if len(sys.argv) < 3:
    sys.exit("usage: extract_subtree.py <treefile> <mrca1> <mrca2>")

tree_file_name = sys.argv[1]
mrca1 = sys.argv[2].strip()
mrca2 = sys.argv[3].strip()

with open(tree_file_name) as treefile:
    tree = newick3.parse(treefile.readline())

print newick3.to_string(phylo3.get_mrca(tree, [mrca1, mrca2]))