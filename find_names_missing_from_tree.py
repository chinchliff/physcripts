#!/usr/bin/env python

if __name__ == '__main__':

    import newick3, os, sys

    if len(sys.argv) < 3:
        sys.exit("usage: find_names_missing_from_tree.py <datafile> <treefile>")

    data_file_name = sys.argv[1]
    tree_file_name = sys.argv[2]

    names = set()
    with open(tree_file_name, "r") as treefile:
        tree = newick3.parse(treefile.readline())

    for l in tree.leaves():
        names.add(l.label)

    with open(data_file_name, "r") as datafile:
        for line in datafile:
            parts = line.split()
            if len(parts) > 1 and parts[0] not in names:
                print parts[0]