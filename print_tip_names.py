#!/usr/bin/env python

if __name__ == '__main__':
    
    import newick3, phylo3, sys

    if len(sys.argv) < 2:
        print "usage: print_tip_names <treefile>"
        sys.exit()

    treefname = sys.argv[1]
    treefile = open(treefname, "r")

    for line in treefile:

        tree = newick3.parse(line)
        for name in [t.label for t in tree.leaves()]:
            print name
