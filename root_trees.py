#!/usr/bin/env python
"""Will root a set of newick trees using a line-delimited list of taxon names to be included in the outgroup."""

if __name__ == "__main__":

    import phylo3, newick3
    import sys

    if len(sys.argv) < 3:
        print __doc__
        print "usage: roottrees <treesfile> <outgroupsfile>"
        sys.exit(0)

    treesfname = sys.argv[1]
    outgroupsfname = sys.argv[2]

    treesfile = open(treesfname,"r")

    outgroupsfile = open(outgroupsfname,"r")
    outgroup_names = [line.strip() for line in outgroupsfile.readlines()]

    rooted_trees = []
    for line in treesfile:

        tree = newick3.parse(line)

        outgroup = phylo3.getMRCA(tree, outgroup_names)
        rooted_tree = phylo3.reroot(tree, outgroup)

        rooted_trees.append(rooted_tree)

    outfile = open(treesfname.rsplit(".tre",1)[0]+".rooted.tre","w")
    for tree in rooted_trees:
        outfile.write(newick3.to_string(tree)+";\n")
