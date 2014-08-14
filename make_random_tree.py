#!/usr/bin/env python

import newick3, phylo3, argparse, sys, random

if __name__ == "__main__":

    description = "make a random tree with X tips"

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-x", "--number-of-tips", type=int, nargs=1, required=False, help="The number of tips. If this is specified, tip labels will be numbers. Either this or a set of input tip labels must be specified.")

    parser.add_argument("-l", "--tip-labels", nargs='*', required=False, help="A set of tip labels for the random tree. If this is not set, then the number of tips must be specified, and tip labels will be numbers.")

    parser.add_argument("-m", "--max-branch-length", type=float, nargs=1, required=False, help="The upper bound for creating randomized branch lengths. If not set, branch lengths will be zero. This argument cannot be used in conjuction with the TREE_AGE argument.")

#    parser.add_argument("-a", "--tree_age", type=float, nargs=1, required=False, help="The age of the tree. If this is set, then an ultrametric tree will be produced, with the age of the root node equal to the value of this argument. This argument cannot be used in conjunction with MAX_BRANCH_LENGTH.")

    args = parser.parse_args()
    
    if (args.number_of_tips != None and args.tip_labels != None) or (args.number_of_tips == None and args.tip_labels == None):
        sys.exit("You must specify either the number of tips or the tip labels, but not both")

    tip_labels = args.tip_labels if args.tip_labels != None else [str(i) for i in range(1, args.number_of_tips[0]+1)]

#    if (args.max_branch_length != None and args.tree_age != None):
#        sys.exit("You may only specify one of the MAX_BRANCH_LENGTH or TREE_AGE arguments, not both")

    max_brlen = args.max_branch_length[0] if args.max_branch_length != None else 0
    make_brlens = True if max_brlen > 0 else False

#    root_age = args.tree_age[0] if args.tree_age != None else 0
#    make_ultrametric = True if root_age > 0 else False
    
    nodes = []
    for l in tip_labels:
        t = phylo3.Node()
        t.label = l
        t.istip = True
        if make_brlens:
            t.length = random.random() * max_brlen
        nodes.append(t)
    
    while len(nodes) > 1:
        random.shuffle(nodes)
        new_node = phylo3.Node()
        new_node.add_child(nodes.pop())
        new_node.add_child(nodes.pop())
        if make_brlens:
            new_node.length = random.random() * max_brlen
        nodes.append(new_node)
    
    print(newick3.to_string(nodes[0], use_branch_lengths=make_brlens)+";")