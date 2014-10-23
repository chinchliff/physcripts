#!/usr/bin/env python

import argparse, newick3, phylo3, random

MIN_BRANCH_LENGTH = 0.01

if __name__ == '__main__':

    description = 'add species randomly to an ultrametric chronogram, ' \
                  'while maintaining ultrametricity.'

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-t', '--input-tree', type=argparse.FileType('r'), nargs=1, \
        required=True, help='The chronogram to which species should be added.')

    parser.add_argument('-n', '--names', type=argparse.FileType('r'), nargs=1, \
        required=True, help='The list of names to be added to the tree.')

    # allow a min branch length to be specified. it must be parseable as a float
    parser.add_argument('-b', '--min-branch-length', type=float, nargs=1, \
        required=False, help='The minimum branch length to be used.')

    # record a boolean value of True if this argument is set
    parser.add_argument('-s', '--include-stem', action='store_true', \
        required=False, help='Pass this argument to allow newly added species to be '
                             'attached to the root of the tree.')

    args = parser.parse_args()
    
    # attempt to parse the input tree
    try:
        tree = newick3.parse(args.input_tree[0])
    except Exception as e:
        print("There was a problem parsing the input tree: " + e.message)
        exit(1)
    
    # extract the names from the names file
    names = [n.strip() for n in args.names[0]]

    # use the user-specified min branch length if specified
    min_branch_length = args.min_branch_length[0] if args.min_branch_length is not None \
                                                  else MIN_BRANCH_LENGTH

    # assign the function that will be used to gather the nodes--the iternodes function
    # will include the root, but the descendants function will not
    get_nodes = phylo3.Node.iternodes if args.include_stem else phylo3.Node.descendants

    for name in names:

        # get one node at random from the tree
        n = random.sample(list(get_nodes(tree)), 1)[0]

        # create a new internal node to be the parent of the tip we're about to add
        new_parent = phylo3.Node()

        if n.parent == None: # root case: new parent becomes the new root of the tree
            tree = new_parent

        else: # non-root case: randomly place new parent on the branch leading to n
            p = n.parent
            p.remove_child(n)
            p.add_child(new_parent)

        # reattach n to the new parent
        new_parent.add_child(n)

        # update branch lengths
        new_parent.length = random.uniform(min_branch_length, n.length - min_branch_length)
        n.length = n.length - new_parent.length

        # add the new tip as a child of the new parent
        new_tip = phylo3.Node()
        new_tip.label = name
        new_tip.istip = True
        new_tip.length = new_parent.depth - new_parent.length
        new_parent.add_child(new_tip)
            
    print(newick3.to_string(tree))