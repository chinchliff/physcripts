#!/usr/bin/env python

import argparse, math, newick3, phylo3, sys

def recur_balanced(p, max_d, max_t, _d=0):
    
    global _j
    
    if max_t < 2 or _d > max_d:
        return

    _d += 1
    n_max_t = [math.floor(max_t/2),]
    n_max_t.append(max_t - n_max_t[0])
    for i in range(2):
        c = phylo3.Node()
        p.add_child(c)
        if n_max_t[i] < 2:
            c.istip = True
            c.label = "T"+str(_j)
            _j += 1
        else:
            recur_balanced(c, max_d, n_max_t[i], _d)

def get_balanced_tree(n_tips_per_tree, br_length_function, mean_br_length):

    global _j
    _j = 1
    
    tree_depth = math.floor(math.log(n_tips_per_tree,2))

    tree = phylo3.Node()
    recur_balanced(tree, tree_depth, n_tips_per_tree)
    
    return assign_branch_lengths(tree, br_length_function, mean_br_length)

def get_pectinate_tree(n_tips_per_tree, br_length_function, mean_br_length):

    tree = phylo3.Node()
    n = tree
    for i in range(n_tips_per_tree - 1):
        c = phylo3.Node()
        t = phylo3.Node()
        t.istip = True
        t.label = "T"+str(i+1)
        n.add_child(c)
        n.add_child(t)
        n = c
        if i == n_tips_per_tree - 2:
            c.istip = True
            c.label = "T"+str(i+2)

    return assign_branch_lengths(tree, br_length_function, mean_br_length)
    
def assign_branch_lengths(tree, br_length_function, mean_br_length):

    # NOTE: assumes a bifurcating tree
    for n in tree.iternodes(order=phylo3.POSTORDER):

        if n.istip:
        
            # check if this is a cherry
            p = n.parent
            for s in p.children:
                if s != n:
                    if s.istip and hasattr(s, "height"):
                        # a cherry, other tip has an assigned length, copy it
                        n.length = s.length
                    else: # not a cherry or other tip has no length, assign a length
                        n.length = br_length_function(mean_br_length)
            n.height = 0 # tips are always height 0

        else: # not a tip
        
            # collect lengths and heights of children
            l = (n.children[0].length, n.children[1].length)
            h = (n.children[0].height, n.children[1].height)

            # normalize the lengths of this node's children so the
            # cumulative distance to the tips is equal on both sides
            if h[0] + l[0] > h[1] + l[1]:
                m = 0
                n.children[1].length += h[0] + l[0] - (h[1] + l[1])
            else:
                m = 1
                n.children[0].length += h[1] + l[1] - (h[0] + l[0])
                
            # assign length and height to this node
            n.length = br_length_function(mean_br_length)
            n.height = n.children[m].height + n.children[m].length

    return tree

def validate_tree_type(t):
    if t == 'pectinate':
        return get_pectinate_tree
    elif t == 'balanced':
        return get_balanced_tree
    else:
        sys.exit('invalid tree type: ' + t)

if __name__ == '__main__':

    parser = argparse.ArgumentParser('')
    
    parser.add_argument('-n', '--number-of-tips', type=int, required=True, \
        help='The number of tips to be in the final tree')
    
    parser.add_argument('-b', '--mean-branch-length', type=float, required=True, \
        help='The length to be used for internal branches')
    
    parser.add_argument('-t', '--tree-generator-method', type=validate_tree_type, required=True, \
        help='The type of tree to generate. Must be either "pectinate" or "balanced".')
    
    args = parser.parse_args()
    
    sys.setrecursionlimit(args.number_of_tips+100)
    
    # currently just supports constant rate
    t = args.tree_generator_method(args.number_of_tips, lambda x: x, args.mean_branch_length)

    print(newick3.to_string(t) + ';')