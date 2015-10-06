#!/usr/bin/env python3

import argparse, newick3, phylo3, sys, os
import figtree_blocks
from StringIO import StringIO

BLACK = "000000"

def get_bin(value, bins):

    '''Check the incoming bins and return the first bin that includes the value.
    Assumes that bins are sorted in descending order.'''

    is_highest_bin = True
    for bin in bins:
    
        # the first bin is the highest bin
        if is_highest_bin:
            is_highest_bin = False
            
            # upper val is inclusive on highest bin
            if value > bin["upper"]:
            
                # value is higher than upper val of highest bin
                return None

        else:
            # when not on the highest bin, upper val is exclusive
            if value >= bin["upper"]:
                continue

        # lower val on all bins is inclusive
        if value < bin["lower"]:
            continue

        else:
            # returns first matching bin for the value
            return bin

    # value is lower than lower val of lowest bin, OR outside of defined bin values (if bins are discontinuous)
    return None

def make_newick_with_color_strings(node):
    
    '''Return a newick string encoding the subtree below the passed node, with FigTree-format color
    labels defined by the color_string property on the subtree's nodes.'''

    if not node.istip:

        substrings = []
        for child in node.children:
            substrings.append(make_newick_with_color_strings(child))

        return "(" + ",".join(substrings) + ")[&!color=#" + node.color_string + "]:" + str(node.length)

    else:

        return node.label + "[&!color=#" + node.color_string + "]:" + str(node.length)

def write_tree(outfile_name, tree):

    '''Write the tree, painted by the color bins, in the figtree format, to the specified outfile. Calls
    make_newick_with_color_strings, which assumes that the tree nodes contain color_string attributes
    that will be populated into the output.'''

    outfile = open(outfile_name,"w")
    outfile.write("#NEXUS\nbegin taxa\n\tdimensions ntax=" + str(len(tree.leaves())) + ";\n\ttaxlabels\n")
    for tip in tree.leaves():
        outfile.write("\t"+tip.label+"\n")
    outfile.write(";\nend;\n\nbegin trees;\n\ttree tree1 = [&R]\t")

    outfile.write(make_newick_with_color_strings(tree))

    outfile.write(";\nend;\n\n")
    outfile.write(figtree_blocks.circular_sorted_tree)
    outfile.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Annotate a labeled tree with figtree-format branch color labels by binning node values into ranges associated with colors.")

    parser.add_argument("-t", "--tree", type=open, nargs=1, required=True, help="A newick tree file with labeled internal nodes. If a node-values file is supplied, then the node labels in the tree will be used to match to the rows of the node-values file. Otherwise, the node labels themselves will be interpreted as the values to be used for tree painting.")

    parser.add_argument("-c", "--color-bins", type=open, nargs=1, required=True, help="A comma-separated list containing min and max values for the bins to which node values will be assigned, and the color which each branch falling into that bin should be painted.")

    parser.add_argument("-d", "--node-values", type=open, nargs=1, help="A list containing N columns of comma-separated values corresponding to nodes in the tree, where the first column is the node label of some internal node in the tree, and the other columns are values assigned to that node. If a node-values file is supplied, the --label argument must used to specify which value column(s) in the node-values file to be used for tree painting.")

    parser.add_argument("-l", "--label", nargs="+", help="One or more column labels from the node-values file whose values will be used for tree-painting. If multiple values are specified, one output tree will be produced for each column of values, painted according to the color scheme specified by the color-bins file. If no node-values file has been designated, all --label values will be ignored.")
    
    parser.add_argument("-v", "--tip-values", type=open, nargs=1, help="A list containing 2 columns of comma-separated values, where the first item on each line is the tip label and the second is the value to be assigned to the the tip. This may be used to assign values to tip nodes in a tree whose internal nodes are assigned values with node labels.")

    parser.add_argument("-o", "--output-dir", nargs=1, help="A location to which painted trees should be saved. If not specified they will be saved to the current working directory.")
    
    parser.add_argument("-n", "--output-prefix", nargs=1, help="A prefix to be attached to output files. If not specified, no prefix will be used (i.e. preexisting output file will be overwritten).")

    args = parser.parse_args()
    
    target_dir = os.path.abspath(".")
    if args.output_dir != None:
        target_dir = os.path.abspath(args.output_dir)
    print("will save output to: " + target_dir)
    
    output_prefix = args.output_prefix[0] if args.output_prefix is not None else ""

    # import the color bins. see the example color bins file. INPUT MUST BE SORTED.
    color_bins = []
    last_bin = None
    for line in args.color_bins[0].readlines():

        parts = line.split(",")
        cur_bin = {"upper": float(parts[0]), "lower": float(parts[1]), "color": parts[2].strip()}
        
        # validate that the bin polarity is correct
        if cur_bin["upper"] < cur_bin["lower"]:
            sys.exit("ERROR: upper bin cutoffs must be greater than lower bin cutoffs.\n" \
                     "offending bin: " + str(cur_bin))

        # validate that bins are sorted
        if last_bin != None and (cur_bin["upper"] > last_bin["lower"] or last_bin["upper"] < cur_bin["lower"]):
            sys.exit("ERROR: color bins must be provided in descending order and cannot overlap.\n" \
                     "high bin " + str(last_bin) + "\n" \
                     "low bin " + str(cur_bin))
    
        color_bins.append(cur_bin)
        last_bin = cur_bin
    args.color_bins[0].close()

    # get the data to be used
    data = {}
    use_node_labels_as_data = True
    if args.node_values != None:
        use_node_labels_as_data = False
        first_line = True
        for line in args.node_values[0].readlines():
#            print parts
            parts = line.strip().split(",")
        
            if first_line:
                first_line = False
                column_labels = parts
                continue
            
            if len(parts) > 1:
#                print parts[1:]
                data[parts[0]] = dict(zip(column_labels[1:],[float(p) for p in parts[1:]]))
        args.node_values[0].close()

    # load the tree
    tree = None
    while tree == None and line != "":
        line = args.tree[0].readline()
        try:
            tree = newick3.parse(StringIO(line))
        except AttributeError:
            continue
    if tree == None:
        sys.exit("Could not find a tree in: " + args.tree[0].name)
    args.tree[0].close()
    
    # now we will paint the branches
    if args.node_values != None and args.label != None: # use values from the node-values file
        
        for column_label in args.label:

            for node in tree.iternodes():
                if node.label in data:
                    this_bin = get_bin(data[node.label][column_label], color_bins)
                    if this_bin != None:
                        node.color_string = this_bin["color"]
                    else:
                        sys.exit("could not assign the value '" + str(data[node.label][column_label]) + "' to a bin")
                else:
                    node.color_string = BLACK
        
            outfile_name = "tree_painted_by_" + column_label + ".tre"
            outfile_name = output_prefix + "_" + outfile_name if len(output_prefix) > 0 else outfile_name
            outfile_path = target_dir+"/"+outfile_name
            print("writing painted tree for data column '" + column_label + "' to: " + outfile_path)
            write_tree(outfile_path, tree)
                
    else: # use values from the internal node labels themselves
        
        # check to see if we have tip values designated (since we won't use tip labels as paintable values) 
        tip_values = None
        if args.tip_values != None:
            tip_values = {}
            for line in args.tip_values[0].readlines():
                name, value = line.split(",")
                value = float(value)
                tip_values[name] = value
        
        for node in tree.iternodes():
            if len(node.children) > 0:
                value = float(node.label)
            elif tip_values != None and node.label in tip_values:
                value = tip_values[node.label]
            else:
                node.color_string = BLACK
                continue
            
            this_bin = get_bin(value, color_bins)
            if this_bin != None:
                node.color_string = this_bin["color"]
            else:
                sys.stderr.write("could not assign the value '" + node.label + "' to a bin")

        outfile_name = "tree_painted_by_node_labels.tre"
        outfile_name = output_prefix + "_" + outfile_name if len(output_prefix) > 0 else outfile_name
        outfile_path = target_dir+"/"+outfile_name
        print("writing painted tree based on node labels to: " + outfile_path)
        write_tree(outfile_name, tree)
