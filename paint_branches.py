#!/usr/bin/env python

import argparse, newick3, phylo3, sys, os
import figtree_blocks

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

    '''Write the tree, painted by the color bins, in the figtree format, to the specified outfile. Assumes
    the tree nodes contain color_string attributes that will be populated into the output.'''

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

    parser.add_argument("-t", "--tree", type=file, nargs=1, required=True, help="A newick tree file with labeled internal nodes. If a node-values file is supplied, then the node labels in the tree will be used to match to the rows of the node-values file. Otherwise, the node labels themselves will be interpreted as the values to be used for tree painting.")

    parser.add_argument("-c", "--color-bins", type=file, nargs=1, required=True, help="A comma-separated list containing min and max values for the bins to which node values will be assigned, and the color which each branch falling into that bin should be painted.")

    parser.add_argument("-d", "--node-values", type=file, nargs=1, help="A list containing N columns of comma-separated values corresponding to nodes in the tree, where the first column is the node label of some internal node in the tree, and the other columns are values assigned to that node. If a node-values file is supplied, the --label argument must used to specify which value column(s) in the node-values file to be used for tree painting.")

    parser.add_argument("-l", "--label", nargs="+", help="One or more column labels from the node-values file whose values will be used for tree-painting. If multiple values are specified, one output tree will be produced for each column of values, painted according to the color scheme specified by the color-bins file. If no node-values file has been designated, all --label values will be ignored.")
    
    parser.add_argument("-v", "--tip-values", type=file, nargs=1, help="A list containing 2 columns of comma-separated values, where the first item on each line is the tip label and the second is the value to be assigned to the the tip. This may be used to assign values to tip nodes in a tree whose internal nodes are assigned values with node labels.")

    parser.add_argument("-o", "--output-dir", nargs=1, help="A location to which output files should be saved. If not specified they will be saved to the current working directory.")

    args = parser.parse_args()
    
    target_dir = os.path.abspath(".")
    if args.output_dir != None:
        target_dir = os.path.abspath(args.output_dir)
    print("will save output to: " + target_dir)

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
            parts = line.strip().split(",")
        
            if first_line:
                first_line = False
                column_labels = parts
                continue
            
            if len(parts) > 1:
                data[parts[0]] = dict(zip(column_labels[1:],[float(p) for p in parts[1:]]))
        args.node_values[0].close()

    # load the tree
    tree = None
    while tree == None and line != "":
        line = args.tree[0].readline()
        try:
            tree = newick3.parse(line)
        except AttributeError:
            continue
    if tree == None:
        sys.exit("Could not find a tree in: " + args.tree[0].name)
    args.tree[0].close()
    
    if args.node_values != None and args.label != None: # use values from the node-values file
        
        for column_label in args.label:

            for node in tree.iternodes():
                if node.label in data:
                    this_bin = get_bin(data[node.label][column_label], color_bins)
                    if this_bin != None:
                        node.color_string = this_bin["color"]
                    else:
                        raise "could not assign the value '" + str(data[node.label][column_label]) + "' to a bin"
                else:
                    node.color_string = "000000" # black
        
            print("writing painted tree for data column '" + column_label + "'")            
            outfile_name = target_dir+"/tree_painted_by_" + column_label + ".tre"
            write_tree(outfile_name, tree)
                
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
#                print "using internal node label as value"
                value = float(node.label)
            elif tip_values != None and node.label in tip_values:
#                print "looking for " + node.label + " in tip labels"
                value = tip_values[node.label]
            else:
#                print "could not find a value for node " + node.label
                node.color_string = "000000" # black
                continue
            
            this_bin = get_bin(value, color_bins)
            if this_bin != None:
                node.color_string = this_bin["color"]
            else:
                sys.stderr.write("could not assign the value '" + node.label + "' to a bin")

        print("writing painted tree based on node labels")
        outfile_name = target_dir+"/tree_painted_by_node_labels.tre"
        write_tree(outfile_name, tree)

exit() ############################################## below here lies nothing (except cannibalized code)

def make_newick_painted_dec_single_locus(node, locusname, dec_color, ind_color):

    if not node.istip:

        substrings = []
        for child in node.children:
            substrings.append(make_newick_painted_dec_single_locus(child, locusname, dec_color, ind_color))

        if locusname in node.dec_loci:
            c = dec_color
        else:
            c = ind_color

        return "(" + ",".join(substrings) + ")[&!color=#" + c + "]:" + str(node.length)

    else:

        if locusname in node.dec_loci:
            c = dec_color
        else:
            c = ind_color

        return node.label + "[&!color=#" + c + "]:" + str(node.length)

def make_newick_painted_by_dec_loci(node, max_count):

    if not node.istip:

        substrings = []
        for child in node.children:
            substrings.append(make_newick_painted_by_dec_loci(child, maxcount))

        n_dec_loci = len(node.dec_loci)
        if n_dec_loci < max_count:
            c = color[len(node.dec_loci)]
        else:
            c = color[max_count]

        return "(" + ",".join(substrings) + ")[&!color=#" + c + "]:" + str(node.length)

    else:

        n_dec_loci = len(node.dec_loci)
        if n_dec_loci < max_count:
            c = color[len(node.dec_loci)]
        else:
            c = color[max_count]

        return node.label + "[&!color=#" + c + "]:" + str(node.length)

if len(sys.argv) < 4:
    print "usage: decisivate_branches <treefile> <autophy_metadata_file> <colorsfile>"
    sys.exit(0)

#treefname = "test.tre"
#samplingfname = "test.loci.txt"

# get infile names
treefname = sys.argv[1]
samplingfname = sys.argv[2]
colorsfname = sys.argv[3]

# open tree
treefile = open(treefname,"r")
tree = newick3.parse(treefile.readline())

# get locus sampling info
samplingfile = open(samplingfname,"r")
loci = {}
curlocus = ""
first = True
for line in samplingfile:
    if first:
        first = False
        print "assuming first line of locus file is column labels, it will be skipped"
        continue

    if len(line.strip()) == 0:
        continue

    parts = line.split(",")
    tipname = parts[0]
    locusname = parts[1]
    if locusname not in loci.keys():
        loci[locusname] = []
    
    loci[locusname].append(tipname)        

# get color info
colordata = open(colorsfname,"r")
color = {}
first = True
maxcount = 0
for line in colordata:
    if first:
        first = False
        print "assuming first line of colors file is column labels, it will be skipped"
        continue

    toks = line.split(",")
    key = int(toks[0])
    c = toks[5].strip()
    color[key] = c
    maxcount = key

print """finding the loci that could inform the parent branch of each node"""
dec_loci_for_branches = {}
for node in tree.iternodes():

    node.dec_loci = []

    if len(node.children) > 2:
        print "node " + node + " " + node.label + " has more than 2 children!"
        continue

    if not node.istip and node.parent != None:
        left_child = node.children[0]
        right_child = node.children[1]

        for n in node.parent.children:
            if n != node:
                sister = n

        assert(sister)

    # this is brute force method. there should be a more efficient way
    # using recursion from tips to keep track of branches whose child branches
    # we've already validated

    for locus_name, locus_exemplars in loci.iteritems():

        if node.parent == None:
            continue

        if node.istip:
            if node.label in locus_exemplars:
                node.dec_loci.append(locus_name)
            continue

        # have to have one exemplar of each child and one child of sister
        left_ok = False
        right_ok = False
        sister_ok = False

        for lchild in left_child.leaves():        
            if lchild.label in locus_exemplars:
                left_ok = True
                break

        if left_ok:
            for rchild in right_child.leaves():
                if rchild.label in locus_exemplars:
                    right_ok = True
                    break

        if right_ok:
            for schild in sister.leaves():
                if schild.label in locus_exemplars:
                    sister_ok = True
                    break

        if left_ok and right_ok and sister_ok:
            node.dec_loci.append(locus_name)

# circular, sorted tree
figtree_block = "begin figtree;\n\tset appearance.backgroundColorAttribute=\"Default\";\n\tset appearance.backgroundColour=#-1;\n\tset appearance.branchColorAttribute=\"User selection\";\n\tset appearance.branchLineWidth=2.0;\n\tset appearance.branchMinLineWidth=0.0;\n\tset appearance.branchWidthAttribute=\"Fixed\";\n\tset appearance.foregroundColour=#-16777216;\n\tset appearance.selectionColour=#-2144520576;\n\tset branchLabels.colorAttribute=\"User selection\";\n\tset branchLabels.displayAttribute=\"Branch times\";\n\tset branchLabels.fontName=\"sansserif\";\n\tset branchLabels.fontSize=8;\n\tset branchLabels.fontStyle=0;\n\tset branchLabels.isShown=false;\n\tset branchLabels.significantDigits=4;\n\tset layout.expansion=0;\n\tset layout.layoutType=\"POLAR\";\n\tset layout.zoom=53;\n\tset legend.attribute=null;\n\tset legend.fontSize=10.0;\n\tset legend.isShown=false;\n\tset legend.significantDigits=4;\n\tset nodeBars.barWidth=4.0;\n\tset nodeBars.displayAttribute=null;\n\tset nodeBars.isShown=false;\n\tset nodeLabels.colorAttribute=\"User selection\";\n\tset nodeLabels.displayAttribute=\"Node ages\";\n\tset nodeLabels.fontName=\"sansserif\";\n\tset nodeLabels.fontSize=8;\n\tset nodeLabels.fontStyle=0;\n\tset nodeLabels.isShown=false;\n\tset nodeLabels.significantDigits=4;\n\tset nodeShape.colourAttribute=\"User selection\";\n\tset nodeShape.isShown=false;\n\tset nodeShape.minSize=10.0;\n\tset nodeShape.scaleType=Width;\n\tset nodeShape.shapeType=Circle;\n\tset nodeShape.size=4.0;\n\tset nodeShape.sizeAttribute=\"Fixed\";\n\tset polarLayout.alignTipLabels=false;\n\tset polarLayout.angularRange=0;\n\tset polarLayout.rootAngle=0;\n\tset polarLayout.rootLength=100;\n\tset polarLayout.showRoot=true;\n\tset radialLayout.spread=0.0;\n\tset rectilinearLayout.alignTipLabels=false;\n\tset rectilinearLayout.curvature=0;\n\tset rectilinearLayout.rootLength=100;\n\tset scale.offsetAge=0.0;\n\tset scale.rootAge=1.0;\n\tset scale.scaleFactor=1.0;\n\tset scale.scaleRoot=false;\n\tset scaleAxis.automaticScale=true;\n\tset scaleAxis.fontSize=8.0;\n\tset scaleAxis.isShown=false;\n\tset scaleAxis.lineWidth=1.0;\n\tset scaleAxis.majorTicks=1.0;\n\tset scaleAxis.origin=0.0;\n\tset scaleAxis.reverseAxis=false;\n\tset scaleAxis.showGrid=true;\n\tset scaleBar.automaticScale=true;\n\tset scaleBar.fontSize=10.0;\n\tset scaleBar.isShown=true;\n\tset scaleBar.lineWidth=1.0;\n\tset scaleBar.scaleRange=0.0;\n\tset tipLabels.colorAttribute=\"User selection\";\n\tset tipLabels.displayAttribute=\"Names\";\n\tset tipLabels.fontName=\"sansserif\";\n\tset tipLabels.fontSize=8;\n\tset tipLabels.fontStyle=0;\n\tset tipLabels.isShown=true;\n\tset tipLabels.significantDigits=4;\n\tset trees.order=true;\n\tset trees.orderType=\"increasing\";\n\tset trees.rooting=false;\n\tset trees.rootingType=\"User Selection\";\n\tset trees.transform=false;\n\tset trees.transformType=\"cladogram\";\nend;"

# standard
#figtree_block = "begin figtree;\n\tset appearance.backgroundColorAttribute=\"Default\";\n\tset appearance.backgroundColour=#-1;\n\tset appearance.branchColorAttribute=\"User selection\";\n\tset appearance.branchLineWidth=1.0;\n\tset appearance.branchMinLineWidth=0.0;\n\tset appearance.branchWidthAttribute=\"Fixed\";\n\tset appearance.foregroundColour=#-16777216;\n\tset appearance.selectionColour=#-2144520576;\n\tset branchLabels.colorAttribute=\"User selection\";\n\tset branchLabels.displayAttribute=\"Branch times\";\n\tset branchLabels.fontName=\"sansserif\";\n\tset branchLabels.fontSize=8;\n\tset branchLabels.fontStyle=0;\n\tset branchLabels.isShown=false;\n\tset branchLabels.significantDigits=4;\n\tset layout.expansion=0;\n\tset layout.layoutType=\"RECTILINEAR\";\n\tset layout.zoom=0;\n\tset legend.attribute=null;\n\tset legend.fontSize=10.0;\n\tset legend.isShown=false;\n\tset legend.significantDigits=4;\n\tset nodeBars.barWidth=4.0;\n\tset nodeBars.displayAttribute=null;\n\tset nodeBars.isShown=false;\n\tset nodeLabels.colorAttribute=\"User selection\";\n\tset nodeLabels.displayAttribute=\"Node ages\";\n\tset nodeLabels.fontName=\"sansserif\";\n\tset nodeLabels.fontSize=8;\n\tset nodeLabels.fontStyle=0;\n\tset nodeLabels.isShown=false;\n\tset nodeLabels.significantDigits=4;\n\tset nodeShape.colourAttribute=\"User selection\";\n\tset nodeShape.isShown=false;\n\tset nodeShape.minSize=10.0;\n\tset nodeShape.scaleType=Width;\n\tset nodeShape.shapeType=Circle;\n\tset nodeShape.size=4.0;\n\tset nodeShape.sizeAttribute=\"Fixed\";\n\tset polarLayout.alignTipLabels=false;\n\tset polarLayout.angularRange=0;\n\tset polarLayout.rootAngle=0;\n\tset polarLayout.rootLength=100;\n\tset polarLayout.showRoot=true;\n\tset radialLayout.spread=0.0;\n\tset rectilinearLayout.alignTipLabels=false;\n\tset rectilinearLayout.curvature=0;\n\tset rectilinearLayout.rootLength=100;\n\tset scale.offsetAge=0.0;\n\tset scale.rootAge=1.0;\n\tset scale.scaleFactor=1.0;\n\tset scale.scaleRoot=false;\n\tset scaleAxis.automaticScale=true;\n\tset scaleAxis.fontSize=8.0;\n\tset scaleAxis.isShown=false;\n\tset scaleAxis.lineWidth=1.0;\n\tset scaleAxis.majorTicks=1.0;\n\tset scaleAxis.origin=0.0;\n\tset scaleAxis.reverseAxis=false;\n\tset scaleAxis.showGrid=true;\n\tset scaleBar.automaticScale=true;\n\tset scaleBar.fontSize=10.0;\n\tset scaleBar.isShown=true;\n\tset scaleBar.lineWidth=1.0;\n\tset scaleBar.scaleRange=0.0;\n\tset tipLabels.colorAttribute=\"User selection\";\n\tset tipLabels.displayAttribute=\"Names\";\n\tset tipLabels.fontName=\"sansserif\";\n\tset tipLabels.fontSize=8;\n\tset tipLabels.fontStyle=0;\n\tset tipLabels.isShown=true;\n\tset tipLabels.significantDigits=4;\n\tset trees.order=false;\n\tset trees.orderType=\"increasing\";\n\tset trees.rooting=false;\n\tset trees.rootingType=\"User Selection\";\n\tset trees.transform=false;\n\tset trees.transformType=\"cladogram\";\nend;"


print """writing painted tree for all loci"""
outfile = open("painted_by_dec_loci.tre","w")
outfile.write("#NEXUS\nbegin taxa\n\tdimensions ntax=" + str(len(tree.leaves())) + ";\n\ttaxlabels\n")
for tip in tree.leaves():
    outfile.write("\t"+tip.label+"\n")
outfile.write(";\nend;\n\nbegin trees;\n\ttree tree1 = [&R]\t")

outfile.write(make_newick_painted_by_dec_loci(tree, maxcount))
outfile.write(";\nend;\n\n")
outfile.write(figtree_block)

outfile.close()

counts_outfile = open("decisive_branches_per_locus.csv","w")
counts_outfile.write("locus,prop_brs_decisive,num_brs_decisive\n")

outdir = "individual_locus_decisiveness_trees/"
try:
    os.mkdir(outdir)
except OSError:
    pass

nbrs = 0
countall = True
for locusname in loci.keys():
    print "writing painted tree for " + locusname

    outfile = open(outdir+"painted_by_dec_" + locusname + ".tre","w")
    outfile.write("#NEXUS\nbegin taxa\n\tdimensions ntax=" + str(len(tree.leaves())) + ";\n\ttaxlabels\n")
    for tip in tree.leaves():
        outfile.write("\t"+tip.label+"\n")
    outfile.write(";\nend;\n\nbegin trees;\n\ttree tree1 = [&R]\t")

    outfile.write(make_newick_painted_dec_single_locus(tree, locusname, color[20], color[0]))
    outfile.write(";\nend;\n\n")
    outfile.write(figtree_block)

    outfile.close()

    if nbrs > 0:
        countall = False

    nbrs_decisive = 0
    for n in tree.iternodes():
        if locusname in n.dec_loci:
            nbrs_decisive += 1

        if countall:
            nbrs += 1
    
    counts_outfile.write(",".join([locusname, str(float(nbrs_decisive)/nbrs), str(nbrs_decisive)]) + "\n")

counts_outfile.close()

branch_counts = {}
for n in tree.iternodes():
    ndecloci = len(n.dec_loci)
    if ndecloci == 0:
        continue
    
    if ndecloci not in branch_counts.keys():
        branch_counts[ndecloci] = 0

    branch_counts[ndecloci] += 1

branch_counts_outfile = open("branch_counts_by_n_dec_loci.csv","w")
branch_counts_outfile.write("nloci_decisive,nbranches\n")
for n, c in branch_counts.iteritems():
    branch_counts_outfile.write(str(n) + "," + str(c) + "\n")
branch_counts_outfile.close()
