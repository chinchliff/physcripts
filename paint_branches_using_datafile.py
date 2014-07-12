#!/usr/bin/env python

import newick3, phylo3, sys, os
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

if __name__ == "__main__":

    if len(sys.argv) < 5:
        print("usage: paint_branches.py <labeledtreefile> <nodescoresfile> <colorbinsfile> <columnlabel1[,columnlabel2,...]> [target_dir=<targetdir>]>")
        sys.exit(0)

    tree_file_path = sys.argv[1] 
    node_scores_file_path = sys.argv[2]
    color_bins_file_path = sys.argv[3] 
    labels_to_use = sys.argv[4].split(",")

    target_dir = os.path.abspath(".")
    if len(sys.argv) > 5:
        parts = sys.argv[5].split("=")
        if parts[0] == "target_dir":
            target_dir = os.path.abspath(parts[1])

    # import the color bins. see the example color bins file. INPUT MUST BE SORTED.
    color_bins = []
    with open(color_bins_file_path, "r") as color_bins_file:

        last_bin = None
        for line in color_bins_file.readlines():

            parts = line.split(",")
            cur_bin = {"upper": float(parts[0]), "lower": float(parts[1]), "color": parts[2].strip()}

            # validate that bins are sorted
            if last_bin != None and (cur_bin["upper"] > last_bin["lower"] or last_bin["upper"] < cur_bin["lower"]):
                print("ERROR: color bins must be provided in descending order and cannot overlap.")
                print "high bin " + str(last_bin)
                print "low bin " + str(cur_bin)
                sys.exit(0)
        
            color_bins.append(cur_bin)
            last_bin = cur_bin

    # get the data to be used
    data = {}
    with open(node_scores_file_path) as node_scores_file:
    
        first_line = True
        for line in node_scores_file:
            parts = line.strip().split(",")
        
            if first_line:
                first_line = False
                column_labels = parts
                continue
            
            if len(parts) > 1:
                data[parts[0]] = dict(zip(column_labels[1:],[float(p) for p in parts[1:]]))

    with open(tree_file_path, "r") as tree_file:
        tree = newick3.parse(tree_file.readline())
    
        for column_label in labels_to_use:

            for node in tree.iternodes():
                if node.label in data:
                    this_bin = get_bin(data[node.label][column_label], color_bins)
                    if this_bin != None:
                        node.color_string = this_bin["color"]
                    else:
                        raise "could not assign the value '" + str(data[node.label][column_label]) + "' to a bin"
                else:
                    node.color_string = ""
        
            print("writing painted tree for data column '" + column_label + "'")
        
            outfile = open(target_dir+"/tree_painted_by_" + column_label + ".tre","w")
            outfile.write("#NEXUS\nbegin taxa\n\tdimensions ntax=" + str(len(tree.leaves())) + ";\n\ttaxlabels\n")
            for tip in tree.leaves():
                outfile.write("\t"+tip.label+"\n")
            outfile.write(";\nend;\n\nbegin trees;\n\ttree tree1 = [&R]\t")

            outfile.write(make_newick_with_color_strings(tree))

            outfile.write(";\nend;\n\n")
            outfile.write(figtree_blocks.circular_sorted_tree)
            outfile.close()


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
