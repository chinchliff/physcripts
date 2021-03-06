#!/usr/bin/env python

if __name__ == '__main__':

    import newick3, phylo3, sys

    #if len(sys.argv) < 4:
    if len(sys.argv) < 3:
        print "usage: prunetips <namestoprunefile> <treefile>" #<outfile>"
        sys.exit()

    badnamesfname = sys.argv[1]
    badnamesfile = open(badnamesfname, "r")
    badnames = [name.strip() for name in badnamesfile.readlines()]

    #print badnames

    treefname = sys.argv[2]
    treefile = open(treefname, "r")
    #outfname = sys.argv[3]
    #outfile = open(outfname,"w")

    logfile = open("prunetips.log","w")

    for line in treefile:

        tree = newick3.parse(line)
    #    print "in tree: " + newick3.to_string(tree) + ";"

        # prune the bad tips
        for tip in tree.leaves():

            name_ok = False

            if tip.label in badnames:
                leftlabel = tip.label
    #            print "removing " + leftlabel
    #            logfile.write("removing " + leftlabel)
                tip.prune(logfile)
            else:
                name_ok = True

        # remove any leftover empty tip nodes    
        if not name_ok:
            for n in tree.descendants():

                nc = n
                while (not nc.istip) and len(nc.children) == 0:
     #               print "pruning an empty tip"
                    logfile.write("pruning an empty tip")
                    np = nc.parent
                    nc.prune()

                    if np:
                        # prepare for next step back
                        nc = np

                    # if we hit the root of the tree, move on to next descendant
                    else:
                        break

        print newick3.to_string(tree) + ";"

        # write the pruned tree
    #    outfile.write(newick3.to_string(tree) + ";\n")
    #    outfile.flush()

    #outfile.close()
