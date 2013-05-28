#!/usr/bin/env python
import sys, sqlite3, newick3, phylo3

if __name__ == "__main__":

    if len(sys.argv) < 5:
#        print "usage: testmonophylyagainsttree <phlawddb> rank=<rank> [include=\"<tx1>,<tx2>,...\"] " \
#            "[includefile=<file>] [exclude=\"<tx3>,<tx4>,...\"]"
        print "usage: testmonophylyagainsttree <phlawddb> <treefile> rank=<rank> [include=\"<tx1>,<tx2>,...\"] [includefile=<file>]"
        sys.exit(0)

    # process command line args
    dbname = sys.argv[1]
    treefname = sys.argv[2]

    # initializing parameters
    target_rank = ""
    cladenames = []
#    exclude_names = []

    for arg in sys.argv[3:]:
        argname, argval = arg.split("=")

        if argname == "rank":
            target_rank = argval.strip()

        elif argname == "include":
            cladenames = [n.strip() for n in argval.split(",")]

        elif argname == "includefile":
            includenamesfile = open(argval,"r")
            cladenames = [n.strip() for n in includenamesfile.readlines()]
            includenamesfile.close()

#        elif argname == "exclude":
#            exclude_names = [n.strip() for n in argval.split(",")]

    assert(len(cladenames) > 0)
    assert(target_rank != "")

    test_tree = newick3.parse(open(treefname,"r"))

    print "will assess monophyly for taxa of rank '" + target_rank + "'"

    con = sqlite3.connect(dbname)
    cur = con.cursor()

    included_taxa = {}
    for name in cladenames:
        cur.execute("SELECT name, left_value, right_value FROM taxonomy WHERE name_class == 'scientific name' AND name LIKE ?",(name,))
        for curname, leftval, rightval in cur.fetchall():
            print "including " + name
            included_taxa[name] = (leftval, rightval)
    
    out_prefix = "_".join(included_taxa)
    if len(out_prefix) > 40:
        out_prefix = out_prefix[0:40] + "_etc"
    out_prefix = out_prefix + "_" + target_rank
    outfile = open(out_prefix + ".monophyly.csv","w")
    outfile.write("taxon,monophyletic,offending_names\n")

    for incl_name, rl_values in included_taxa.iteritems():

        print "now searching " + incl_name

        cur.execute("SELECT name, id, left_value, right_value FROM taxonomy WHERE name_class == 'scientific name' " \
                    "AND node_rank == ? AND left_value > ? AND right_value < ?",(target_rank, rl_values[0], rl_values[1]))
        results = cur.fetchall()

        for tax_name, tax_id, tax_leftval, tax_rightval in results:
        
            cur.execute("SELECT name FROM taxonomy WHERE name_class == 'scientific name' AND left_value > ? AND right_value < ?", \
                            (tax_leftval,tax_rightval))

            taxo_mrca_names = [n[0] for n in cur.fetchall()]

            if len(taxo_mrca_names) < 1:
                print tax_name + " has no children, will be skipped"
                continue

            print tax_name

            tree_mrca = phylo3.getMRCA(test_tree, taxo_mrca_names)

            if tree_mrca == None:
                print "    could not find any descendants in tree, will be skipped"
                continue

            tree_mrca_names = [n.label for n in tree_mrca.descendants() if n.label != None]
            
            monophyly = False
            offending_names = set(tree_mrca_names) - set(taxo_mrca_names)
            if len(offending_names) == 0:
                monophyly = True

            outfile.write(",".join([tax_name, str(monophyly), " | ".join(offending_names)]) + "\n")
            outfile.flush()

    outfile.close()
