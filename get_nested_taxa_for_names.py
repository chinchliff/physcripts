#!/usr/bin/env python

"""uses a phlawd sqlite3 taxonomy database to retrieve all the names associated with each taxon 
in the provided list. <targetrank> specifies the rank of the taxa to be returned."""

import sys, sqlite3

if len(sys.argv) < 5:
    print("usage: get_nested_taxa_for_names.py <targetrank> <taxonnamesfile> <gbdb> <outfile>")
    sys.exit(0)

target_rank = sys.argv[1].strip();

higher_taxon_names = [s.strip() for s in open(sys.argv[2]).readlines()]

conn = sqlite3.connect(sys.argv[3])
c = conn.cursor()

outfile = open (sys.argv[4],"w")
outfile.write("higher_taxon,child_name\n")

for name in higher_taxon_names:
    if len(name) > 0:
        print name
        results = c.execute("SELECT left_value, right_value FROM taxonomy WHERE name == ? ", (name,))
        res = c.fetchone()
        if (res != None):
            lval, rval = res
#        print(lval, rval)
        else:
            print("Name " + name + " is not in the taxonomy")
            outfile.write(",".join([name,"(none)"]) +"\n")
            continue

        results = c.execute("SELECT name FROM taxonomy WHERE left_value > ? AND " + \
            "right_value < ? AND name_class == 'scientific name' AND node_rank == ?", (lval, rval, target_rank))

        for child_name in results.fetchall():
            outfile.write(",".join([name,child_name[0]]) +"\n")

outfile.close()
