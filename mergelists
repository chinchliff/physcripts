#!/usr/bin/env python

import sys

"""Merges two linebreak-delimited text lists, excluding duplicate entries"""

if __name__ == "__main__":

    if len(sys.argv) < 4:
        print "usage: mergelists <list1> <list2> <outfile>"

    fname1 = sys.argv[1]
    fname2 = sys.argv[2]
    outfname = sys.argv[3]

    file1 = open(fname1, "rb")
    file2 = open(fname2, "rb")
    outfile = open(outfname, "wb")

    a = set([line.strip() for line in file1.readlines()])
    b = set([line.strip() for line in file2.readlines()])

    allitems = a.union(b)

    for item in allitems:
	outfile.write(item + "\n")

    outfile.close()
