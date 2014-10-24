#!/usr/bin/env python
'''filters a phylip alignment; can accept *either* a set of names to accepted or a set to be excluded, and saves the output to another file'''

if __name__ == "__main__":

    import os
    import sys

    if len(sys.argv) < 4:
        print "usage: filter_phylip.py <infile> [excluded=<excludedtaxafile> | accepted=<acceptednamesfile>] <outfile>"
        sys.exit(0)

    infile = open(sys.argv[1],"r")
    outfile = open(sys.argv[3],"w")
    
    listtype, namesfile = sys.argv[2].split("=")
    

    if listtype == "excluded" or listtype == "accepted":
        names = [n.strip() for n in open(namesfile,"r").readlines()]
    else:
        print "unrecognized type for names list; please use 'excluded' or 'accepted'"
        sys.exit(0)

    while True:
        testline = infile.readline()
        try:
            ntax, nsites = [int(n) for n in testline.split()]
            break
        except IndexError:
            continue

    saved = list()

    ntax = 0
    for name, seq in [line.split() for line in infile]:

        if listtype == "excluded":
            if name in names:
                continue

        elif listtype == "accepted":
            if name not in names:
                continue

        saved.append(name + " " + seq + "\n")
        ntax += 1

    outfile.write(str(ntax) + " " + str(nsites) + "\n")

    for line in saved:
        outfile.write(line);

    outfile.close()
