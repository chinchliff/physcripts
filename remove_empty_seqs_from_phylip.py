#!/usr/bin/env python

if __name__ == '__main__':

    import sys

    if len(sys.argv) < 2:
        print "usage: removeemptyseqsfromphylip <phylipfile>"
        exit(0)

    infile = open(sys.argv[1],"r")

    savedseqs = {}
    nsites = None

    first = True
    i = 0
    n = 0
    for line in infile:
        if first:
            first = False
        continue

        i += 1
        if i % 500 == 0:
            print i

        try:
            name, seq = [q.strip() for q in line.split(" ")]
        except IndexError:
            continue

        empty = True
        for char in seq:
            if char != "-" and char != "N":
            empty = False
            break
    
        if not empty:
            savedseqs[name] = seq
            n += 1

            if nsites == None:
                nsites = len(seq)

    outfname = sys.argv[1].rsplit(".phy",1)[0] + ".noemptyseqs.phy"
    outfile = open(outfname,"w")
    outfile.write(str(n) + " " + str(nsites) + "\n")

    for name, seq in savedseqs.iteritems():
        outfile.write(" ".join([name,seq]) + "\n")
    outfile.close()
