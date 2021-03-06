#!/usr/bin/env python

if __name__ == "__main__":

    import sys
    
    if len(sys.argv) < 2:
        print "usage: makesamplingmatrixfromalignment <phylipaln>"
        sys.exit(0)

    infname = sys.argv[1]
    infile = open(infname, "r")

    outfile = open(infname.rsplit(".phy",1)[0] + ".sitewisesamplingmatrix.csv","w")

    empty_chars = set(["N","-"])

    first = True
    for line in infile:
        
        if first:
            first = False
            ntax, nsites = [int(n) for n in line.strip().split(" ")]
            outfile.write("taxon," + ",".join(["site_"+str(n) for n in range(0,nsites)]) + "\n")
            continue

        taxname, seq = line.strip().split(" ")
        cols = []
        for char in seq:
            if char in empty_chars:
                cols.append("0")
            else:
                cols.append("1")

        outfile.write(",".join([taxname,] + cols) + "\n")

    outfile.close()
