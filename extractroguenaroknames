#!/usr/bin/env python

import os
import sys

if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        print "usage extractroguenaroknames.py <roguenarokresultsfile>"
        sys.exit(0);

    roguesfile = open(sys.argv[1],"r")
    outfile = open(roguesfile.name+".names","w")
    for line in roguesfile:
        for name in line.split()[2].split(","):
            outfile.write(name + "\n")

    outfile.close()
    print "rogue names written to "+outfile.name
