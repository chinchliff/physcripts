#!/usr/bin/env python3

import os, sys

if len(sys.argv) < 2:
    print("usage: nexus2fasta_dirty <nexusfile> <fastaoutputfilename>")
    sys.exit(0)

with open(sys.argv[1], "r") as infile:
    with open(sys.argv[2],"w") as outfile:
        recording = False
        for line in infile:
            if line.strip().upper() == "MATRIX":
                recording = True
            elif line.strip() == ";":
                recording = False
                break
            
            if recording:
                parts = [s.strip() for s in line.split()]
                if len(parts) == 2:
                    outfile.write(">" + parts[0] + "\n" + parts[1] + "\n")
