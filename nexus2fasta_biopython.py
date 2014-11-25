#!/usr/bin/env python

import sys
from Bio import AlignIO

"""Uses the BioPython parser to convert a nexus file to a fasta"""

if __name__ == "__main__":

    if len(sys.argv) < 3:
        print "usage: nex2fasta <nexusfile> <outfile>"

    ifile = sys.argv[1]
    iformat = "nexus"

    ofile = sys.argv[2]
    oformat = "fasta"

    ihandle = open(ifile,"r")
    ohandle = open(ofile,"w")

    alignment = AlignIO.parse(ihandle,iformat)
    AlignIO.write(alignment,ohandle,oformat)
