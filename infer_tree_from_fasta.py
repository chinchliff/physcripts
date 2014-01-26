#!/usr/bin/env python

import os, sys

if len(sys.argv) < 2:
    print "infer_tree_from_fasta.py <fastafile> [--noalign] [--noclean]"
    sys.exit(0)

fastafile = sys.argv[1]
alnname = fastafile.rsplit(".fasta",1)[0]

doalign = True
docleannames = True
if len(sys.argv) > 2:
    for arg in sys.argv[2:]:
        if arg == "--noalign":
	    doalign = False
        elif arg == "--nocleannames":
            docleannames = False

if docleannames:
    cleanerargs = ["clean_gb_names_simple_fasta.py", fastafile]
    os.system(" ".join(cleanerargs))
    cleanedfasta = fastafile.rsplit(".fasta",1)[0] + ".cleaned.fasta"
    fastatoalign = cleanedfasta
else:
    fastatoalign = fastafile

if doalign:
    alignedfasta = fastatoalign + ".aligned"
    mafftargs = ["mafft", fastatoalign, ">", alignedfasta]
    os.system(" ".join(mafftargs))
    fastatoconvert = alignedfasta
else:
    fastatoconvert = fastatoalign

phylipfile = alnname + ".phy"
fasta2phylipargs = ["fasta2phylip", fastatoconvert, phylipfile]

os.system(" ".join(fasta2phylipargs))

raxmlargs = ["raxmlHPC-PTHREADS-SSE3", \
                "-s", phylipfile, \
                "-n", alnname, \
                "-p", "123", \
                "-m", "GTRCAT", \
                "-T", "4"]

os.system(" ".join(raxmlargs))

figtreeargs = ["figtree", "RAxML_bestTree." + alnname]
os.system(" ".join(figtreeargs))
