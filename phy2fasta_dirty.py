import sys
import os
from Bio import SeqIO

cwd = os.getcwd() + os.sep

try:
    fpath = sys.argv[1]
except IndexError:
    exit("The specified file could not be opened.")

print "Reading: " + fpath

try:
	ifile = open(fpath,"rU")
except OSError:
    exit("There was a problem opening the specified directory. Are you sure it exists? Quitting.")

ofile = open(sys.argv[1] + ".fasta","w")


lines = ifile.readlines()

for line in lines:
    bits = line.split()
    ofile.write(">" + bits[0] + "\n")
    ofile.write(bits[1] + "\n")

ofile.close()
