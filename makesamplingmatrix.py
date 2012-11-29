#!/usr/bin/env python
"""
makesamplingmatrix.py
version 0.1
cody hinchliff
2011.3.30

This script accesses a directory, and traverses all FASTA files in it, recording the names of all taxa present
in each file. Then it creates a tab-delimited file containing a matrix where the rows represent the taxa and the
columns the FASTA files. The intended use is for a directory containing a set of FASTA files each corresponding
to a single locus, and containing homologous sequences of that locus for different taxa. The script will record
a 1 in the resulting matrix if a taxon is present in a locus file, or a 0 if not.

Key point: the script does not intelligently differentiate FASTA files from other types, and it will attempt
to parse any file in the directory. For this reason, you should remove all other files before you run the script.

To call it, at the command prompt just type:

>python makesamplingmatrix.py [pathtodirectory]/

It will create (or overwrite!) a file in the passed directory called 'sampling_matrix.txt' that may be opened
in any conventional spreadsheet or text-editor app. This file should be in the proper format for use in the
Decisivator application.

This script requires BioPython to be installed.

"""

import sys
import os
from Bio import SeqIO

cwd = os.getcwd() + os.sep

try:
    specdir = sys.argv[1]
    if specdir[0] == '/':
        dirpath = specdir
    else:
        dirpath = cwd + specdir
except IndexError:
	print "No directory specified, using current working directory."
	dirpath = ""

if dirpath[len(dirpath) - 1] != "/":
	dirpath += "/"

print "Reading files from: " + dirpath

try:
    myfiles = os.listdir(dirpath)
except OSError:
    exit("There was a problem opening the specified directory. Are you sure it exists? Quitting.")

taxa = dict()
loci = list()

for curfile in myfiles:
	try:
		seqfile = file(dirpath + curfile, 'rU')
		seqfilepathtokens = seqfile.name.split("/")
		locusname = seqfilepathtokens[len(seqfilepathtokens) - 1].split(".")[0]

		seqs = SeqIO.parse(seqfile,"fasta") # change this for other sequence file types, such as nexus, etc.

		if locusname not in ("", "sampling_matrix"):
			loci.append(locusname)
			
		for seq in seqs:
			if not seq.id in taxa:
				taxa[seq.id] = list()
		
			taxa[seq.id].append(locusname)
	
	except IOError:
		print "Error: there was a problem reading the FASTA file: " + dirpath + "/" + curfile

outfname = dirpath + "sampling_matrix.txt"
outfile = file(outfname,"w")

taxlist = taxa.keys()
taxlist.sort()
firstline = "\t"

lastlocus = loci[(len(loci) - 1)]

for locus in loci:
    firstline += ("\t" + locus.strip())

outfile.write(firstline+"\n")

for taxon in taxlist:

	line = taxon + "\t"
	for locus in loci:	

		if locus in taxa[taxon]:		
			line += "1"
		else:
			line += "0"
	
		if locus != lastlocus:
			line += "\t"
		else:
			line += "\n"
	
	outfile.write(line)

outfile.close()

exit("The file %s was created successfully." % outfname)