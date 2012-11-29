"""
	counttaxa.py
	version 0.1
	cody hinchliff
	2009.10.01
	
	edited 2011.11.6
	
	This script accesses a directory, and counts the number of sequences in any and all FASTA
	files found in that directory, binning them by genus, and records the counts into a file
	called 'counts_by_genus.csv' that may be opened in any conventional spreadsheet or
	text-editor app. This file will be overwritten if it already exists.
	
	To call it, at the command prompt just type:
	
	>python count_seqs_by_genus.py [pathtodirectory]/
	
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

print "Reading files from: " + dirpath

try:
    myfiles = os.listdir(dirpath)
except OSError:
    exit("There was a problem opening the specified directory. Are you sure it exists? Quitting.")

alltaxa = list()
seqcountbyfile = dict()

for currfname in myfiles: 

	skip = False
	try:
		myfile = file(dirpath + "/" + currfname, 'rU')
		records = [record for record in SeqIO.parse(myfile,"fasta")]
		print "Reading: " + currfname
	except IOError:
		print "Error: there was a problem opening the file: " + currfname
		skip = True
			
	if not skip:
		if len(records) <= 0:
			print "Warning: no records found in: " + currfname + ". Skipping this file."
		else:
			seqcount = seqcountbyfile[currfname] = dict()
			
			for record in records:
				currgenus = record.id.split("_")[0]
				
				try:
					seqcount[currgenus] += 1
				except KeyError:
					seqcount[currgenus] = 1
				
				try:
					alltaxa.index(currgenus)
				except ValueError:
					alltaxa.append(currgenus)


savedcounts = dict()
for thisfile in seqcountbyfile:
    if len(seqcountbyfile[thisfile]) != 0:
        savedcounts[thisfile] = seqcountbyfile[thisfile]

seqcountbyfile = savedcounts

outfname = dirpath + "/counts_by_genus.csv"
outfile = file(outfname,"w")

alltaxa.sort()
firstline = "file"

for taxon in alltaxa:
    firstline += ("," + taxon.strip())

outfile.write(firstline+"\n")

for thisfile in seqcountbyfile:
	
    currcount = seqcountbyfile[thisfile]
    print "Recording counts from: " + thisfile
	
    currline = thisfile
    for taxon in alltaxa:
		
        currline += ","
        if taxon in currcount:
            currline += str(currcount[taxon])
        else:
            currline += "0"
	
    outfile.write(currline+"\n")

