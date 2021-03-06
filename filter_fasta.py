#!/usr/bin/env python
"""Input files are expected to be in fasta format. The script will traverse all files
in the input dir, so the input dir should contain only fasta files. The taxon list
should be a line-delimited text file containing the names of tips as they
correspond to those in the fasta alignments."""

import sys
import os

def filter(ifname, taxa_approved):  #, ofname, taxa_approved):

	ofname = ifname + "." + "filtered"
	print "\nfiltering " + ifname + " > " + ifname + ".filtered"

	tempofname = ofname + ".temp"
	infile = open(ifname, "rb")
	outfile = open(tempofname, "wb")
	
	ntax = 0
	OnValidTaxon = False

	firstline = True
	for line in infile:
		if line[0] == '>':
			cur_taxon = line.split()[0].strip('_')
			if cur_taxon.strip('>') in taxa_approved:
				outfile.write(cur_taxon + "\n")
				ntax += 1
				OnValidTaxon = True
			else:
				OnValidTaxon = False
		elif OnValidTaxon:
			outfile.write(line)			
	
	outfile.close()
		
	if ntax == 0:
		os.remove(tempofname)
		print "No acceptable taxa found. No output saved."
	else:
		os.rename(tempofname,ofname)
		print "Retained " + str(ntaxfound) + " taxa."



if __name__ == "__main__":

	if len(sys.argv) < 3:
		print "usage: filterfasta <fastafile> <acceptednamesfile>"
		sys.exit(0)

	ipath = sys.argv[1]

	taxlistpath = sys.argv[2]
	taxlist = open(taxlistpath, "rU")
	global taxa_approved
	taxa_approved = [taxon.strip() for taxon in taxlist.readlines()]

	if os.path.isdir(ipath):
		usedir = True

		if ipath[len(ipath) - 1] != '/':
			ipath += '/'

		ifiles = os.listdir(ipath)
		for ifname in ifiles:
			if ifname[0] != '.':
				filter(ipath + ifname, taxa_approved) # ipath + ofname, taxa_approved) #ifname)
	else:
		filter(ipath, taxa_approved)
#		print "single case not yet implemented"
