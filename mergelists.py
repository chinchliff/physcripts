#!/usr/bin/env python2.6
#
#  merge_txt_lists.py
#  
#
#  Created by Cody on 4/28/11.
#  Copyright (c) 2011 The Evergreen State College. All rights reserved.
#

fname1 = "ndhf_taxa.txt"
fname2 = "rbcl_taxa.txt"

outfname = "has_ndhf_or_rbcl_taxa.txt"

file1 = open(fname1, "rb")
file2 = open(fname2, "rb")
outfile = open(outfname, "wb")

taxa = [line.strip() for line in file1.readlines()]

for line in file2:
	taxon = line.strip()
	if taxon not in taxa:
		taxa.append(taxon)

taxa.sort()

for taxon in taxa:
	outfile.write(taxon + "\n")

outfile.close()