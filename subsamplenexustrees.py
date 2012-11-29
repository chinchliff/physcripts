#!/usr/bin/env python
#
#  subsample_trees_nex.py
#  
#  Created by Cody on 8/26/11.
#

import sys
import os

if len(sys.argv) < 4:
	exit("Usage: \n\n./sunsample_trees_nex.py [infile] [outfile] [tree sampling frequency]\n")

try:
	relpath = sys.argv[1]
	fullpath = os.getcwd() + os.sep + relpath
	datafile = open(fullpath,'rU')
except IOError:
	exit("\nThe input file '%s' could not be found.\n" % sys.argv[1])

try:
	outfile = open(sys.argv[2],'wB')
except IOError:
	exit("\nThe output file '%s' could not be opened. Please try another name.\n" % sys.argv[2])

try:
	linefreq = int(sys.argv[3])
except ValueError:
	exit("\nYou cannot use '%s' as a line sampling frequency. Please specify an integer number.\n" % sys.argv[3]) 

i = 0
nsaved = 0
for line in datafile:
	if line.split()[0] == 'tree':
#		print "found a tree"
		if i % linefreq == 0:
#			print "saving tree from line %s" % i
			outfile.write(line)
			nsaved += 1;
	else:
		outfile.write(line)
	
	i += 1

outfile.close()
exit("Saved %s trees to the file '%s'." % (nsaved,outfile.name))