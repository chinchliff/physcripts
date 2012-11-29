#!/Library/Frameworks/Python.framework/Versions/2.6/bin/python

from Bio import AlignIO
from Bio.Alphabet import generic_dna
import os
import sys

def convert(ifname,iformat,ofname,oformat):
	
	ihandle = open(ifname,"rU")
	ohandle = open(ofname,"wb")

	alignment = AlignIO.read(ihandle,iformat,alphabet=generic_dna)
	
	AlignIO.write(alignment,ohandle,oformat)

	print "Alignment length %i" % alignment.get_alignment_length()

############

ipath = sys.argv[1]
iformat = sys.argv[2]

if os.path.isdir(ipath):
	usedir = True

	if ipath[len(ipath) - 1] != '/':
		ipath += '/'
			
	oformat = sys.argv[3]

	ifiles = os.listdir(ipath)	
	for ifname in ifiles:
		if ifname[0] != '.':

			ofname = ifname + "." + oformat
			
			print "\nconverting " + ifname + " > " + ofname
			
			try:
				convert(ipath + ifname, iformat, ipath + ofname, oformat)
			except ValueError:
				print "skipped due to errors - mismatched input format?"
				os.remove(ipath + ofname)

else:
	ofname = sys.argv[3]
	oformat = sys.argv[4]

	convert(ipath, iformat, ofname, oformat)

