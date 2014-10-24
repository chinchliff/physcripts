#!/usr/bin/env python

import sys
import os
import time
#import MySQLdb

from Bio import SeqIO
from Bio import AlignIO
from Bio import Alphabet #import generic_dna
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# user libraries, available from bitbucket python module repo 
import dnaalignment

#print "\n\nEnter the path to the file:"
#filename = raw_input("> ")

infname = sys.argv[1]
cwd = os.getcwd() + os.sep

try:
        infpath = cwd + infname
        myfile = file(infpath, 'rU')
        mydata = SeqIO.parse(myfile,"fasta")
        print "\nReading from: %s" % infpath

except IOError:
        exit("Error: there was a problem opening the input file: " + infpath)

try:
        fnamebits = os.path.splitext(infname) 
        outfname = """%s_dupsconsensed.fasta""" % (fnamebits[0])
        outfpath = cwd + outfname
        outfile = file(cwd + outfname,"w")
        print "Writing to: %s \n" % outfpath

except IOError:
        exit("Error: there was a problem opening the output file: " + outfpath)

# Vars to accumulate metadata for seqs we've processed, to be used during/after dup checking
checkedspp = list()
tempfnames = list()
odata = list()
cdata = list(mydata)

# Variables to hold temporary values during dup checking
appending = False
dups = list()

# for every record (iterate by record index)
for curindex in range(len(cdata)-1):

    compdata = list(cdata) # make a list copy for iteration
    currecord =  compdata.pop(curindex)

    ## if the species of the current record hasn't yet been seen, then look for all other records
    ## of the same species, and append them to the duplicates list
    if currecord.id not in checkedspp:
        for comprecord in compdata:
            if currecord.id == comprecord.id:
                dups.append(comprecord)
                appending = True

    # if we found duplicates
    if appending:
        dups.append(currecord) # add the record we just checked against to the list of dups 
        print "%s: \n\t%i sequences found" % (currecord.id,len(dups)) 

        tempfname = "%s%s%f_TEMP.fasta" % (cwd, currecord.id, time.time())
        temphandle = file(tempfname, "w")
        tempfnames.append(tempfname)

        # save the records to a temp file so AlignIO can access them                
        SeqIO.write(dups, temphandle, "fasta")
        temphandle.close()

        # read all records from tempfile
        duphandle = file(tempfname, "rU")
        alignment = AlignIO.read(duphandle,"fasta")

        cons = ConAlign.consensus(alignment)
        
#        print "the consensus is %s" % cons.consensus_seq
        print "\t%i polymorphic columns consensed\n" % cons.ncols_polymorphic
        
        os.remove(tempfname)

        conrecord = SeqRecord(Seq(cons.consensus_seq,Alphabet.Gapped(Alphabet.generic_dna)),id=currecord.id, \
                    name=currecord.id+" consensus sequence",description=currecord.id+" consensus sequence")

        odata.append(conrecord)

    elif currecord.id not in checkedspp:
        odata.append(currecord)

    # reset temp vars
    dups = list()
    appending = False
    
    checkedspp.append(currecord.id) # add the species we checked to our ignore list for future checks

SeqIO.write(odata, outfile, "fasta")
