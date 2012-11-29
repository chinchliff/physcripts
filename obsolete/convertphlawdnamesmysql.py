#!/usr/bin/env python

"""
convertnames.py
version 0.1
cody hinchliff
2009.10.01

This script will convert the dbids in PHLAWD output FASTA files to the appropriate species names. 
It requires the file 'convertnames.cfg' to be present in the same directory. This file should
specify the name, userid, and password associated with the database from which PHLAWD is drawing
sequences.

To call it, at the command prompt just type

>python convertnames.py [pathtoFASTAfile]

The script will create a file in the same directory as the input file, called [inputfile]_withnames.FASTA

"""

import sys
import os
import cfg_file
import MySQLdb
from Bio import SeqIO

#print "\n\nEnter the path to the file:"
#filename = raw_input("> ")

infname = sys.argv[1]
cwd = os.getcwd() + os.sep

try:
        infpath = cwd + infname
        myfile = file(infpath, 'rU')
        mydata = SeqIO.parse(myfile,"fasta")

        print "Reading from: " + infpath
except IOError:
        exit("Error: there was a problem opening the input file: " + infpath)

try:
        fnamebits = os.path.splitext(infname) 
        outfname = """%s_withnames.fasta""" % (fnamebits[0])
        outfpath = cwd + outfname
        outfile = file(cwd + outfname,"w")

        print "Writing to: " + outfpath
except IOError:
        exit("Error: there was a problem opening the output file: " + outfpath)


config = cfg_file.reader("convertnames.cfg")

convertdbids = False
trimnames = False
deleteindetspp = False
removecfs = False
deletecfspp = False

for p in config.parameters:

        if p.name == "convertdbids" and p.value.lower() == "true":
                convertdbids = True
        elif p.name == "db":
                thedb = p.value
        elif p.name == "user":
                theuser = p.value
        elif p.name == "password":
                thepasswd = p.value

        elif p.name == "trimnames" and p.value.lower() == "true":
                trimnames = True
        elif p.name == "removecfs" and p.value.lower() == "true":
                removecfs = True
        elif p.name == "deleteindetspp" and p.value.lower() == "true":
                deleteindetspp = True
        elif p.name == "deletecfspp" and p.value.lower() == "true":
                deletecfspp = True


if convertdbids:
        try:
                thedb
                theuser
                thepasswd

        except NameError:
                exit("There seems to be a problem with the config file. You need to specify the values:\n\n" \
                     "db = [dbname]\nuser = [dbuser]\npassword = [password]")
                
        db=MySQLdb
        db=MySQLdb.connect(user=theuser,passwd=thepasswd,db=thedb)
        c = db.cursor()

savedrecords = list()
for record in mydata:

        rename = True
        save = True

        if convertdbids:
                dbid = record.id.strip().split()[0]
                c.execute("""SELECT taxon_name.name FROM taxon_name WHERE taxon_name.taxon_id = '"""+dbid+"""' AND name_class = 'scientific name'""")
                
                try:
                        newname = c.fetchall()[0][0]
                except IndexError:
                        rename = False
                        print "Could not locate a corresponding record for taxon id: " + dbid
        else:
                newname = record.id

        if rename:
                badchars = "!@#$%^&*()~=+{}|[]:;'<>?,./\\\n\t\r\" "
                for char in badchars:
                        newname = newname.replace(char,"_")

                while len(newname.split("__")) > 1:
                        newname = newname.replace("__","_")

                if trimnames:
                        namebits = newname.strip("_").split("_")

                        try:
                                if namebits[1] == "cf" or namebits[1] == "aff":
                                        if deletecfspp:
                                                save = False
                                        elif removecfs:
                                                namebits[1] = namebits[2]
                                        
                                if namebits[1] == "sp" and deleteindetspp:
                                        save = False

                                newname = namebits[0] + "_" + namebits[1]

                        except IndexError:
                                pass

                record.id = newname

        if save:
                savedrecords.append(record)

SeqIO.write(savedrecords,outfile,"fasta")
outfile.close()


