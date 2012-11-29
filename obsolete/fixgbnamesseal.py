#/usr/bin/env python

"""
Use this script to strip all the extra crap that genbank put into
taxon names. The script reads .seal files and finds the lines containing
taxon names. When a taxon name line is encountered, the script builds
a new name and replaces the line.

All output is saved to a new file, the original file is not touched.

Note, right now this will cut out all info except Genus and species...
This is a problem when we have subspecies. Need to think about that."""

import sys
import os

print "\n\nEnter the path (name) of the file:"
filename = raw_input('> ')

filepath = os.getcwd() + os.sep + filename
print "\nChecking for file at: " + filepath

myfile = file(filepath,'rU')
mylines = myfile.readlines()

outfile = file(filename + "_fixed.seal","w")

fixnextname = False;

for line in mylines:

#    print "Fix next name = "
#    print fixnextname
    if fixnextname == True:

        namestringpos = line.find("Name=")

        if namestringpos != -1:
            fixnextname = False

            thisname = line[namestringpos+5:len(line)]
#            print "this is the one called:"
#            print thisname

            gbidentifierpos = thisname.find("|gb|")

            if gbidentifierpos != -1:
                gbnumstart = gbidentifierpos + 4
                gbnumstop = thisname.find(".")
                gbnum = thisname[gbnumstart:gbnumstop]
#                print gbnum

                words = line.split()
                thisgenus = words[1]
                thisspecies = words[2]

                newname = thisgenus[0] + "_" + thisspecies + "_" + gbnum

#                print newname

                outfile.write("Name=\"" + newname + "\";\n")

            else:
                outfile.write(line)

        else:
            outfile.write(line)

    else:

#        print line
        if line.find("ID='PSeq'") != -1:
#            print "found it"
            fixnextname = True

        outfile.write(line)
