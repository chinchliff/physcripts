#!/usr/bin/env python
 
## USE THIS FILE AS EXAMPLE CODE TO WORK WITH EXCEL/CSV FILES!

import os
import csv

print "\n\nEnter the path (name) of the file:"
filename = raw_input('> ')

modelfilepath = os.getcwd() + os.sep + filename
print "\nChecking for scores logfile at: " + modelfilepath

datafile = file(modelfilepath,'rU')
datafile.readlines()
obsnewlines = datafile.newlines

datafile.seek(0)
data = datafile.readlines()

print "\nAttempting to parse data file..."

cleandata = list()
dontsave = True
i = 1
for line in data:
    words = line.split()
    if len(words) != 0:
        if dontsave == False:
            if words[0] == "-ln":
                cleandata.append([i, float(words[2]) * -1])
                dontsave = True
                i=i+1
        elif words[0] == "Tree":
            dontsave = False

outfile = open(filename + "_extracted.csv","wb")
csvwriter = csv.writer(outfile, dialect='excel')

if csvwriter.writerow(["model","lnL score"]): 
	print "\nOk! Writing data to " + outfile.name

for row in cleandata:
    # for some reason this only works in the regular python compiler, not ipython
    csvwriter.writerow(row)

print "\nWrite successful. Open " + outfile.name + " to access scores.\n\n"
