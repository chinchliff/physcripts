#!/usr/bin/env python

if __name__ == "__main__":

    import sys
    import os
    import re
    
    if len(sys.argv) < 4:
        print "usage: make_constraint_tree.py <starttreefile> <tipsetsfile> <outfile>"
        sys.exit(0)

    starttreepath = sys.argv[1]
    tipsetspath = sys.argv[2]
    outfilepath = sys.argv[3]

    starttreefile = open(starttreepath,"r")
    tipsetsfile = open(tipsetspath,"r")

    # read in the tip sets from the file
    tipsets = dict()
    name = ""
    for line in tipsetsfile:
        if line[0] != "\t" and line[0] != " ":
            curname = line.strip()
            tipsets[curname] = list()
        else:
            tipsets[curname].append(line.strip())

    startstring = ""
    while startstring == "":
        startstring = starttreefile.readline().strip()

    possible_names = re.findall(r"[^,:()]+",startstring)
    newstr = startstring
    for name in possible_names:
        if name in tipsets.keys():

            constraint_tip_string = ",".join(tipsets[name])
            if len(tipsets[name]) > 1:
                constraint_tip_string = "(" + constraint_tip_string + ")"

            m = re.search(r"(?<=[(,])"+name+r"(?=[,:)])",newstr)
            newstr = newstr[:m.start()] + constraint_tip_string + newstr[m.end():]

    outfile = open(outfilepath,"w")
    outfile.write(newstr)
    outfile.close()
