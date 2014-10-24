#!/usr/bin/env python

if __name__ == '__main__':

    import os, sys, re

    if len(sys.argv) < 2:
        print "clean_gb_names_simple_fasta.py <fastafile>"
        sys.exit(0)

    fastafile = sys.argv[1]
    alnname = fastafile.rsplit(".fasta",1)[0]
    alignedfasta = alnname + ".cleaned.fasta"

    outfile = open(alignedfasta,"w")
    i = 0
    for line in open(fastafile, "rU").readlines():
        line = line.strip()
        if len(line) > 0:
            if line[0] == ">":
                parts = line.split("|")
                if len(parts) > 3:
                    label = parts[4].strip()
                else:
                    label = line.strip()
                newlabel = re.sub(r"[\s!@#$%^&*()\[\]{}:\",./\\<>?|]+","_",label)[:40]+"_"+str(i)
                outfile.write(">"+newlabel+"\n")
                i += 1
            else:
                outfile.write(line+"\n")
