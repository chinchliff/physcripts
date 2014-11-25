#!/usr/bin/env python

import newick3, os, phylo3, png, sys

size_scalar = 4
breakup = 8

if __name__ == "__main__":

    if len(sys.argv) != 4:
        print "usage: make_matrix_fig.py infile.csv infile.tre outfile.png"
        sys.exit(0)

    infile = open(sys.argv[1],"r")
    infiletree = open(sys.argv[2],"r")
    tree = newick3.parse(infiletree.readline())

    order_map = {}
    count = 0
    for i in tree.leaves():
        order_map[i.label] = count
        count += 1

    print str(len(order_map)) + " tips in tree"
    infiletree.close()

    first = True
    sampling_array = ['']*len(order_map)

    # set colors
    palette=[(0xff,0xff,0xff), (0x00,0x00,0x00)]

    for line in infile:
        if first == True:
            first = False
            continue
        lineparts = line.strip().split(",")
        tipname = lineparts[0]
        if tipname not in order_map:
            continue

        # duplicate it
        cols_str = ""
        for col in lineparts[1:]:
#            for k in range(size_scalar):
            cols_str += col
        sampling_array[order_map[tipname]] = "".join(cols_str)
        #array[order[st[0]]] = "".join(st[1:])
        #array.append("".join(st[1:]))

    print len(sampling_array)
    start = 0
    end = len(sampling_array)
#    for i in range(breakup):
#        nend = end
#        start = nend
#        end = min(start + (len(array) / breakup), len(array)) + 1
    print start, end
    #size increase in blocks
    order2 = []
    for j in sampling_array[start:end]:
#        for k in range(size_scalar):
        order2.append(j)
    #array = map(lambda x: map(int, x), array)
    tarray = map(lambda x: map(int, x), order2)
    infile.close()
    outfile = open(sys.argv[3]+"."+str(i)+".png","wb")
    #w = png.Writer(len(array[0]),len(array),greyscale=True,bitdepth=1)
    

    # write the image
    w = png.Writer(len(tarray[0]),len(tarray),palette=palette,bitdepth=1)
    w.write(outfile,tarray)
    outfile.close()
