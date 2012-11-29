import random, sys, os, re, copy

if __name__ == "__main__":

    if len(sys.argv) < 4:
        print "usage: make_bs_alignments.py <alignment> <n_replicates> [<partfile>] [<randseed>]"
        sys.exit(1)

    aln_filename = sys.argv[1]
    n_replicates = int(sys.argv[2])

    print aln_filename

    if len(sys.argv) > 3:
        part_filename = sys.argv[3]
    elif len(sys.argv) > 4:
        random.seed(sys.argv[4])
    else:
        random.seed()

    # read the alignment into a dict
    aln = dict()
    aln_file = open(aln_filename,"rb")
    on_first_line = True
    for line in aln_file:
        toks = line.split()
        if len(toks) > 1:
            if on_first_line == True:
                on_first_line = False
            else:
                aln[toks[0]] = toks[1]

    # extract summary stats and taxon list
    taxa = aln.keys()
    ntaxa = len(taxa)
    ncols = len(aln[taxa[0]])

    # map the partitions to a dict we will use for lookup operations
    partmap = dict()
    part_file = open(part_filename,"rb")
    for line in part_file:
        toks = [t.strip() for t in re.split(r"[,=\r]+",line)]
        ptype = toks[0]
        name = toks[1]
        bounds = [b.strip() for b in toks[2].split("-")]
        start = int(bounds[0])
        end = int(bounds[1])
        partmap[name] = (ptype, start, end)

    partnames = partmap.keys()

    # create a dict for fast column-partition lookups
    colmap = dict()
    for pname in partnames:
        start = partmap[pname][1]
        end = partmap[pname][2]
        for i in range(start, end + 1):
            colmap[i] = pname

    # progress bar parameters
    p_bar_length = 50
    p_interval = ncols / p_bar_length
    report_intervals = [p_interval * i for i in range(p_bar_length)] + [ncols,]

    # build bootstrap alignments
    count = 1
    while count <= n_replicates:
        # bs alignment will be held in nested dict: level 1 keys are part names; level 2 keys are taxon names
        new = dict(zip(partnames,map(copy.deepcopy,[dict(zip(taxa,[""]*ntaxa))]*len(partnames))))
        print "replicate " + str(count)
        first_column = True
        j = 0
        for i in range(ncols):
            if i == report_intervals[j]:
                # update the progress bar if necessary
                if i != ncols:
                    bar_complete = "|" + ("=" * ((i/p_interval) - 1)) + ">"
                    bar_incomplete = (" " * (p_bar_length - len(bar_complete))) + "|"
                    sys.stdout.write("\r" + bar_complete + bar_incomplete + " " + str((float(j)/p_bar_length) * 100) + "%")
                    sys.stdout.flush()
                else:
                    sys.stdout.write("\r" + ("=" * p_bar_length) + "|" + " complete")
                j += 1
            # get random column, look up which partition it corresponds to
            r = random.randint(0,ncols-1)
            partname = colmap[r+1]
            for t in taxa:
                if first_column:
                    # write this column to the corresponding partition in the alignment dict
                    new[partname][t] = aln[t][r]
                else:
                    new[partname][t] += aln[t][r]
            first_column = False
            i += 1

        # prepare to write output alignment and partition files
        aln_outfilename = os.path.basename(aln_filename).rsplit(".")[0] + ".bs." + str(count) + ".phy"
        part_outfilename = os.path.basename(part_filename).rsplit(".")[0] + ".bs." + str(count) + ".part"

        aln_outfile = open(aln_outfilename,"wb")
        aln_outfile.write("%s %s\n" % (ntaxa, ncols))

        # concatenate sequences across partitions, write to phylip file
        for taxon in taxa:
            seq = ""
            for pname in partnames:
                seq += new[pname][taxon]
            aln_outfile.write("%s %s\n" % (taxon, seq))
        aln_outfile.close()

        # calc start/end of partitions, write raxml partition file
        part_outfile = open(part_outfilename,"wb")
        last = 0
        for pname in partnames:
            ptype = partmap[pname][0]
            start = last + 1
            end = last + len(new[pname][taxa[0]])
            line = ("%s, %s = %s-%s\n" % (ptype, pname, start, end))
            part_outfile.write(line)
            last = end
        part_outfile.close()

        print "\n"
        count += 1

    print "done."
