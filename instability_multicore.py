#!/usr/bin/env python
#
#  Created by Cody Hinchliff on 10/13/11.

"""This script will calculate I^s scores, as described in (Hinchliff, C. E. and E. H. Roalson. 2012. Using supermatrices for phylogenetic inquiry: an example using the sedges. Systematic Biology). It requires a set of trees sharing a common set of tips, to be input as a newick file (though any format readable by dendropy should be trivial to use, just change the format in the appropriate line). It outputs a comma-delimited table containing the raw instability scores (the numerator from the right side of the equation in the referenced paper), as well as the scaled I^s scores. Taxa that move more have higher scores.
	
	This script requires the dendropy Python module, which itself has some other dependencies. It has been tested on Mac OS v. 10.6 and 10.7, Debian Squeeze, and Ubuntu Linux 12.04.
	
	This script will generate a directory for temporary files. These can become large for large trees and may be safely deleted once the script has finished executing."""

import os
import sys
import dendropy as dp
import cPickle as pickle
from copy import copy
from multiprocessing import Process, Array, Queue
from ctypes import c_char_p
from Queue import Empty

def calc_I_partial(dm_filename_array,taxon_labels_array,freq,offset,q):

	I_partial = dict()

	dmfiles = list(dm_filename_array)
	taxa = list(taxon_labels_array)
	t = len(taxa)

	for taxlabel in taxon_labels_array:
		I_partial[taxlabel] = 0

	# for every distance matrix x
	for i, file_x in enumerate(dmfiles):

		# only do calcs for every freq'th tree
		if (i + offset) % freq == 0 and i < len(dmfiles) - 1:
			print "Thread %s: currently comparing dists against %s" % (offset + 1, file_x)
			print "\t" + str((len(dmfiles)-(i+1))) + " matrix comparisons to be made against this matrix"

			dm_x = pickle.load(open(file_x,"rb"))
			
			# sum differences in distance for taxa t1, t2 between matrix x and all other matrices y
			for file_y in dmfiles[i+1:]:
				# read matrix y from file
				dm_y = pickle.load(open(file_y,"rb"))

				# for every taxon t1
				for j, t1 in enumerate(taxa):

					# compare t1 to to every other taxon t2
					for t2 in taxa[0:j] + taxa[j+1:t]:			

						# find distances for t1, t2 from matrices x and y
						D_ijx = dm_x[t1][t2]
						D_ijy = dm_y[t1][t2]
						
						# add distance difference to sum
						I_partial[t1] += abs(D_ijx - D_ijy)
						
	q.put(I_partial)

if len(sys.argv) < 3:
	print "usage: ./instability_multicore.py <treefile> <n_threads> [logfilename]"
	sys.exit(0)

# attempt to open the treefile
try:
	infname = sys.argv[1]
	treefile = open(infname,"rb")
except IOError:
	exit("The file '%s' could not be opened." % infname)

# second argument defines number of threads to use
try:
	nthreads = int(sys.argv[2])
except ValueError:
	exit("The number of threads must be an integer")

# third argument provides an optional logfile name 
try:
	loglabel = sys.argv[3]
except IndexError:
	loglabel = sys.argv[1]

# attempt to open logfile
try:
	logfname = loglabel+".instability.csv"
	logfile = open(logfname,"wb")
except IOError:
	"You must specify a different name for the log file."

# there is some question of whether it makes sense to use brlens
# for now we just set them all equal to one
make_edge_lengths_uniform = True

# define a taxon set to be shared by all trees
taxa = dp.TaxonSet()

# create a directory to hold temp files
try:
	wdir = os.path.dirname(os.path.abspath(sys.argv[1]))
	temppath = sys.argv[1]+"_temp"
	os.mkdir(temppath)
except OSError:
	pass

# read trees from file
trees = dp.TreeList(stream=treefile, schema='newick', taxon_set=taxa)
ntrees = len(trees)
print "\nread %s trees from the input file." % ntrees

if ntrees < 100:
	print "\n*** note: at least 100 trees is recommended for reliable precision. ***"

rfraction = 5
nrandtrees = int(ntrees / rfraction)
print "\none out of every %s trees (total %s) will have its taxa randomly assigned, and these\n" \
	  "randomized trees will be used to estimate average expected pairwise taxon distance." % (rfraction, nrandtrees)

#if useHD:
fnamelist = list()
randfnamelist = list()

j = 0
# calc distance matrices for all trees
for i, curtree in enumerate(trees):

	if make_edge_lengths_uniform:
		for n in curtree.get_node_set():
			n.edge_length = 1
	
	# get patristic distance matrix for each tree
	dmatrix = dp.treecalc.PatristicDistanceMatrix(curtree)

	# save into a dict for less overhead and easier pickling
	dmdict = dict()
	for t1 in taxa:
		dmdict[t1.label] = dict()
		for t2 in taxa:
			dmdict[t1.label][t2.label] = dmatrix(t1,t2)

	tempfile = open(temppath + "/%s.temp" % i, "wb")
	pickle.dump(dmdict,tempfile)
	fnamelist.append(os.path.abspath(tempfile.name))
	tempfile.close()
	
	# randomize each rfraction-th tree to use for estimating max theoretical distances
	if i % rfraction == 0:
		randtree = dp.Tree(curtree)
		randtree.randomly_assign_taxa(create_required_taxa=False)
		randdmatrix = dp.treecalc.PatristicDistanceMatrix(randtree)

		# save into a dict for less overhead and easier pickling
		randdmdict = dict()
		for t1 in taxa:
			randdmdict[t1.label] = dict()
			for t2 in taxa:
				randdmdict[t1.label][t2.label] = dmatrix(t1,t2)
		
		randtempfile = open(temppath + "/random_%s.temp" % j, "wb")
		pickle.dump(randdmdict,randtempfile)
		randfnamelist.append(os.path.abspath(randtempfile.name))
		randtempfile.close()

		j += 1
	
	print "stored distance matrix for tree %s" % i

t = len(taxa)
del trees

# initialize multiprocessing objects
queue = Queue()
randdmfilenames_array = Array(c_char_p,len(randfnamelist))
for i, fname in enumerate(randfnamelist):
	randdmfilenames_array[i] = fname

taxlabels_array = Array(c_char_p,len(taxa))
for i, taxlabel in enumerate([taxon.label for taxon in taxa]):
	taxlabels_array[i] = taxlabel

# calculate distances from randomized trees
print "\nestimating RANDOMIZED distances:"
processes = [Process(target=calc_I_partial, args=(randdmfilenames_array,taxlabels_array,nthreads,i,queue)) for i in range(nthreads)]

for p in processes:
	p.start()

I_random = dict()
retrieved = 0
while retrieved < nthreads:
	try:
		I_partial = queue.get()
		if I_partial.__class__ is dict:
			for key, val in I_partial.iteritems():
				try:
					I_random[key] += val
				except KeyError:
					I_random[key] = val
			retrieved += 1
	except Empty:
		pass

for p in processes:
	p.join()

#print I_random

I_avg_est = float(sum(I_random.values())) / (len(I_random) * ((nrandtrees * (nrandtrees - 1)) / 2))
print "\n\naverage expected distance estimated by random taxon placement is %s\n" % I_avg_est

# calculate distances from observed trees 
print "\ncalculating instability scores:"

dmfilenames_array = Array(c_char_p,len(fnamelist))
for i, fname in enumerate(fnamelist):
	dmfilenames_array[i] = fname

processes = [Process(target=calc_I_partial, args=(dmfilenames_array,taxlabels_array,nthreads,i,queue)) for i in range(nthreads)]

for p in processes:
	p.start()

I_sum = dict()
retrieved = 0

while retrieved < nthreads:
	try:
		I_partial = queue.get()
		if I_partial.__class__ is dict:
			for key, val in I_partial.iteritems():
				try:
					I_sum[key] += val
				except KeyError:
					I_sum[key] = val
			retrieved  += 1
	except Empty:
		pass

for p in processes:
	p.join()

#print I_sum

# start recording scores
logfile.write("Leaf,Mean per-tree distance,Score,,Avg expected distance = %s\n" % I_avg_est)
logfile.flush()

# generate relative values by scaling average per-tree score to expected maximum for this value
I_score = dict()
for t in I_sum:
#	I_score[t] = float(I_sum[t]) / ntrees / I_avg_est
	I_score[t] = float(I_sum[t]) / (I_avg_est * ((ntrees * (ntrees - 1)) / 2))

	# record final values
	logfile.write(str(t) + "," + str(I_sum[t] / ntrees) + "," + str(I_score[t]) + "\n")

logfile.close()

print "done."