#!/usr/bin/env python2.6
#
#  instability.py
#
#  Created by Cody on 2/16/11.
#

# import python modules
import os
import sys
import dendropy as dp
import cPickle as pickle
from copy import copy
from math import factorial as fac

def calc_I(taxa,distance_matrix_list=None,file_list=None):

	I = dict()

	if distance_matrix_list is None and file_list is None:
		exit("Calc function requires matrices or temp files")

	elif distance_matrix_list is not None and file_list is None:
		dmlist = distance_matrix_list
		# for every taxon t1
		for i, t1 in enumerate(taxa):
			I[t1] = 0
			print str(t1)
			# compare t1 to to every other taxon t2
			for t2 in taxa[0:i] + taxa[i+1:t]:			
				# for every distance matrix x
				for j, dm_x in enumerate(dmlist):
					# find distance for t1, t2 from x
					D_ijx = dm_x(t1,t2)
					# sum differences in distance for taxa t1, t2 between matrix x and all other matrices y
					for dm_y in dmlist[j+1:]:
						# find distance difference for (t1, t2), between matrices x and y
						D_ijy = dm_y(t1,t2)
						# add distance difference to sum
						I[t1] += abs(D_ijx - D_ijy)
			print I[t1]
	
	elif distance_matrix_list is None and file_list is not None:
		dmfiles = file_list

		for taxon in taxa:
			I[taxon.label] = 0

		# for every distance matrix x
		for i, file_x in enumerate(dmfiles):
			print "Currently comparing dists against " + file_x
			print "\t" + str((len(dmfiles)-(i+1))) + " matrix comparisons to be made against this matrix"
			# read matrix x from file
			dm_x = pickle.load(open(file_x,"rb"))

			# sum differences in distance for taxa t1, t2 between matrix x and all other matrices y
			for file_y in dmfiles[i+1:]:
				# read matrix y from file
				dm_y = pickle.load(open(file_y,"rb"))

				# for every taxon t1
				for j, taxon1 in enumerate(taxa):
					t1 = taxon1.label

					# compare t1 to to every other taxon t2
					for taxon2 in taxa[0:j] + taxa[j+1:t]:			
						t2 = taxon2.label

						# taxon sets don't match after pickling, so we have to look up taxa in post-pickled matrices using labels
						tax_x1 = dm_x.taxon_set.get_taxon(label=t1)
						tax_x2 = dm_x.taxon_set.get_taxon(label=t2)
						tax_y1 = dm_y.taxon_set.get_taxon(label=t1)
						tax_y2 = dm_y.taxon_set.get_taxon(label=t2)

						# find distances for t1, t2 from matrices x and y
						D_ijx = dm_x(tax_x1,tax_x2)
						D_ijy = dm_y(tax_y1,tax_y2)
						
						# add distance difference to sum
						I[t1] += abs(D_ijx - D_ijy)
						
	else:
		exit("Calc function accepts a list in memory or a list of file, but not both")
	
	return I

# this should eventually be set via command line or config file
useHD = True

# attempt to open the treefile
try:
	infname = sys.argv[1]
	treefile = open(infname,"rb")
except IndexError:
	exit("You must provide the path to a treefile.")
except IOError:
	exit("The file '%s' could not be opened." % infname)

# second command line argument provides an optional logfile name 
try:
	loglabel = sys.argv[2]
except IndexError:
	loglabel = sys.argv[1]

# attempt to open logfile
try:
	logfname = loglabel+".instability.csv"
	logfile = open(logfname,"wb")
except IOError:
	"You must specify a different name for the log file."

# optional: set to true when edges in treefile have no lengths
make_edge_lengths_uniform = True

# define a taxon set to be shared by all trees
taxa = dp.TaxonSet()

# create a directory to hold temp files (for low-memory option, currently not used)
if useHD:
	try:
		wdir = os.path.dirname(os.path.abspath(sys.argv[1]))
		temppath = wdir+"/instability_temp"
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

if useHD:
	fnamelist = list()
	randfnamelist = list()
else:
	dmlist = list()
	randdmlist = list()

j = 0
# calc distance matrices for all trees
for i, curtree in enumerate(trees):

	# if edges have no lengths, set all edge lengths to one 
	if make_edge_lengths_uniform:
		for n in curtree.get_node_set():
			n.edge_length = 1
	
	# find all patristic distances for each tree:

	dmatrix = dp.treecalc.PatristicDistanceMatrix(curtree)

	# These matrices are very large. A potentially necessary improvement would be
	# to save each matrix to a file instead of appending them to a list. We would
	# then need to read the files from disk for every taxon pair to perform the
	# summary calculations.
	#
	# Because we need to compare each matrix to each other, it doesn't really make
	# sense to save more than one matrix per file, so the number of file reads would
	# be (ntrees)! / 2(ntrees - 2)!. It shouldn't too bad of a speed decrease if it
	# allows us do the task with minimal memory overhead.

	if useHD:
		tempfile = open(temppath + "/%s.temp" % i, "wb")
		pickle.dump(dmatrix,tempfile)
		fnamelist.append(os.path.abspath(tempfile.name))
		tempfile.close()
	else:
		dmlist.append(dmatrix)
	
	# randomize each rfraction-th tree to use for estimating max theoretical distances
	if i % rfraction == 0:
		randtree = dp.Tree(curtree)
		randtree.randomly_assign_taxa(create_required_taxa=False)
		randdmatrix = dp.treecalc.PatristicDistanceMatrix(randtree)
		
		if useHD:
			randtempfile = open(temppath + "/random_%s.temp" % j, "wb")
			pickle.dump(randdmatrix,randtempfile)
			randfnamelist.append(os.path.abspath(randtempfile.name))
			randtempfile.close()
		else:
			randdmlist.append(randdmatrix)

		j += 1
	
	print "stored distance matrix for tree %s" % i

t = len(taxa)

# calculate distances from randomized trees
print "\nestimating RANDOMIZED distances:"
if useHD:
	I_random = calc_I(taxa,file_list = randfnamelist)
else:
	I_random = calc_I(taxa,distance_matrix_list = randdmlist)

I_avg_est = sum(I_random.values()) / len(I_random) / nrandtrees
print "\n\naverage expected per-tree distance estimated by random taxon placement is %s\n" % I_avg_est

# calculate distances from observed trees 
print "\ncalculating instability scores:"
if useHD:
	I_sum = calc_I(taxa,file_list = fnamelist)
else:
	I_sum = calc_I(taxa,distance_matrix_list = dmlist)

# start recording scores
logfile.write("Leaf,Mean per-tree distance,Score,,Avg expected distance = %s\n" % I_avg_est)

# generate relative values by scaling average per-tree score to expected maximum for this value
I_score = dict()
for t in I_sum:
	I_score[t] = float(I_sum[t]) / ntrees / I_avg_est

	# record final values
	logfile.write(str(t) + "," + str(I_sum[t] / ntrees) + "," + str(I_score[t]) + "\n")

logfile.close()

print "done."