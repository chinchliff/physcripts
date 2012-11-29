#!/usr/bin/env python2.3
#
#  count_4trees.py
#  
#
#  Created by Cody on 2/11/11.
#  Copyright (c) 2011 The Evergreen State College. All rights reserved.
#

import copy
import dendropy
import sys

def combine_2(list_items):
	# make a list containing all 2-item combinations for a given list of values
	
	combs = list()
	l = list_items
	
	while len(l) > 0:
		cur = l.pop()
		for next in l:
			combs.append(cur+next)

	return combs

def make_pairsets(list_items_1,list_items_2):
	# given two lists, make a list of 2-item tuples containing
	# every possible pair consisting of one item from each list 

	pairs = list()
	l1 = list_items_1

	while len(l1) > 0:
		i1 = l1.pop()
		l2 = copy.copy(list_items_2)

		while len(l2) > 0:
			i2 = l2.pop()
			pairs.append((i1,i2))
	
	return pairs

def uniquify(the_list_to_uniquify):
	# remove duplicate entries from a given list
	
	all = the_list_to_uniquify
	uniques = list()
	while len(all) > 0:

		cur = all.pop()

		if cur not in all:
			uniques.append(cur)

	return uniques

# attempt to open the treefile
try:
	tree = dendropy.Tree(stream=open(sys.argv[1]), schema = "nexus")
except IndexError:
	exit("You must provide the path to a treefile.")
except IOError:
	exit("The file '%s' could not be opened." % sys.argv[1])

# get a list of all internal nodes
internal_nodes = tree.internal_nodes()

# for each internal node, get a list of taxa to the right/left of node. these are the splits
splits = list()
for n in tree.internal_nodes():
	child_leaves = n.leaf_nodes()
	a = [n.taxon.label for n in child_leaves]

	b = list()
	for m in tree.leaf_nodes():
		if m not in child_leaves:
			b.append(m.taxon.label)
	
	if len(a) > 1 and len(b) > 1:
		splits.append((a,b))

# format the splits into a nice list for logging
strsplits = ""
for s in splits:
	strsplits += str(s[0]) + " | " + str(s[1]) + "\n"

print "The current tree contains the following splits: \n" + strsplits

all = list()
for split in splits:

	# calculate all 2-taxon combinations possible on each side of the current split 
	c1 = combine_2(split[0])
	c2 = combine_2(split[1])

	# make all possible pairs of the 2-taxon combinations, these are the 4-taxon subtrees 
	pairs = make_pairsets(c1,c2)

	# append the subtrees to the master list
	all += pairs

# remove redundant subtrees
uniques = uniquify(all)

print "the number of unique four-taxon subtrees is %s:" % str(len(uniques))
print uniques