#!/usr/bin/env python
#  Created by Cody on 10/30/11.

import sys
import os
import subprocess
import dendropy as dp

#treefname = "RAxML_bootstrap.has_ndhf_or_rbcl_stable_tips"
#outgroup_taxon_names = ["Diplasia karatifolia", "Hypolytrum nemorum", "Hypolytrum longifolium", "Hypolytrum bullatum", "Mapania meditensis", "Mapania cuspidata", "Scirpodendron ghaeri", "Lepironia articulata", "Chrysitrix dodii", "Chrysitrix capensis", "Chorizandra enodis", "Chorizandra cymbaria", "Paramapania parvibractea", "Mapania paradoxa", "Mapania macrophylla"]	

def root_trees():
	global treefname

	taxa = dp.TaxonSet()
	treelist = dp.TreeList()
	treelist.read_from_path(treefname, schema="newick", taxon_set=taxa)

	global 	outgroup_taxon_names
	outgroup_taxa = list()

	for name in outgroup_taxon_names:
		for t in taxa:
			print t.label
			if t.label == name:
				outgroup_taxa.append(t)

	print outgroup_taxa

	for tree in treelist:
		rootnode = tree.mrca(taxa=outgroup_taxa)
		tree.reroot_at_edge(rootnode.edge, length1 = rootnode.edge_length / 2 , length2 = rootnode.edge_length / 2, update_splits = True)
		tree.print_plot()

	outfile = open(treefname + ".rooted", "wb")
	treelist.write(outfile,schema="newick", edge_lengths = True)
	rooted_trees_fname = outfile.name
	outfile.close()

#########################################################

def collect_treefile(treefile_name):

	treefile = open(treefile_name, "rb")
	
	trees = dict()
	for line, index in enumerate(trees):
		trees[index] = line
	
	return trees

def collect_treedir(treefile_dir):

	trees = dict()
	for treefilename in os.listdir(treefile_dir):
		thisfile = open(treefile_dir + "/" + treefilename)
		trees[treefilename] = thisfile.readline()
	
	return trees

#########################################################

def calc_pathd8(trees,seqlen,additional_lines):

	global temppath

	try:
		os.mkdir(temppath)
	except OSError:
		pass

	try:
		os.mkdir(temppath + "/infiles")
	except OSError:
		p = os.listdir(temppath + "/infiles")
		if len(p) > 0:
			for fname in p:
				os.unlink(temppath + "/infiles/" + fname)

	try:
		os.mkdir(temppath + "/outfiles")
	except OSError:
		p = os.listdir(temppath + "/outfiles")
		if len(p) > 0:
			for fname in p:
				os.unlink(temppath + "/outfiles/" + fname)

	for name, tree in trees.iteritems():
		infile = open(temppath + "/infiles/" + name + ".pathd8","wb")
		infile.write("Sequence length = " + str(seqlen) + "; \n" + tree + "\n" + additional_lines)
		infile.close()

	infiles = os.listdir(temppath + "/infiles")

	for infname in infiles:
		infpath = temppath + "/infiles/" + infname
		outfpath = temppath + "/outfiles/" + infname + ".out"
		subprocess.call(["./PATHd8", infpath, outfpath])
		print infname

#########################################################

def extract_pathd8_trees(save_as_files):
	
	global temppath
	
	if save_as_files == True:
		pathd8_treedir_path = temppath + "/pathd8_trees" 
		try:
			os.mkdir(pathd8_treedir_path)
		except OSError:
			pass
	else:
		pathd8_treefile = open(temppath + "/pathd8_trees.tre","wb")
	
	for fname in os.listdir(temppath + "/outfiles"):
		found_tree = False
		file = open(temppath + "/outfiles/" + fname, "rb")
		print file.name
		for line in file:
			if not found_tree:
				if line[0:13] == "d8 tree    : ":
					tree = line[13:]
					
					if save_as_files == True:
						d = open(pathd8_treedir_path + "/" + fname + ".d8.tre","wb")
						d.write(tree)
						d.close()					
					else:
						pathd8_treefile.write(tree)
					
					found_tree = True
	
	if save_as_files != True:
		pathd8_treefile.close()


temppath = "/Users/cody/Desktop/cyp_ltt_plots/besttrees_rooted"
additional_lines = "mrca: Carex_pulicaris, Mapania_meditensis, minage=87;"
seqlen = 16016

trees = collect_treedir(temppath + "/treefiles")
calc_pathd8(trees, seqlen, additional_lines)
extract_pathd8_trees(save_as_files=True)
