#!/usr/bin/env python
#
#  bootstrap_manip.py
#  
#
#  Created by Cody on 6/15/11.
#  Copyright (c) 2011 The Evergreen State College. All rights reserved.

import dendropy

trees = dendropy.TreeList.get_from_path("/Users/cody/Desktop/allopt.tre","newick")

mapanioideae_names = ("Mapania meditensis", "Mapania cuspidata", "Scirpodendron bogneri", "Hypolytrum bullatum", "Hypolytrum longifolium", \
					  "Hypolytrum nemorum", "Diplasia karatifolia", "Paramapania parvibractea", "Mapania macrophylla", "Lepironia articulata", \
					  "Chorizandra cymbaria", "Chorizandra enodis", "Chrysitrix dodii", "Chrysitrix capensis")

rootedtrees = list()
for tree in trees:

	mapanioideae = tree.mrca(taxon_labels=mapanioideae_names)
	newbrlen = mapanioideae.edge_length / 2
	tree.reroot_at_edge(mapanioideae.edge,length1=newbrlen,length2=newbrlen,update_splits=True)
	
	rootedtrees.append(tree)

rootedtrees_tlist = dendropy.TreeList(rootedtrees)
rootedtrees_tlist.write_to_path("/Users/cody/Desktop/allopt.rooted.tre","newick")