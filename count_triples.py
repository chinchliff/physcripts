#!/usr/bin/env python2.3
#
#  count_triples.py
#  
#
#  Created by Cody on 2/13/11.
#  Copyright (c) 2011 The Evergreen State College. All rights reserved.
#

import numpy as np

h = 5
w = 2
matr = np.empty((h,w), dtype=bool)

matr[0,0] = True
matr[1,0] = True
matr[2,0] = True
matr[3,0] = True
matr[4,0] = False
matr[0,1] = True
matr[1,1] = True
matr[2,1] = True
matr[3,1] = False
matr[4,1] = True

print matr

# identify satisfied triples by column
col = 0
while col < w:
	row = 0
	while row < h:
		triples = list()
		if matr[row,col] == True:
			for 