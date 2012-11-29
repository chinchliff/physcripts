#!/usr/bin/env python2.3
#
#  combinations.py
#  
#
#  Created by Cody on 2/16/11.
#  Copyright (c) 2011 The Evergreen State College. All rights reserved.
#

def combine2(list_items):
	# make a list containing all 2-item combinations for a given list of values
	
	combs = list()
	l = list_items
	
	while len(l) > 0:
		cur = l.pop()
		for next in l:
			combs.append(cur+next)

	return combs