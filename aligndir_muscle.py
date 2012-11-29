#!/usr/bin/env python2.6

from subprocess import call
import os

cwd = "/Users/cody/Desktop/its_test"
files = os.listdir(cwd)

for file in files:
	try:
		call("muscle -in %s/%s -out %s/%s.aligned" % (cwd, file, cwd, file), shell = True)
	except OSError:
		pass
