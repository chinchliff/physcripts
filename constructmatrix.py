import sys
import os
import cfg_file
from Bio import SeqIO

class marker:
    def __init__(self, mname, mfile):
        self.name = mname
        self.path_to_file = mfile

def compare_valid_taxa(x, y):
    if x.__class__ is marker.__class__ and y.__class__ is marker.__class__:
        if x.prop_valid_taxa > y.prop_valid_taxa:
            return 1
        elif x.prop_valid_taxa == y.prop_valid_taxa:
            return 0
        else: # x < y
            return -1

config = cfg_file.reader("construct_matrix.cfg")

use_list_file = False
list_file = None
markers = list()

adding_marker = False

for p in config.parameters:

    if not adding_marker:
        if p.name == "sampling_rank":
            sampling_rank = p.value
        elif p.name == "start_node":
            start_node = p.value
        elif p.name == "path_to_list_file":
            path_to_list_file = p.value
        elif p.name == "use_list_file" and p.value.lower() == 'true':
            use_list_file = True

        elif p.name == "marker":
            adding_marker = True
            newmarkername = p.value
    else:
        adding_marker = False
        if p.name == "file":
            markers.append(marker(newmarkername,p.value))

## Now get a list of all taxa of rank -> sampling_rank. For now,
## we're just feeding it a pre-built list. In the future we'll need
## an option to pull the list from the local database that PHLAWD
## is using

if use_list_file:
    list_file = open(path_to_list_file, "rU")
    taxon_list = list_file.readlines()

    alltaxa = list()

    for taxon in taxon_list:
        alltaxa.append(taxon.strip())

num_taxa_total = len(alltaxa)

#### else:

for marker in markers:

    try:
        marker.fasta_file = open(marker.path_to_file,'rU')
        marker.data = list(SeqIO.parse(marker.fasta_file,"fasta")) # make this a list so we can iterate over it more than once
    except IOError:
        print "There was a problem opening the file at:\n" + marker.path_to_file

    marker.valid_taxa = list()

###### For now, simplifying this by just looking for genera, which can be
## gleaned from taxon names. In future, need to look up each organism and
## find the appropriate taxon of rank sampling_rank to append to the list

    for seqRecord in marker.data:
        genus = seqRecord.name.split('_')[0]
        
        if genus in alltaxa and genus not in marker.valid_taxa:
            marker.valid_taxa.append(genus)

######

    marker.num_valid_taxa = len(marker.valid_taxa)
    marker.prop_valid_taxa = float(marker.num_valid_taxa) / num_taxa_total

print "\n-----------------------------------------------------\n" \
      "Characterizing data sets:\n"

for marker in markers:
    print marker.name
    print "\tNo. valid taxa found: %i" % marker.num_valid_taxa
    print "\tProp. valid taxa represented: %.2f%% " % marker.prop_valid_taxa

# sort markers in decreasing order by the proportion of valid taxa they contain
markers.sort(compare_valid_taxa,reverse=True)

# make a list of all taxa present in highest ranking dataset, which we will whittle down
common_taxa = [sequence.name for sequence in markers[0].data]
start_count = len(common_taxa)

print "\n-----------------------------------------------------\n" \
      "Removing taxa not present across all datasets:\n\n" \
      "Starting with %i taxa in %s\n" % ( start_count, markers[0].name )

last_count = start_count

for marker in markers:
    t = [sequence.name for sequence in marker.data] # get list of taxa in this dataset
    for taxon in list(common_taxa): # copy the list so we can remove items from the original
        if taxon not in t:
            common_taxa.remove(taxon) # remove taxa from common_taxa if not in this dataset

    cur_count = len(common_taxa)

    print "Deleted %i taxa after comparing against %s (retaining %i taxa)" % \
          ( last_count - cur_count, marker.name, cur_count)

    last_count = cur_count

print "\nTaxa present in all datasets: "

for t in common_taxa:
    print "\t" + t

print "\n"







    
