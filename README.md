## Basic documentation
    
This directory contains a bunch of scripts to automate common (or not so common) tasks, mostly targeted at bioinformatics and phylogenetic systematic research. To install the scripts, clone the git repo into a local directory and add that directory to your PATH and PYTHONPATH environment variables. E.g.:

```bash
# clone the scripts repo
cd ~/ && git clone https://github.com/chinchliff/physcripts

# add it to paths
echo 'export PATH=/Users/cody/physcripts:$PATH' >> ~/.bash_profile
echo 'export PYTHONPATH=/Users/cody/physcripts:$PYTHONPATH' >> ~/.bash_profile

# load new environment variables
source ~/.bash_profile
```

## Running examples

The descriptions below include a variety of example calls to the scripts. The input files used in these examples are included in the `test_files` directory. The example calls are intended to be run from within the physcripts directory, but you can run the elsewhere by if you modify the calls to indicate the correct path to the example input files.

## Scripts

What follows is a (complete) list of the available Python scripts, along with (highly incomplete) descriptive documentation about each one. More information can usually be found in the form of comments within the scripts (in addition, of course, to the hopefully somewhat self-documenting code itself).



---

### Add tips to a chronogram

`add_species_to_chronogram`

Adds names from a file into to an tree, while maintaining the depths for all affected internal nodes, and ensuring that newly added tips will have depths equal to the depth of their sister clades. In other words, this will maintain the time-relations in the chronogram as new tips are added.
```bash

add_species_to_chronogram.py -d ~/phylo/data/gbpln_2014_04_25.db -n test_files/names_and_synonyms.txt 
```

`clean_gb_names_simple_fasta`

`clean_phlawd_spp_tip_names`

`excise_knuckles`

---

### Get information about names from a PHLAWD database

`extract_names_and_ids_from_phlawd_db`

Given a set of names, extract information about those names from the specified PHLAWD database.
```bash

extract_names_and_ids_from_phlawd_db.py -t test_files/tree_with_50_tips.tre -n test_files/names_50.txt -s -b 0.1 
```

`extract_names_from_roguenarok_output`

`figtree_blocks`

`fill_tips_from_taxonomy`

`filter_fasta`

Usage:

./filter_fasta.py <path to input dir> <path to accepted taxon list>

Input files are expected to be in fasta format. The script will traverse all files
in the input dir, so the input dir should contain only fasta files. The taxon list
should be a line-delimited text file containing the names of tips as they
correspond to those in the fasta alignments.

`filter_list`

`filter_phylip`

filters a phylip alignment; can accept *either* a set of names to accepted or a set to be excluded, and saves the output to another file

`find_names_missing_from_tree`

`fix_gb_names_in_fasta`

Cleans the namestrings in fasta files downloaded from ncbi, saving only minimal name information with no special characters. Optionally will use an interactive prompt to allow the user to rename things with names in unrecognized formats

`get_nested_taxa_for_names`

uses a phlawd sqlite3 taxonomy database to retrieve all the names associated with each taxon 
in the provided list. <targetrank> specifies the rank of the taxa to be returned.

`make_bootstraps`

`make_constraint_tree`

`make_partitions_from_pb_settings`

`make_random_tree`

`make_readme`

Generate simple markdown documentation for the scripts in this directory

`make_sampling_matrix_from_alignment`

`make_sampling_matrix_from_fastas`

This script accesses a directory, and traverses all FASTA files in it, recording the names of all taxa present
in each file. Then it creates a tab-delimited file containing a matrix where the rows represent the taxa and the
columns the FASTA files. The intended use is for a directory containing a set of FASTA files each corresponding
to a single locus, and containing homologous sequences of that locus for different taxa. The script will record
a 1 in the resulting matrix if a taxon is present in a locus file, or a 0 if not.

Key point: the script does not intelligently differentiate FASTA files from other types, and it will attempt
to parse any file in the directory. For this reason, you should remove all other files before you run the script.

It will create (or overwrite!) a file in the passed directory called 'sampling_matrix.txt' that may be opened
in any conventional spreadsheet or text-editor app. This file should be in the proper format for use in the
Decisivator application.

This script requires BioPython.

`make_tip_list_by_genus`

Reads a newick tree with phlawd-style tip names -- <ncbi_id>_<genus>_<spepithet> -- and creates a line-delimited list of the generic names in the tree and their constituent species that is written to outfile. Names are just extracted and parsed with regex search, does not use the tree structure at all.

`match_taxa_in_fasta`

`merge_lists`

`newick3`

Classes and methods for performing basic operations on Newick trees.
usage: nexus2fasta_dirty <nexusfile> <fastaoutputfilename>
usage: nexus2fasta_dirty <nexusfile> <fastaoutputfilename>

`paint_branches`

`phylip2fasta`

`phylo3`

Classes and methods for performing basic operations on phylogenetic trees.

`print_tip_names`

`prune_long_tips`

`prune_tips`

`root_tree_against_master`

Place the root of a target tree in a position determined by the relationships present in a given master tree. The following conditions must be met:

1. The master tree contains all the taxa in the target tree.
2. Some bipartition X = A|B exists in the master tree such that for some bipartition Y = C|D in the target, C is a subset of A and D is a subset of B

If this case is satisfied, the target tree will be rooted at Y.

`root_trees`

Will root a set of newick trees using a line-delimited list of taxon names to be included in the outgroup.

`strip_node_labels`

`subsample_edge_quartets`

Consider each node in the rooted tree to identify a bipartition, which is represented in
the tree as the outgoing edge connecting the node to its parent. In a fully bifurcating tree,
each node connected to this edge (the child and the parent) will have two other connected edges.
If we unroot the tree and consider the outgoing direction to be away from the bipartition of
interest, then leaf sets descendant from these outgoing edges form the four sets of the quartet
induced by the bipartition. We perform a number of replicated tests consisting of randomly drawing
one taxon from each of these four quartets, reconstructing the topology for the four random tips
using the sequence data with sequence data from the original alignment, and we record the topology
for each replicate. The resulting topology sets are used to calculate the ICA score for the 
bipartition.

`test_monophyly_against_tree`

`transpose_fastas`

transpose_fastas.py
version 0.1
cody hinchliff
2012.7.21

This script accesses a directory, and traverses all FASTA files in it, recording them into a dict object. It then
writes new fasta files to a directory named 'inverted' within the provided dir. One inverted fasta file is created
for each sequence id found in the original set of fastas; each of these files contains the sequences associated
with this sequence id, labeled with ids representing the original fasta files whence they came.
	
For example, if the script is passed a directory containing fasta files corresponding to loci, containing sequences
labeled with taxon names, the inverted directory will contain fasta files corresponding to taxon names, containing
sequences labeled with locus names (drawn from the filenames of the original fastas).
