#!/usr/bin/env python3

if __name__ == '__main__':

    import os, sys, re

    if len(sys.argv) < 2:
        print("usage: match_taxa_in_fasta.py <taxonlist> <fastafile>")
        sys.exit(0)

    taxa = {}
    with open(sys.argv[1], "r") as taxon_list_file:
        for t in taxon_list_file.readlines():
            tax_name = re.sub("[.!~@#$%^&*()+`{}\[\]:\";'<>?|,.\\/ _]+", '_', t.strip()).strip("_")
            taxa[tax_name] = ""

    with open(sys.argv[2], "r") as fasta_file:
        for line in fasta_file.readlines():
            if line.strip()[0] == ">":
                name = re.sub("[.!~@#$%^&*()+`{}\[\]:\";'<>?|,.\\/ _]+", '_', line.strip(">")).strip().strip("_")
                if name in taxa:
                    print(name)
