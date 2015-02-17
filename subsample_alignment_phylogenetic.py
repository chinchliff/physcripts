#!/usr/bin/env python

'''Subsample an alignment using a tree.'''

from phyaln import Alignment, PhylogeneticSubsampler

def read_rate_file(f):
    rates = {}
    for parts in (l.split('=') for l in f):
        if len(parts) == 2:
            name, rate = (p.strip() for p in parts)
            rates[name] = float(rate)
    return rates
    
if __name__ == '__main__':

    import argparse, newick3 
    
    parser = argparse.ArgumentParser(description=__doc__)
    
    parser.add_argument('-a', '--alignment', type=open, required=True, \
        help='the location of the extended phylip alignment file to be subsampled.')

    parser.add_argument('-r', '--rates', type=open, required=True, \
        help='the evolutionary rates (substitution/time) underlying the models for the partitions.')
    
    parser.add_argument('-t', '--tree', type=open, required=True, \
        help='the tree.')

    parser.add_argument('-q', '--partitions', type=open, required=True, \
        help='the location of the raxml partitions file corresponding to the alignment to be subsampled.')
            
    parser.add_argument('-x', '--random-seed', type=int, required=False, \
        help='an integer seed for the random number generator function')
    
    parser.add_argument('-n', '--output-label', required=False, default='', \
        help='a label to be attached to output files')

    args = parser.parse_args()
    
    a = Alignment(args.alignment, args.partitions)
    t = newick3.parse(args.tree)
    r = read_rate_file(args.rates)
    
    s = PhylogeneticSubsampler(alignment=a, tree=t, rates=r)
    
    s.subsample()

    s.write_subsampled_output(args.output_label)
    
    s.report_sampled_partitions()
    
    print('files have been written to: ' + s.output_label + '.sampling_matrix.txt, ' + s.output_label + '.subsampled.phy\n' \
              'sampling proportion is ' + str(s.get_sampling_proportion()))
