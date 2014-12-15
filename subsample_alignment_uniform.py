#!/usr/bin/env python

'''Subsample an alignment uniformly at random.'''

DEFAULT_LABEL = 'subsampled_uniform'

class Alignment():

    def __init__(self, alignment, partitions):
        self._parse_partitions(partitions)
        self._parse_alignment(alignment)

    def _parse_partitions(self, partitions_file):
        if partitions_file is None:
            self._parts = None
            return

        self._parts = {}
        self._parts_names = []
        for line in partitions_file:
            toks = [t.strip() for t in re.split(r'[,=]+',line)]
            ptype = toks[0]
            name = toks[1]
            bounds = [b.strip() for b in toks[2].split("-")]
            start = int(bounds[0])
            end = int(bounds[1])
            self._parts[start] = { 'name': name, 'type': ptype, 'start': start, 'end': end, 'data': {}, }
        
    def _parse_alignment(self, alignment_file):
    
        self.base_name = os.path.basename(alignment_file.name).rsplit('.phy',1)[0]
    
        ntax = 0
        self.ncols = 0
        self._taxa = []
        on_first_line = True
        for line in alignment_file:
            toks = [l.strip() for l in line.split()]
            if len(toks) > 1:
                if on_first_line == True:
                    on_first_line = False
                    ntax = int(toks[0])
                    self.ncols = int(toks[1])
                    if self._parts is None:
                        self._parts = {1: { 'name': 'all', 'type': 'DNA', 'start': 1, 'end': int(self.ncols), 'data': {}, }}
                else:
                    name = toks[0]
                    seq = toks[1]
                    self._taxa.append(name)
                    for s in self._parts.keys():
                        e = self._parts[s]['end']
                        self._parts[s]['data'][name] = seq[s-1:e]
            else:
                raise IndexError("too many items on line '" + line + "' in alignment") 
    
        assert(ntax == len(self._taxa))
    
    def random_part(self):
        return list(self._parts.keys())[random.randrange(len(self._parts))]

    def ntaxa(self):
        return len(self._taxa)

    def nparts(self):
        return len(self._parts)

    def taxa(self):
        return self._taxa
        
    def parts(self):
        return self._parts.keys()

class Subsample():

    def __init__(self, alignment):
        import random
        self.alignment = alignment
        self.init_subsample()
    
    def init_subsample(self):
        self._taxa_sampled_per_part = dict(zip(self.alignment.parts(), [0,] * self.alignment.nparts()))
        self._parts_sampled_per_taxon = dict(zip(self.alignment.taxa(), [0,] * self.alignment.ntaxa()))
        self._sample_bitmap = {}
        
    @staticmethod
    def validate_p(p):
        p = float(p)
        if not p > 0 and p < 1:
            raise TypeError('sampling proportion d must satisfy 0 < d < 1')
        return p

    def set_random_seed(self, r):
        if r is not None:
            random.seed(r)
        else:
            random.seed()

    @staticmethod
    def _uniform_probability_function(p): #, t, s):
        return random.random() < p
        
    def subsample(self, p):
    
        from copy import deepcopy as copy
    
        p = Subsample.validate_p(p)
        self.init_subsample()
        
        k = 0 # total count of sampled sites
        for t in self.alignment.taxa():
            if t not in self._sample_bitmap:
                self._sample_bitmap[t] = {}

            for s in self.alignment.parts():
                if s not in self._sample_bitmap[t]:
                    self._sample_bitmap[t][s] = {}
                
                if Subsample._uniform_probability_function(p): # could be replaced with generic probability function
                    self._sample_bitmap[t][s] = True
                    self._taxa_sampled_per_part[s] += 1
                    self._parts_sampled_per_taxon[t] += 1
                    k += 1
                else:
                    self._sample_bitmap[t][s] = False

            # open a random site, to be sure we have at least one sampled partition for all taxa
            if (self._parts_sampled_per_taxon[t] < 1):
                r = self.alignment.random_part()
                if self._sample_bitmap[t][r] != True:
                    self._parts_sampled_per_taxon[t] += 1
                self._sample_bitmap[t][r] = True

        self.sampling_proportion = float(k) / (self.alignment.nparts() * self.alignment.ntaxa())


    # iterators over the dicts to return data sorted by sampling proportion
    def taxa_sorted(self):
        return sorted(self.alignment.taxa(), \
            cmp = lambda p, q: cmp(self._parts_sampled_per_taxon[p], \
                                   self._parts_sampled_per_taxon[q]), reverse=True)

    def parts_sorted(self):
        return sorted(self.alignment.parts(), 
            cmp = lambda p, q: cmp(self._taxa_sampled_per_part[p], \
                                   self._taxa_sampled_per_part[q]), reverse=True)


    def write_subsampled_output(self, l):
        self.output_label = self.alignment.base_name + '.' + (DEFAULT_LABEL + '.' + l if len(l) > 0 else DEFAULT_LABEL)
        
        with open(self.output_label + '.sampling_matrix.tsv', 'w') as s:
            s.write('\t'+'\t'.join([self.alignment._parts[n]['name'] for n in self.parts_sorted()])+'\n')
            for t in self.taxa_sorted():
                s.write(t)
                for p in self.parts_sorted():
                    s.write('\t')
                    if self._sample_bitmap[t][p]:
                        s.write('1')
                    else:
                        s.write('0')
                s.write('\n')

        with open(self.output_label + '.phy', 'w') as a:
            a.write('{} {}\n'.format(self.alignment.ntaxa(), self.alignment.ncols))
            for t in self.taxa_sorted():
                a.write(t + ' ')
                for p in self.parts_sorted():
                    if self._sample_bitmap[t][p]:
                        a.write(self.alignment._parts[p]['data'][t])
                    else:
                        a.write('-' * len(self.alignment._parts[p]['data'][t]))
                a.write('\n')
    
if __name__ == '__main__':

    import argparse, copy, operator, os, random, re, sys 
    
    parser = argparse.ArgumentParser(description=__doc__)
    
    parser.add_argument('-a', '--alignment', type=open, required=True, \
        help='the location of the extended phylip alignment file to be subsampled.')
    
    parser.add_argument('-p', '--sampling-proportion', type=float, required=True, \
        help='a decimal value d | 0 < d < 1 to use as the sampling proportion. NOTE: no matter ' \
             'what value of d is provided, at least one sequence will be retained for every ' \
             'partition in the alignment.')
    
    parser.add_argument('-q', '--partitions', type=open, required=False, \
        help='the location of the raxml partitions file corresponding to the alignment to be subsampled.')
    
    parser.add_argument('-x', '--random-seed', type=int, required=False, \
        help='an integer seed for the random number generator function')
    
    parser.add_argument('-n', '--output-label', required=False, default='', \
        help='a label to be attached to output files')

    args = parser.parse_args()
    
    a = Alignment(args.alignment, args.partitions)
    s = Subsample(a)
    s.set_random_seed(args.random_seed)
    s.subsample(args.sampling_proportion)

    s.write_subsampled_output(args.output_label)
    
    print('files have been written to: ' + s.output_label + '.sampling_matrix.txt, ' + s.output_label + '.subsampled.phy\n' \
              'sampling proportion is ' + str(s.sampling_proportion))