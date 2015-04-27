#!/usr/bin/env python

class Alignment():

    def __init__(self, alignment, partitions):
        self._parse_partitions(partitions)
        self._parse_alignment(alignment)

    def _parse_partitions(self, partitions_file):

        import re

        self._parts_by_name = {}

        if partitions_file is None:
            return

        for line in partitions_file:
            toks = [t.strip() for t in re.split(r'[,=]+',line)]
            ptype = toks[0]
            name = toks[1]
            bounds = [b.strip() for b in toks[2].split("-")]
            start = int(bounds[0])
            end = int(bounds[1])
            p = Partition(alignment = self, name = name, type = ptype, start = start, end = end)
            self._parts_by_name[p.name] = p
        
    def _parse_alignment(self, alignment_file):
    
        import os
    
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
                    if len(self.partitions()) < 1:
                        p = Partition(alignment = self, name = 'all', type = 'DNA', start = 1, end = int(self.ncols))
                        self._parts_by_name[p.name] = p
                else:
                    name = toks[0]
                    seq = toks[1]
                    self._taxa.append(name)
                    for p in self.partitions():
                        p.set_seq_for_taxon(name, seq[p.start-1:p.end])
            else:
                raise IndexError("too many items on line '" + line + "' in alignment") 
    
        assert(ntax == len(self._taxa))
    
    def partitions(self): # should return an iterator over the partitions in this alignment
        return self._parts_by_name.values()
    
    def get_partition(self, p):
        return self._parts_by_name[p]
    
    def random_part(self): # should return a partition object
        return list(self._parts_by_name.keys())[random.randrange(len(self._parts_by_name))]

    def ntaxa(self):
        return len(self._taxa)

    def nparts(self):
        return len(self._parts_by_name)

    def taxon_labels(self):
        return self._taxa
        
    def partition_labels(self):
        return self._parts_by_name.keys()

class Partition():

    def __init__(self, alignment = None, name = '', type = '', start = 0, end = 0):
        self.alignment = alignment
        self.name = name
        self.type = type
        self.start = start
        self.end = end
        self._seqs_by_taxa = {}

    def set_seq_for_taxon(self, t, seq):
        self._seqs_by_taxa[t] = seq
    
    def get_seq_for_taxon(self, t):
        return self._seqs_by_taxa[t]

class _Subsampler():

    DEFAULT_LABEL = 'subsampled'

    def __init__(self, alignment):
        import random
        self.alignment = alignment
        self._init_subsample()
    
    def _init_subsample(self):

        self._sample_bitmap = {}
        self._num_taxa_sampled_by_part = {}
        self._num_parts_sampled_by_taxon = {}

        # iterate over taxa first since they are the level 0 keys in the sample bitmap
        for t in self.alignment.taxon_labels():
            self._num_parts_sampled_by_taxon[t] = 0
            self._sample_bitmap[t] = {}

        # then add the partitions
        for p in self.alignment.partition_labels():
            self._num_taxa_sampled_by_part[p] = 0
            for t in self.alignment.taxon_labels():
                self._sample_bitmap[t][p] = False
        
        self.k = 0
    
    def get_sampling_proportion(self):
        return float(self.k) / (self.alignment.nparts() * self.alignment.ntaxa())
    
    def _validate(self, t, p):
        assert t in self._sample_bitmap
        assert p in self._sample_bitmap[t]
        assert t in self._num_parts_sampled_by_taxon
        assert p in self._num_taxa_sampled_by_part
    
    def set_random_seed(self, s):
        random.seed(s)

    def sample(self, t, p):
        self._validate(t,p)
        self._sample_bitmap[t][p] = True
        self._num_taxa_sampled_by_part[p] += 1
        self._num_parts_sampled_by_taxon[t] += 1
        self.k += 1
    
    def unsample(self, t, p):
        self._validate(t,p)
        self._sample_bitmap[t][p] = False
        self._num_taxa_sampled_by_part[p] -= 1
        self._num_parts_sampled_by_taxon[t] -= 1
        self.k -= 1
    
    def is_sampled(self, t, p):
        self._validate(t,p)
        return self._sample_bitmap[t][p]
    
    def number_parts_sampled_for(self, t):
        return self._num_parts_sampled_by_taxon[t]
    
    def number_taxa_sampled_for(self, p):
        return self._num_taxa_sampled_by_part[p]
    
    # iterators over the dicts to return data sorted by sampling proportion
    def taxon_labels_sorted(self):
        return sorted(self.alignment.taxon_labels(), \
            cmp = lambda p, q: cmp(self._num_parts_sampled_by_taxon[p], \
                                   self._num_parts_sampled_by_taxon[q]), reverse=True)

    def partition_labels_sorted(self):
        return sorted(self.alignment.partition_labels(), 
            cmp = lambda p, q: cmp(self._num_taxa_sampled_by_part[p], \
                                   self._num_taxa_sampled_by_part[q]), reverse=True)

    def write_subsampled_output(self, label=''):
        self.output_label = self.alignment.base_name + '.' + \
            (self.DEFAULT_LABEL + '.' + label if len(label) > 0 else self.DEFAULT_LABEL)
        
        with open(self.output_label + '.sampling_matrix.tsv', 'w') as s:
            s.write('\t'+'\t'.join(self.partition_labels_sorted())+'\n')
            for t in self.taxon_labels_sorted():
                s.write(t)
                for p in self.partition_labels_sorted():
                    s.write('\t' + ('1' if self.is_sampled(t,p) else '0'))
                s.write('\n')

        with open(self.output_label + '.phy', 'w') as a:
            a.write('{0} {1}\n'.format(self.alignment.ntaxa(), self.alignment.ncols))
            for t in self.taxon_labels_sorted():
                a.write(t + ' ')
                for p in self.partition_labels_sorted():
                    pp = self.alignment.get_partition(p)
                    if self.is_sampled(t,p):
                        a.write(pp.get_seq_for_taxon(t))
                    else:
                        a.write('-' * (pp.end - pp.start))
                a.write('\n')

class SimpleSubsampler(_Subsampler):

    def __init__(self, alignment):
        _Subsampler.__init__(self, alignment)
    
    def subsample(self, probability_function):
        self._init_subsample()
        
        for t in self.alignment.taxon_labels():
            for p in self.alignment.partition_labels():
                if probability_function():
                    self.sample(t, p)

            if (self.number_parts_sampled_for(t) < 1):
                p = self.alignment.random_part().name
                self.sample(t, p)

class PhylogeneticSubsampler(_Subsampler):

    IMPOSSIBLY_HIGH_RATE = 1000000
    NESTED_SAMPLING_RATE = 0.01
    
    def __init__(self, alignment, tree, rates):
        _Subsampler.__init__(self, alignment)
        self.set_tree(tree)
        self.set_rates(rates)

        for p in self.alignment.partitions():
            p.sampled_nodes = set()
        
        for n in self.tree.iternodes():
            n.sampled_partitions = set()
    
    def set_tree(self, t):
        leaf_labels = set([l.label for l in t.leaves()])
        alignment_labels = set(self.alignment.taxon_labels())
        
        a = leaf_labels.difference(alignment_labels)
        b = alignment_labels.difference(leaf_labels)
        
        # might want to allow a tree with a superset of labels in the alignment
        if not len(a) == 0 and len(b) == 0:
            raise ValueError('Taxon mismatch between tree and alignment. Taxa [{}] are in ' \
                             'the tree but not in alignment, and/or [{}] are in the alignment ' \
                             'but not in the tree.'.format(', '.join(a),', '.join(b)))
        self.tree = t
    
    def set_rates(self, rates):
        rates_labels = set(rates.keys())
        parts_labels = set(self.alignment.partition_labels())
        
        a = rates_labels.difference(parts_labels)
        b = parts_labels.difference(rates_labels)
        
        # might want to allow a rate set for a superset of partitions in the alignment
        if not len(a) == 0 and len(b) == 0:
            raise ValueError('Name mismatch between specified rates and partitions. Rates are ' \
                             'specified for partitions [{}] which are not known, and/or known ' \
                             'partitions [{}] do not have specified rates.'.format(', '.join(a),', '.join(b)))
        self.rates = rates        

        lowest_rate = self.IMPOSSIBLY_HIGH_RATE
        for name, r in rates.iteritems():
            if r < lowest_rate:
                self.slowest_partition = name
                lowest_rate = r

    @staticmethod
    def has_parent_sampled_for(n, p):
#        print 'start : ' + '['+','.join([l.label for l in n.leaves()])+']'
        while not n.is_root:
            n = n.parent
#            print('['+','.join([l.label for l in n.leaves()])+']')
            if p in n.sampled_partitions:
                return True
        return False
        
    @staticmethod
    def get_sampling_frequency(n, p):
        # return the sampling frequency for this this rate at this node, based on:
        #   - the depth of of the node
        #   - the rate of the partition
        #   - the number of nodes already sampled for this partition
        
        return 0.5 # dummy

    def set_sampled_partition(self, n, p):
        self.record_sampled_node_on_partition(n, p)
        self.record_sampled_partition_on_node(n, p)
        
    def record_sampled_node_on_partition(self, n, p):
        self.alignment.get_partition(p).sampled_nodes.add(n)

    def record_sampled_partition_on_node(self, n, p):
        n.sampled_partitions.add(p)

    def subsample(self):
        
        import phylo3, random
        self._init_subsample()
        
        for n in self.tree.iternodes(phylo3.PREORDER):
        
            if n.is_tip:
                continue
                
            if n.is_root:
                self.set_sampled_partition(n, self.slowest_partition)
                continue
            
            for p in self.alignment.partition_labels():
                s = self.get_sampling_frequency(n, self.rates[p])
            
                if s > random.random():
                    if not self.has_parent_sampled_for(n, p):
                        self.set_sampled_partition(n, p)
                    else:
                        if self.NESTED_SAMPLING_RATE > random.random():
                            self.set_sampled_partition(n,p)
        
    def report_sampled_partitions(self):
        for n in self.tree.iternodes():
            print('['+','.join([l.label for l in n.leaves()])+']')
            print('    '+','.join(n.sampled_partitions))