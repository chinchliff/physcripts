from Bio import SeqIO
from Bio import AlignIO
from Bio import Alphabet #import generic_dna
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class consensus:

    """
    The consenus class accepts a BioPython alignment object and calculates a consensus sequence
    with robust ambiguity coding to represent polymorphic columns. The resulting sequence
    is stored in the self.consensus_seq property variable.

    Properties:

        self.alignment          : the alignment to be consensed
        self.consensus_seq      : the calculated consensus sequence
        self.ncols_polymorphic  : the number of polymorphic columns encountered (does not
                                  include columns containing an N or X but no more than one
                                  other nucleotide coding)

    Methods:

        __init__        : calculates a consensus from the provided alignment and populates
                          all property variables

    """

    def __init__(self, alignment):
        
        self.alignment = alignment
        self.ncols_polymorphic = 0

        class base:
            def __init__(self, name, codes):
                self.name = name
                self.codes = list(codes)
                self.ispresent = False

        # set up ambiguity codes and base objects
        nA = base('Adenine',['A','M','R','W','V','H','D'])
        nC = base('Cytosine',['C','M','S','Y','V','H','B'])
        nG = base('Guanine',['G','R','S','K','V','B','D'])
        nT = base('Thymine',['T','W','Y','K','H','B','D'])

        bases = [nA,nC,nG,nT]

        # this will hold the consensus sequence
        conseq = str()

        for column in range(self.alignment.get_alignment_length()):

            # for each colum in the alignment, make a list of all the codes present
            allchars = [c for c in (self.alignment.get_column(column))]
            clist = [c for c in allchars if c != allchars[0]] + [allchars[0]]

            # if we only have one code in this column, nothing is more is necessary. save it and we're done
            if len(clist) == 1:
                acode = clist[0]

            else:
                # if we have multiple codes, first we strip out meaningless 'N' and 'X' chars
                symbols = [symbol for symbol in clist if symbol not in ('N','X')]
                nsymbols_incl_gaps = len(symbols)

                # strip gap chars. at some point this should get the appropriate gap codings from the alphabet
                symbols_excl_gaps = [symbol for symbol in symbols if symbol != '-']
                nsymbols_excl_gaps = len(symbols_excl_gaps)

                if nsymbols_excl_gaps == 1:
                    # if we found a gap char, increment ncols_polymorphic, even though the rest of the codes
                    # in this column are identical. (in this case, the polymorphism indicates an indel).
                    if nsymbols_incl_gaps > nsymbols_excl_gaps:
                        self.ncols_polymorphic += 1
                    
                    acode = symbols_excl_gaps[0]
                    
                else: # we still have multiple codes, so prepare to determine what bases they represent
                    nA.ispresent = False
                    nC.ispresent = False
                    nG.ispresent = False
                    nT.ispresent = False

                    # increment polymorphism count
                    self.ncols_polymorphic += 1

                    # for each base, see if any of our codes corresponding to it.
                    # if we find a match, mark that base as present and move on to the next base
                    for base in bases:
                        for c in base.codes:
                            if c in symbols:
                                base.ispresent = True
                                break

                    A = nA.ispresent
                    C = nC.ispresent
                    G = nG.ispresent
                    T = nT.ispresent

					# The following ambiguity coding is used
##                  M = A, C
##                  R = A, G
##                  W = A, T
##                  S = C, G
##                  Y = C, T
##                  K = G, T
##                  V = A, C, G
##                  H = A, C, T
##                  D = A, G, T
##                  B = C, G, T

                    # record the appropriate ambiguity code for the combination of bases we found 
                    if A:
                        if C:
                            if G:
                                if T:
                                    acode = 'N'
                                else:
                                    acode = 'V'
                            elif T:
                                acode = 'H'
                            else:
                                acode = 'M'
                        elif G:
                            if T:
                                acode = 'D'
                            else:
                                acode = 'R'
                        elif T:
                            acode = 'W'
                    elif C:
                        if G:
                            if T:
                                acode = 'B'
                            else:
                                acode = 'S'
                        elif T:
                            acode = 'Y'
                    else:
                         acode = 'K'

            # we're done with this column
            conseq += acode

        # all columns consensed, now save consensus seq and exit
        self.consensus_seq = conseq
        return None
