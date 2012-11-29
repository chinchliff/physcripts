from Bio import AlignIO

ifile = "comb.nex"
iformat = "nexus"

ofile = "comb.fasta"
oformat = "fasta"

ihandle = open(ifile,"rU")
ohandle = open(ofile,"w")

alignment = AlignIO.parse(ihandle,iformat)

AlignIO.write(alignment,ohandle,oformat)

#print "Alignment length %i" % alignment.get_alignment_length()

#for record in alignment :
#    print record.seq, record.id
