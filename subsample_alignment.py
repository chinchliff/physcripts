#! python
#
#  subsample_alignment.py
#
#  Created by Cody on 11/18/11.


partitions = {"psbBTNH": [13965, 16263],
			  "rps4": [24160, 25260],
			  "matK": [7570, 9070],
			  "rpoC2": [17678, 21561],
			  "rps16": [21562, 22609],
			  "atpB": [6083, 7569],
			  "ndhF": [11886, 13964],
			  "rbcL": [16264, 17677]}

infile = "/Users/cody/Desktop/Angiosperms/first_nexus_from_joseph/original_alignment/angiosperms.fasta"
outfile = "/Users/cody/Desktop/Angiosperms/first_nexus_from_joseph/subsampled_alignment/angiosperms.8cploci.fasta" 
partfile = "/Users/cody/Desktop/Angiosperms/first_nexus_from_joseph/original_alignment/angiosperms.8cploci.part"

alignment_in = open(infile,"rb")
alignment_out = open(outfile,"wb")

part_file = open(partfile,"wb")

part_names = partitions.keys()
part_names.sort()
#print part_names

next_position = 1
for part in part_names:
	start_original = partitions[part][0]
	end_original = partitions[part][1]
	
	difference = end_original - start_original
	
	start_new = next_position
	end_new = next_position + difference
	
	next_position = end_new + 1
	
	part_file.write("DNA, " + part + " = " + str(start_new) + "-" + str(end_new) + "\n")

part_file.close()

for line in alignment_in:
	lineTokens = line.split()
	newline = lineTokens[0] + " "

	raw_data = lineTokens[1]
	for part in part_names:
		start = partitions[part][0]
		end = partitions[part][1]
		newline += raw_data[start - 1:end]

	alignment_out.write(newline + "\n")