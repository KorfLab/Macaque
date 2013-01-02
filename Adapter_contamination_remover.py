#!/usr/bin/env python3
import argparse, re, gzip

def Command_line():
	parser = argparse.ArgumentParser(description="Adapter_contamination_remover.py removes reads from a FASTQ file which are \
									 solely Illumina adapter or PCR primer contamination due to short read sequencing. This \
									 script takes in a file listing kmers which, according to FastQC, are likely to be Illumina \
									 adapter sequences or PCR primer sequences.")
	parser.add_argument('-i', default="New_miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane2_25605.fastq.gz", type=str, help='FASTQ file to be adapter trimmed.', metavar='FASTQFile')
	parser.add_argument('-k', default="New_miRNA_analysis/FASTq_files_new/FastQC_Illumina_Adapter_Kmers.txt", type=str, help='File listing adapter kmers found by FastQC', metavar='AdapterKmers')
	parser.add_argument('-o', default="New_miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane2_25605_adapters_primers_removed_trimmed.fastq.gz", type=str, help='Output filepath for FASTQ file with adapter-containing reads removed', metavar='OutputPath')
	parser.add_argument('-fivep', default=3, type=int, help='Number of bps to cut off the 5\' end.')
	parser.add_argument('-threep', default=24, type=int, help='Position which marks the first bp to be removed on the 3\' end.')
	
	args = parser.parse_args()
	fasta_file = args.i
	kmer_file = args.k
	five_prime_cutoff = args.fivep
	three_prime_cutoff = args.threep
	outpath = args.o
	
	return(fasta_file, kmer_file, five_prime_cutoff, three_prime_cutoff, outpath)

fasta_file, kmer_file, five_prime_cutoff, three_prime_cutoff, outpath = Command_line()

# Prepare kmer list
kmer_list = []
with open(kmer_file) as infile:
	for line in infile:
		line = line.split("\t")
		if line[0] != "Sequence" and line[0].strip() != "":
			kmer_list.append(line[0].strip())
kmer_set = set(kmer_list) # Ensures each kmer is unique. Adapter kmers found in multiple files have all but one instance removed

# Open input and output files
gzipped_in = False
gzipped_out = False
if fasta_file[-3:] == ".gz":
	infile = gzip.open(fasta_file, 'rb')
	gzipped_in = True
else:
	infile = open(fasta_file)

if outpath[-3:] == ".gz":
	outfile = gzip.open(outpath, 'wb')
	gzipped_out = True
else:
	outfile = open(outpath, 'w')

# Parse infile and output reads which do not contain adapter contamination
prev_line = []
output_read = False
for line in infile:
	if gzipped_in == True: line = str(line, encoding='utf8')
	splitline = line.split("\t")
	if splitline[0][:1] == "@":
		# First line for a read
		prev_line = [line]
		output_read = False
	elif splitline[0][:1] == "+" and output_read == True:
		# Third line for a read
		# Only store if second line did not match a contamination kmer
		# and did not contain an ambiguous 'N' nucleotide.
		prev_line.append(line)
	elif len(prev_line) == 3:
		# We have reached the read quality line for read we intend on keeping.
		# Make sure to perform 5' and 3' trimming here as well.
		# Then, output all four lines
		prev_line.append(line[five_prime_cutoff:three_prime_cutoff-1] + "\n")
		if gzipped_out == True:
			for data in prev_line:
				outfile.write(bytes(data,"UTF-8"))
		else:
			for data in prev_line:
				outfile.write(data)
	else:
		stripped_splitline = splitline[0].strip()
		if len(prev_line) == 1 and stripped_splitline not in kmer_set and 'N' not in stripped_splitline:
			# Second line for a read. If it does not match a contamination
			# kmer and does not contain an ambiguous 'N' nucleotide, proceed
			# to write out all four read lines. Note that 5' and 3' trimming
			# occurs in this step.
			output_read = True
			prev_line.append(line[five_prime_cutoff:three_prime_cutoff-1] + "\n")
infile.close()
outfile.close()