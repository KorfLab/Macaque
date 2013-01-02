#!/usr/bin/env python3
import argparse, re, gzip
from collections import Counter

def Command_line():
	parser = argparse.ArgumentParser(description="miRNA_BLAST_Parse.py sorts through tabulated BLAST results and outputs the number of mapped read counts for each miRNA family. When a read maps equally well to >1 gene family, it will be ignored unless >80% of results map to one particular family.")
	parser.add_argument('-i', default="miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane2_25605_adapter_removed_BLAST_results.txt.gz", type=str, help='BLAST results text file to be analyzed for miRNA expression', metavar='BLASTResultsFile')
	parser.add_argument('-o', default="miRNA_analysis/Differential_expression/100302_SOLEXA2_FC6147KAAXX_Lane2_25605_adapter_removed_expression.txt", type=str, help='Output file for miRNA family expression counts', metavar='OutputFile')
	
	args = parser.parse_args()
	i_file = args.i
	o_file = args.o

	return(i_file, o_file)

i_file, o_file = Command_line()
#pattern = r"^(\S+)\s\w+-([a-zA-Z]+-?\d+)[a-zA-Z]?-?[35]?p?"
#recomp = re.compile(pattern)
#pattern2 = r"([a-zA-Z]+)(\d+)"
#recomp2 = re.compile(pattern2)
pattern = r"^\w+-(.*)$"
recomp = re.compile(pattern)

gene_family_count = Counter()
read_gene_family_count = {}
prev_read_name = ""
with gzip.open(i_file, 'rb') as infile:
	for line in infile:
		line = str(line, encoding='utf8')
		line = line.split("\t")
		read_name = line[0]
		match = recomp.match(line[1])
		gene_family = match.group(1)
		#match = recomp.match(line)
		#if match:
		#	read_name, gene_family = match.groups()
		#	if "-" not in gene_family:
		#		match2 = recomp2.match(gene_family)
		#		first_part, second_part = match2.groups()
		#		gene_family = str(first_part) + "-" + str(second_part)
		if read_name == prev_read_name:
			read_gene_family_count[read_name][gene_family] += 1
		else:
			read_gene_family_count[read_name] = Counter({gene_family: 1})
			if prev_read_name != "":
				total_matches = sum([i for i in read_gene_family_count[prev_read_name].values()])
				for gene_family in read_gene_family_count[prev_read_name]:
					if read_gene_family_count[prev_read_name][gene_family] / total_matches >= 0.8:
						gene_family_count[gene_family] += 1
						break
				del read_gene_family_count[prev_read_name]
			prev_read_name = read_name
	
with open(o_file, 'w') as outfile:
	outline = "Gene_family\tReads_mapped\n"
	outfile.write(outline)
	for gene_family, count in sorted(gene_family_count.items(), key = lambda gene_family_count: gene_family_count[1], reverse=True):
		outline = str(gene_family) + "\t" + str(count) + "\n"
		outfile.write(outline)