#!/usr/bin/env python3
import argparse, gzip, os
from collections import Counter

def Command_line():
	parser = argparse.ArgumentParser(description="non_miRNA_BLAST_Parse.py sorts through tabulated BLAST results (whose reads originated from a FASTA file containing reads that did not map to the Macaca mulatta miRNA reference) and outputs the number of mapped read counts for each category of ncRNA. When a read maps equally well to >1 gene family, it will be ignored unless >80% of results map to one particular family. This script will also determine the number of reads from the ncRNA FASTA file which map to the ncRNA reference and output the percentage that map.")
	parser.add_argument('-b', default="miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane2_25605_adapter_removed_non_miRNA_BLAST_RESULTS.txt.gz", type=str, help='BLAST results text file to be analyzed for small RNA expression', metavar='BLASTResultsFile')
	parser.add_argument('-f', default="miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane2_25605_adapter_removed_non_miRNA.fasta.gz", type=str, help='FASTA file that was BLASTed against the reference', metavar='FASTAFile')
	parser.add_argument('-r', default="miRNA_analysis/Macaca_mulatta_v69_ncrna.fa.gz", type=str, help='File containing the sequence name and small RNA type for Macaca mulatta ncRNA sequences', metavar='ReferenceFile')
	parser.add_argument('-o1', default="miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane2_25605_adapter_removed_not_macaque_non_miRNA_expression.txt", type=str, help='Output file for miRNA family expression counts', metavar='ExpressionOutput')
	parser.add_argument('-o2', default="miRNA_analysis/FASTq_files_new/non_miRNA_Reads_mapped.txt", type=str, help='', metavar='ReadsMappedOutput')
	
	args = parser.parse_args()
	blast_file = args.b
	fasta_file = args.f
	ref_file = args.r
	o_file = args.o1
	o2_file = args.o2

	return(blast_file, fasta_file, ref_file, o_file, o2_file)

#blast_file, fasta_file, ref_file, o_file, o2_file = Command_line()
fasta_filelist = ["miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane2_25605_adapter_removed_non_miRNA.fasta.gz",
				  "miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane3_35088_adapter_removed_non_miRNA.fasta.gz",
				  "miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane6_23930_adapter_removed_non_miRNA.fasta.gz",
				  "miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane7_34686_adapter_removed_non_miRNA.fasta.gz",
				  "miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane8_34084_adapter_removed_non_miRNA.fasta.gz",
				  "miRNA_analysis/FASTq_files_new/100402_SOLEXA2_FC614B1AAXX_Lane4_33761WT_adapter_removed_non_miRNA.fasta.gz",
				  "miRNA_analysis/FASTq_files_new/100402_SOLEXA2_FC614B1AAXX_Lane5_33761LB_adapter_removed_non_miRNA.fasta.gz",
				  "miRNA_analysis/FASTq_files_new/100402_SOLEXA2_FC614B1AAXX_Lane6_34249WT_adapter_removed_non_miRNA.fasta.gz",
				  "miRNA_analysis/FASTq_files_new/100402_SOLEXA2_FC614B1AAXX_Lane7_34249LB_adapter_removed_non_miRNA.fasta.gz",
				  "miRNA_analysis/FASTq_files_new/100623_SOLEXA2_0623_FC61YNKAAXX_Lane8_35049_adapter_removed_non_miRNA.fasta.gz"]
blast_filelist = ["miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane2_25605_adapter_removed_non_miRNA_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane3_35088_adapter_removed_non_miRNA_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane6_23930_adapter_removed_non_miRNA_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane7_34686_adapter_removed_non_miRNA_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane8_34084_adapter_removed_non_miRNA_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100402_SOLEXA2_FC614B1AAXX_Lane4_33761WT_adapter_removed_non_miRNA_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100402_SOLEXA2_FC614B1AAXX_Lane5_33761LB_adapter_removed_non_miRNA_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100402_SOLEXA2_FC614B1AAXX_Lane6_34249WT_adapter_removed_non_miRNA_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100402_SOLEXA2_FC614B1AAXX_Lane7_34249LB_adapter_removed_non_miRNA_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100623_SOLEXA2_0623_FC61YNKAAXX_Lane8_35049_adapter_removed_non_miRNA_BLAST_RESULTS.txt.gz"]

ref_file = "miRNA_analysis/Macaca_mulatta_v69_ncrna.fa.gz"
o_file = "miRNA_analysis/FASTq_files_new/non_miRNA_ncRNA_categories.txt"
o2_file = "miRNA_analysis/FASTq_files_new/non_miRNA_reads_mapped.txt"
			
# Open Macaca_mulatta_v69_ncrna.fa, read in each ncRNA's name, and assign that name to one of 7 ncRNA categories:
# ncrna:snRNA, ncrna:snoRNA, ncrna:rRNA, ncrna:misc_RNA, mt_genbank_import:Mt_tRNA, mt_genbank_import:Mt_rRNA, ncrna:miRNA

small_rna_dict = {}
with gzip.open(ref_file, 'rb') as infile:
	for line in infile:
		line = line.decode('utf-8')
		if line[0:1] == ">":
			line = line.split(" ")
			seq_name = line[0][1:]
			small_rna_type = line[1]
			small_rna_dict[seq_name] = small_rna_type

# Begin parsing each pair of FASTA & BLAST data, outputting for each file (o_file) the number of reads mapped to each ncRNA category, and also
# outputting (o2_file) for each file the number of non-miRNA reads which also did not map to the Macaque ncRNA reference

with open(o_file, 'w') as outfile:
	outline = "Filename\tncrna:rRNA\tpct_ncrna:rRNA\tncrna:snoRNA\tpct_ncrna:snoRNA\tncrna:snRNA\tpct_ncrna:snRNA\tmt_genbank_import:Mt_tRNA\tpct_mt_genbank_import:Mt_tRNA\tmt_genbank_import:Mt_rRNA\tpct_mt_genbank_import:Mt_rRNA\tncrna:misc_RNA\tpct_ncrna:misc_RNA\tncrna:miRNA\tpct_ncrna:miRNA\tMixed_BLAST_Hits\tpct_Mixed_BLAST_Hits\n"
	outfile.write(outline)
	with open(o2_file, 'w') as outfile2:
		outline = "FASTA_Filename\tFASTA_non_miRNA_reads\tBLAST_results_reads\tReads_missing_in_BLAST_results\tPct_reads_mapped\n"
		outfile2.write(outline)
		for i in range(len(fasta_filelist)):
			fasta_file = fasta_filelist[i]
			blast_file = blast_filelist[i]
			print("Working with",str(fasta_file))
			
			# Open FASTA file containing reads which did not map to Macaca mulatta miRNA reference. Make a list of all read names, then convert list to a set
			
			fasta_read_list = []
			with gzip.open(fasta_file, 'rb') as infile:
				for line in infile:
					line = line.decode('utf-8')
					line = line.rstrip()
					if line[0:1] == ">":
						fasta_read_list.append(line[1:])
			fasta_read_set = set(fasta_read_list)
			fasta_read_list = []
			
			# Parse BLAST results file for reads that did not map to miRNA but might map to a different ncRNA.
			# Aside from noting the ncRNA category for each read, this script also keeps track of read names in a list (then a set)
			# so that the total number of reads that failed to map to a miRNA as well as all known Macaque ncRNAs can be calculated.
			
			blast_read_list = []
			rna_type_count = Counter()
			read_rna_type_count = {}
			prev_read_name = ""
			with gzip.open(blast_file, 'rb') as infile:
				for line in infile:
					line = line.decode('utf-8')
					line = line.split("\t")
					read_name = line[0]
					blast_read_list.append(read_name)
					seq_name = line[1]
					seq_small_rna_type = small_rna_dict[seq_name]
					if read_name == prev_read_name:
						read_rna_type_count[read_name][seq_small_rna_type] += 1
					else:
						read_rna_type_count[read_name] = Counter({seq_small_rna_type: 1})
						if prev_read_name != "":
							total_matches = sum([i for i in read_rna_type_count[prev_read_name].values()])
							hit_was_found = False
							for seq_small_rna_type in read_rna_type_count[prev_read_name]:
								if read_rna_type_count[prev_read_name][seq_small_rna_type] / total_matches >= 0.8:
									rna_type_count[seq_small_rna_type] += 1
									hit_was_found = True
									break
							if hit_was_found == False:
								rna_type_count["Mixed_BLAST_Hits"] += 1
							del read_rna_type_count[prev_read_name]
						prev_read_name = read_name
			
			blast_read_set = set(blast_read_list)
			blast_read_list = []
			
			# Write out ncRNA category counts for this file
			
			total_c = sum([x for x in rna_type_count.values()])
			rRNA_c = rna_type_count["ncrna:rRNA"]
			rRNA_pct = round(rRNA_c / total_c * 100, 2)
			snoRNA_c = rna_type_count["ncrna:snoRNA"]
			snoRNA_pct = round(snoRNA_c / total_c * 100, 2)
			snRNA_c = rna_type_count["ncrna:snRNA"]
			snRNA_pct = round(snRNA_c / total_c * 100, 2)
			mt_tRNA_c = rna_type_count["mt_genbank_import:Mt_tRNA"]
			mt_tRNA_pct = round(mt_tRNA_c / total_c * 100, 2)
			mt_rRNA_c = rna_type_count["mt_genbank_import:Mt_rRNA"]
			mt_rRNA_pct = round(mt_rRNA_c / total_c * 100, 2)
			miscRNA_c = rna_type_count["ncrna:misc_RNA"]
			miscRNA_pct = round(miscRNA_c / total_c * 100, 2)
			miRNA_c = rna_type_count["ncrna:miRNA"]
			miRNA_pct = round(miRNA_c / total_c * 100, 2)
			mixed_c = rna_type_count["Mixed_BLAST_Hits"]
			mixed_pct = round(mixed_c / total_c * 100, 2)
			outline = "\t".join([os.path.basename(fasta_file), str(rRNA_c), str(rRNA_pct), str(snoRNA_c), str(snoRNA_pct), str(snRNA_c), str(snRNA_pct), str(mt_tRNA_c), str(mt_tRNA_pct), str(mt_rRNA_c), str(mt_rRNA_pct), str(miscRNA_c), str(miscRNA_pct), str(miRNA_c), str(miRNA_pct), str(mixed_c), str(mixed_pct)]) + "\n"
			outfile.write(outline)
			
			# Write out % reads mapping for this file
			
			fasta_read_count = len(fasta_read_set)
			blast_read_count = len(blast_read_set)
			missing_read_set = fasta_read_set - blast_read_set
			missing_read_count = len(missing_read_set)
			pct_mapped = round((fasta_read_count - missing_read_count) / fasta_read_count * 100,2)
			
			outline = "\t".join([os.path.basename(fasta_file), str(fasta_read_count), str(blast_read_count), str(missing_read_count), str(pct_mapped)]) + "\n"
			outfile2.write(outline)