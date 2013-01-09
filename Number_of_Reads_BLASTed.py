#!/usr/bin/env python3
import argparse, gzip, os

def Command_line():
	parser = argparse.ArgumentParser(description="Calculates the number of reads from a FASTQ file that are present in a BLAST results file.")
	parser.add_argument('-f', default="miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane2_25605_adapter_removed.fastq.gz", type=str, help='FASTQ filename.', metavar='FASTQ')
	parser.add_argument('-b', default="miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane2_25605_adapter_removed_BLAST_RESULTS.txt.gz", type=str, help='BLAST ', metavar='BLASTFile')
	parser.add_argument('-o', default="miRNA_analysis/FASTq_files_new/", type=str, help='Output directory', metavar='OutputDir')
	
	args = parser.parse_args()
	fastq_f = args.f
	blast_f = args.b
	outdir = args.o
	out_f = os.path.splitext(blast_f)[0] + "_reads_mapped.txt"
	outdir = os.path.join(outdir,)
	
	return(fastq_f, blast_f, out_f)

#fastq_f, blast_f, out_f = Command_line()
fastq_filelist = ["miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane2_25605_adapter_removed.fastq.gz",
				  "miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane3_35088_adapter_removed.fastq.gz",
				  "miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane6_23930_adapter_removed.fastq.gz",
				  "miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane7_34686_adapter_removed.fastq.gz",
				  "miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane8_34084_adapter_removed.fastq.gz",
				  "miRNA_analysis/FASTq_files_new/100402_SOLEXA2_FC614B1AAXX_Lane4_33761WT_adapter_removed.fastq.gz",
				  "miRNA_analysis/FASTq_files_new/100402_SOLEXA2_FC614B1AAXX_Lane5_33761LB_adapter_removed.fastq.gz",
				  "miRNA_analysis/FASTq_files_new/100402_SOLEXA2_FC614B1AAXX_Lane6_34249WT_adapter_removed.fastq.gz",
				  "miRNA_analysis/FASTq_files_new/100402_SOLEXA2_FC614B1AAXX_Lane7_34249LB_adapter_removed.fastq.gz",
				  "miRNA_analysis/FASTq_files_new/100623_SOLEXA2_0623_FC61YNKAAXX_Lane8_35049_adapter_removed.fastq.gz"]
blast_filelist = ["miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane2_25605_adapter_removed_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane3_35088_adapter_removed_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane6_23930_adapter_removed_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane7_34686_adapter_removed_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100302_SOLEXA2_FC6147KAAXX_Lane8_34084_adapter_removed_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100402_SOLEXA2_FC614B1AAXX_Lane4_33761WT_adapter_removed_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100402_SOLEXA2_FC614B1AAXX_Lane5_33761LB_adapter_removed_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100402_SOLEXA2_FC614B1AAXX_Lane6_34249WT_adapter_removed_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100402_SOLEXA2_FC614B1AAXX_Lane7_34249LB_adapter_removed_BLAST_RESULTS.txt.gz",
				  "miRNA_analysis/FASTq_files_new/100623_SOLEXA2_0623_FC61YNKAAXX_Lane8_35049_adapter_removed_BLAST_RESULTS.txt.gz"]
out_f = "miRNA_analysis/FASTq_files_new/Reads_mapped.txt"
with open(out_f, 'w') as outfile:
	outline = "FASTQ_Filename\tFASTQ_adapter_removed_reads\tBLAST_results_reads\tReads_missing_in_BLAST_results\tPct_reads_mapped\n"
	outfile.write(outline)
	for i in range(len(fastq_filelist)):
		fastq_f = fastq_filelist[i]
		blast_f = blast_filelist[i]
		fastq_list = []
		blast_list = []
		print("Working with",str(fastq_f))
		
		if fastq_f[-3:] == ".gz":
			infile = gzip.open(fastq_f, 'rb')
		else:
			infile = open(fastq_f)
		for line in infile:
			if fastq_f[-3:] == ".gz": line = str(line, encoding='utf8')
			line = line.strip()
			if line[:1] == "@":
				fastq_list.append(line[1:])
		infile.close()
		fastq_set = set(fastq_list)
		fastq_list = []
		
		if blast_f[-3:] == ".gz":
			infile = gzip.open(blast_f, 'rb')
		else:
			infile = open(blast_f)
		for line in infile:
			if blast_f[-3:] == ".gz": line = str(line, encoding='utf8')
			line = line.strip()
			blast_list.append(line[0])
		infile.close()
		blast_set = set(blast_list)
		blast_list = []
		
		fastq_read_count = len(fastq_set)
		blast_read_count = len(blast_set)
		missing_read_set = fastq_set - blast_set
		missing_read_count = len(missing_read_set)
		outline = str(os.path.basename(fastq_f)) + "\t" + str(fastq_read_count) + "\t" + str(blast_read_count) + "\t" + str(missing_read_count) + "\t" + str(round((fastq_read_count - missing_read_count) / fastq_read_count * 100,2)) + "\n"
		outfile.write(outline)