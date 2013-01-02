#!/usr/bin/env python3
# Written by Matthew Porter @ UC Davis, Korf and Comai Labs
import argparse, re, os
from collections import Counter

def Command_line():
    parser = argparse.ArgumentParser(description="Merge_Expression_Files.py combines the reads mapped per gene for all data files in the provided folder into a single table.")
    parser.add_argument('-i', default="miRNA_analysis/Differential_expression/BLAST_against_Macaca_only", type=str, help='Input folder for miRNA expression data files', metavar='Input')
    parser.add_argument('-o', default="miRNA_analysis/Differential_expression/BLAST_against_Macaca_only/miRNA_Expression.txt", type=str, help='Output file for combined miRNA expression data files.', metavar='Output')
    
    args = parser.parse_args()
    indir = args.i
    outfile = args.o
    
    return(indir, outfile)

indir, outfile = Command_line()
pattern = r"(^.*)_adapter_removed_expression.txt$"
recomp = re.compile(pattern)

gene_family = {}
lane_list = []
for dirname, dirnames, filenames in os.walk(indir):
    for f_name in filenames:
        match = recomp.match(f_name)
        if match:
            lane = match.group(1)
            lane_list.append(lane)
            
for lane in lane_list:
    open_file_name = os.path.join(indir, str(lane) + "_adapter_removed_expression.txt")
    with open(open_file_name) as infile:
        for line in infile:
            line = line.strip()
            line = line.split("\t")
            family = line[0]
            count = line[1]
            if line[0] != "Gene_family":
                try:
                    gene_family[family][lane] = str(line[1])
                except:
                    gene_family[family] = {x: str(0) for x in lane_list}
                    gene_family[family][lane] = str(count)

with open(outfile, 'w') as outfile:
    lane_outline = "Gene_family\t" + "\t".join(lane_list) + "\n"
    outfile.write(lane_outline)
    for family in sorted(gene_family):
        outline = str(family) + "\t"
        for lane in lane_list:
            outline += gene_family[family][lane] + "\t"
        outline += "\n"
        outfile.write(outline)