#!/usr/bin/env python3
import re, copy
import pdb

def list_species():
	pattern = r">\S+\s\S+\s(\S+\s\S+)\s(\S{3}-\d+)"
	recomp = re.compile(pattern)
	species_list = []
	with open("miRNA_analysis/mirBase_miRNAs/mature.fa") as infile:
		for line in infile:
			match = recomp.match(line)
			if match:
				species, miRNA = match.groups()
				species_list.append(species)
	
	species_set = set(species_list)
	for entry in sorted(species_set):
		print(entry)

def filter_database():
	pattern = r">\S+\s\S+\s(\S+\s\S+)\s\S+"
	recomp = re.compile(pattern)
	
	seq_dict = {}
	species = ""
	name_line = ""
	with open("miRNA_analysis/mirBase_miRNAs/mature.fa") as infile:
		for line in infile:
			line = line.strip()
			match = recomp.match(line)
			if match:
				species = match.group(1)
				name_line = line
			else:
				seq = line.strip()
				try:
					seq_dict[seq][species] = name_line
				except:
					seq_dict[seq] = {species: name_line}
	
	species_order = []
	with open("miRNA_analysis/mirBase_miRNAs/Filter_species_MacaqueOnly.txt") as infile:
		for line in infile:
			species_order.append(line.strip())
	
	remove_seq_list = []
	for species in species_order:
		for seq, species_dict in seq_dict.items():
			# For when only Macaca mulatta miRNAs may be included in the output
			if species not in species_dict:
				remove_seq_list.append(seq)
			elif species in species_dict and len(species_dict.keys()) > 1:
				species_dict_copy = copy.deepcopy(species_dict)
				for delete_species in species_dict.keys():
					if species not in delete_species:
						del species_dict_copy[delete_species]
				seq_dict[seq] = species_dict_copy
			
			# For when many species may be included in the output, not just Macaca mulatta
			#if species in species_dict and len(species_dict.keys()) > 1:
			#	species_dict_copy = copy.deepcopy(species_dict)
			#	for delete_species in species_dict.keys():
			#		if species not in delete_species:
			#			del species_dict_copy[delete_species]
			#	seq_dict[seq] = species_dict_copy
	
	for seq in remove_seq_list:
		del seq_dict[seq]
	
	with open("miRNA_analysis/mirBase_miRNAs/mature_MacaqueOnly.fa", 'w') as outfile:
		for seq, species_dict in seq_dict.items():
			for species, name_line in species_dict.items():
				outline = "\n".join([name_line,seq]) + "\n"
				outfile.write(outline)

#list_species()
filter_database()