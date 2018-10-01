#!/usr/bin/python

# 03/01/2018

# This program takes the TargetP, PredAlgo and results from a manual search for PTS1 and PTS2 sequences and prints a data table for each PHYLDOG ORTHOGROUP with protein locations (1 or 0) for each subcellular location (chloroplast, mirotchondria, signalling, peroxisome or other) 
# Table Structure:  Chloroplast(1 or 0) / Mitochondria(1 or 0)  / Signalling(1 or 0)  / Peroxisome (1 or 0) / Other(1 or 0) / 
# One table and therefore one file for each orthogroup
# For the algae genes the PredAlgo result is taken unless it does not exist as the sequence is too short (na) in which case the TargetP result is taken. 
# For the PTS search: if a gene has the PTS1 or PTS2 and also a TargetP or PredAlgo localisation then the TargetP/PredAlgo one is taken. 

# The user can choose on the command line which reliability scores (RS) for TargetP are allowed. 1 indicates the strongest prediction, and 5 the lowest. So the value the user puts in is the weakest prediction score that is allowed. For example, selecting '3' will take all the results with a RS of 3, 2, and 1 as being true.   

import sys
import glob
import os
import csv
import fnmatch
import re
from Bio import SeqIO 


def make_orthogroup_proteomes_dict():
	d = {}
	os.chdir("/home/rona/Phytozome10_updated_orthogroups/Phyldog/output_full_dataset")
	for filename in glob.glob('OG00*.recon'):	
		with open(filename) as file:
			orthogroup = (filename[:len(filename)-12])
			d[orthogroup] = []
			for line in file:
				token = line.split('\t')
				if "_" in token[0]:							# Find lines corresponding to gene ids
					d[orthogroup].append(token[0])
					
	return d		
		
def make_targetP_gene_dict(RC):
	d = {}
	os.chdir("/home/rona/subcellular_localisation_prediction/TargetP/TargetP_RESULTS")
	for filename in glob.glob('out.*.txt'):
		with open(filename) as file:              
			for line in file:    		
				token = line.rstrip().split()         # Split up line into tokens, ignores whitespace.		
				if len(token) >= 7 and "_" in token[0]:  # If it is a line containing gene output...
					reliability_class = int(token[7])
					if reliability_class <= int(RC):                  # If it has the RC specified or a lower (better) one
						d[token[0]] = (token[6])      # ...and to dictionary	
					else:
						d[token[0]] = '_'     		   # else give the gene a location of other '_'		
	return d
	
def make_pred_algo_dict():
	d = {}
	os.chdir("/home/rona/subcellular_localisation_prediction/PredAlgo")
	with open("predalgoinputfile331679.fasta_pred.txt") as file:
		next(file)
		for line in file: 
			token = line.rstrip().split()
			if len(token) == 10:
				d[token[0]] = "na"
			else:
				d[token[0]] = token[4]	
	return d
	
def species_used():
	species = []
	os.chdir("/home/rona/Phytozome10_updated_orthogroups")
	with open("SpeciesUsedIDs.txt") as species_used_f:
		for line in species_used_f:
			token1 = line.rstrip().rsplit()
			species.append(token1[0][:-1])	
	return species

	
def make_PTS_dictionary(PTS1_tripeptides, PTS2_sequence, species_used):
	d = {}
	os.chdir("/home/rona/Phytozome10_updated_orthogroups/Proteomes")
	for item in species_used:
		with open("Species%s.fa" % item) as filename:
			for record in SeqIO.parse(filename, "fasta"):
				if find_PTS1(record, PTS1_tripeptides) or find_PTS2(record, PTS2_sequence):
					d[record.id] = "P"	
				else:
					d[record.id] = "_"
	return d
				
def find_PTS1(record, PTS1_tripeptides):
	if record.seq[-3:] in PTS1_tripeptides:
		return True
	if record.seq[-4:-1] in PTS1_tripeptides and record.seq.endswith("*"):
		return True
	return False
	
def find_PTS2(record, PTS2_sequence):
	PTS2_matches = re.findall(PTS2_sequence, (str(record.seq))[0:30])
	if PTS2_matches:
		return True
	else:
		return False
		
		
def create_output_file(orthogroup_genes, outfile, TargetP_dict, PredAlgo_dict, PTS_dict):
	for item in orthogroup_genes:
		if item.startswith(("11", "15", "26", "28", "43")):           # These are the IDs for algae species
			location = find_gene_location(item, 'algae', TargetP_dict, PredAlgo_dict, PTS_dict)
		else:
			location = find_gene_location(item, 'vascular', TargetP_dict, PredAlgo_dict, PTS_dict)
			
		write_to_file(item, location, outfile)		
			
def find_gene_location(gene, plant, TargetP_dict, PredAlgo_dict, PTS_dict):
	if plant == 'algae':
		d = PredAlgo_dict
	if plant == 'vascular':
		d = TargetP_dict
	if d[gene] == 'na':
		d = TargetP_dict
		
	if PTS_dict[gene] == 'P':
		if d[gene] == 'C':
			return('CP')
		if d[gene] == 'M':
			return('MP')
		if d[gene] == 'S' or d[gene] == 'SP':
			return('SP')
		return('P')
			
	if d[gene] == 'C':
		return('C')
	if d[gene] == 'M':
		return('M')
	if d[gene] == 'S' or d[gene] == 'SP':
		return('S')

	return('O')
							
def write_to_file(gene, location, file):
	if location == 'C':				
		file.write(str(gene).ljust(10) + repr(1).rjust(2) + repr(0).rjust(2) + repr(0).rjust(2)+ repr(0).rjust(2) + repr(0).rjust(2) + '\n')
	if location == 'CP':
		file.write(str(gene).ljust(10) + repr(1).rjust(2) + repr(0).rjust(2) + repr(0).rjust(2)+ repr(1).rjust(2) + repr(0).rjust(2) + '\n')	
	if location == 'M':
		file.write(str(gene).ljust(10) + repr(0).rjust(2) + repr(1).rjust(2) + repr(0).rjust(2)+ repr(0).rjust(2) + repr(0).rjust(2) + '\n')
	if location == 'MP':
		file.write(str(gene).ljust(10) + repr(0).rjust(2) + repr(1).rjust(2) + repr(0).rjust(2)+ repr(1).rjust(2) + repr(0).rjust(2) + '\n')		
	if location == 'S':
		file.write(str(gene).ljust(10) + repr(0).rjust(2) + repr(0).rjust(2) + repr(1).rjust(2)+ repr(0).rjust(2) + repr(0).rjust(2) + '\n')
	if location == 'SP':
		file.write(str(gene).ljust(10) + repr(0).rjust(2) + repr(0).rjust(2) + repr(1).rjust(2)+ repr(1).rjust(2) + repr(0).rjust(2) + '\n')	
	if location == 'P':
		file.write(str(gene).ljust(10) + repr(0).rjust(2) + repr(0).rjust(2) + repr(0).rjust(2)+ repr(1).rjust(2) + repr(0).rjust(2) + '\n')	
	if location == 'O':
		file.write(str(gene).ljust(10) + repr(0).rjust(2) + repr(0).rjust(2) + repr(0).rjust(2)+ repr(0).rjust(2) + repr(1).rjust(2) + '\n')
	if location not in ['C', 'CP', 'M', 'MP', 'S', 'SP', 'P', 'O']:
		file.write('ERROR')	
							
def main():								
	print("Usage:")
	print("  script.py       TargetP_reliability_coefficient_cutoff (the score taken as accepted)")
	if len(sys.argv) < 2:
		sys.exit()	
		
	# Make a dictionary that holds every gene(item) that belongs in each orthogroup (key)					
	orthogroup_proteome_dict = make_orthogroup_proteomes_dict()                 # Fill the dictionary 

	# Make a dictionary of the TargetP results for each gene	
	TargetP_dict = make_targetP_gene_dict(sys.argv[1])	
			
	# Make a dictionary of the PredAlgo results for each algae gene
	PredAlgo_dict = make_pred_algo_dict()		
			
	# Get a list of species including in Phyldog to narrow down PTS1 search
	species_list = species_used()
	
	# Make a dictionary of the PTS1 results for each gene
	PTS1_tripeptides = ["SRL", "SRM", "SRI", "ARL", "ARM", "PRL", "SKL", "SKM", "AKL"]
	PTS2_sequence = 'R[L|I]\w\w\w\w\wHL'
	PTS_dict = make_PTS_dictionary(PTS1_tripeptides, PTS2_sequence, species_list)

	# Find final location of each gene
	os.chdir('/home/rona/subcellular_localisation_prediction/TargetP_PredAlgo_PTS1_PTS2/results_by_orthogroup/RC_%s' % sys.argv[1])
	for key in orthogroup_proteome_dict: 
		with open('%s.location.data' % key, "w") as outfile:
			create_output_file(orthogroup_proteome_dict[key], outfile, TargetP_dict, PredAlgo_dict, PTS_dict)	
			
if __name__ == '__main__':
	main()		
		