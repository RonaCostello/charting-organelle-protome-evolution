#!/usr/bin/python

# 03/01/2018

# This program takes the TargetP, PredAlgo and manual PTS (PTS1 & PTS2) search results and prints a data table for each species used with protein locations (1 or 0) for each subcellular location (chloroplast, mirotchondria, signalling, peroxisome or other) 
# Table Structure:  Chloroplast(1 or 0) / Mitochondria(1 or 0)  / Signalling(1 or 0)  / Peroxisome (1 or 0) / Other(1 or 0) / 
# One table and therefore one file for each species 
# For the algae genes the PredAlgo result is taken unless it does not exist as the sequence is too short (na) in which case the TargetP result is taken. 
# For the peroxisomal protein search: if a gene has PTS1 or PTS2 then it is labelled as peroxisomal. If it also has a TargetP or PredAlgo localisation then it is labelled as both e.g. CP, MP, SP

# The user can choose on the command line which reliability scores (RS) for TargetP are allowed. 1 indicates the strongest prediction, and 5 the lowest. So the value the user puts in is the weakest prediction score that is allowed. For example, selecting '3' will take all the results with a RS of 3, 2, and 1 as being true.   

import sys
import glob
import os
import csv
import fnmatch
import re
from Bio import SeqIO 


def make_species_proteomes_dict():
	d = {}
	species = []
	os.chdir("/home/rona/Phytozome10")
	with open("SpeciesUsedIDs.txt") as species_used_f:
		for line in species_used_f:
			token1 = line.rstrip().rsplit()
			species.append(token1[0][:-1])	
			
	os.chdir("/home/rona/Phytozome10/Proteomes")
	for item in species:
		d[item] = []
		with open("Species%s.fa" % item) as f:
			for line in f:
				if line.startswith(">"):
					token2 = line.rstrip().rsplit()
					d[item].append(token2[0][1:])
	return d, species		

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

def make_PTS_dictionary(PTS1_tripeptides, PTS2_sequence, species_used):
	d = {}
	os.chdir("/home/rona/Phytozome10/Proteomes")
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

def output_location_files(species_genes, outfile, TargetP_dict, PredAlgo_dict, PTS_dict):
	for item in species_genes:
		if item.startswith(("11", "15", "26", "28", "43")):           # These are the IDs for algae species
			location = find_gene_location(item, 'algae', TargetP_dict, PredAlgo_dict, PTS_dict)
		else:
			location = find_gene_location(item, 'vascular', TargetP_dict, PredAlgo_dict, PTS_dict)
			
		result = [item, location]
		outfile.writerow(result)		
			
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
	else:
		return('O')
							
def main():								
	print("Usage:")
	print("  script.py       TargetP_reliability_coefficient_cutoff (the score taken as accepted)")
	if len(sys.argv) < 2:
		sys.exit()	
		
	# Make a dictionary that holds every gene(item) that belongs in each species (key)					
	species_proteome_dict, species_list = make_species_proteomes_dict()                 # Fill the dictionary 

	# Make a dictionary of the TargetP results for each gene	
	TargetP_dict = make_targetP_gene_dict(sys.argv[1])	
			
	# Make a dictionary of the PredAlgo results for each algae gene
	PredAlgo_dict = make_pred_algo_dict()		
			
	# Make a dictionary of the PTS1 results for each gene
	PTS1_tripeptides = ["SRL", "SRM", "SRI", "ARL", "ARM", "PRL", "SKL", "SKM", "AKL"]
	PTS2_sequence = 'R[L|I]\w\w\w\w\wHL'
	PTS_dict = make_PTS_dictionary(PTS1_tripeptides, PTS2_sequence, species_list)

	# Find final location of each gene
	os.chdir('/home/rona/subcellular_localisation_prediction/TargetP_PredAlgo_PTS1_PTS2/results_by_species/RC_%s' % sys.argv[1])
	for key in species_proteome_dict: 
		with open('%s.location.data' % key, "w") as outfile:
			w = csv.writer(outfile, delimiter = ',', quotechar = ' ', quoting = csv.QUOTE_MINIMAL)
			output_location_files(species_proteome_dict[key], w, TargetP_dict, PredAlgo_dict, PTS_dict)	
			
if __name__ == '__main__':
	main()		
		
