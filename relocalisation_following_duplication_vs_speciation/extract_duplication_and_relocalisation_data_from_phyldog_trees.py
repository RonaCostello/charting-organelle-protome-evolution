#!/usr/bin/env python

# A program that produces a csv file with each row listing a node in the phyldog orthogroup trees along with data ascertaining whether a relocalisation event happened at that node 
# (and if so, where and in which direction, and the retention score) and whether the node is a duplication, terminal duplication, or speciation node and if so what the duplciation score is. All one file. 

import glob
import csv
from ete3 import Tree
import os
import re

orthogroup_tree_recon_location = ("/home/rona/Phytozome10/Phyldog/output_full_dataset")
species_tree_location = ("/cellar/rona/Phytozome10/Phyldog")
output_location = ("/home/rona/ancestral_state_estimation/relocalisations_duplications/TargetP_PredAlgo_PTS1_PTS2_likelihoods/RC_5")                            #/cellar/rona/relocalisation_following_duplication
result_file = ("Phyldog_relocalisation_duplication_data.csv")
location_of_localisation_predictions = ("/home/rona/subcellular_localisation_prediction/TargetP_PredAlgo_PTS1_PTS2/results_by_orthogroup/RC_5")
location_of_ACE_likelihoods =("/home/rona/ancestral_state_estimation/ancestral_state_likelihoods/TargetP_PredAlgo_PTS1_PTS2_likelihoods/RC_5")

def link_nodes(recon_file):
	d = {}
	with open(recon_file) as file:
		for line in file:
			token = line.rstrip().split()
			if "n" in token[0]:	              # If the line corresponds to a node and not to a gene id
				d[token[0]] = token[1]
	return d

def get_protein_locations(orthogroup):
	d = {}
	with open("%s.location.data" %(orthogroup)) as file:
		for line in file:
			token = line.rstrip().split()
			d[token[0]] = (token[1], token[2], token[3], token[4])  # Dict: geneID = chloro(1/0), mito(1/0), signal(1/0), perox(1/0) 	
	return d
	
def get_ancestral_likelihoods(orthogroup):
	d = {}
	with open("%s.TargetP_PredAlgo_peroxisomal_likelihoods.csv" %(orthogroup)) as csvfile: 
		reader = csv.DictReader(csvfile)
		for row in reader:
			d[row['node']] = (row['likelihoods_chloroplast'], row['likelihoods_mitochondria'], row['likelihoods_secretory'], row['likelihoods_peroxisomal'])
	return d
		
def get_relocalisations(node, node_up, node_ancestral_state_dict, protein_localisation_dict):    # This function fills the list of relocalisations to be returned.
	ancestral_state_changes	= []
	for i in range(4):
		if node_ancestral_state_dict[node.name][i] == 'na' or node_ancestral_state_dict[node.name][i] == 'NA':
			ancestral_state_changes.extend(["none", "na", "na"])
			continue
		if float(node_ancestral_state_dict[node.name][i]) > 0.5 and float(node_ancestral_state_dict[node_up.name][i]) <= 0.5:  # if there is a gain
			ancestral_state_changes.extend(get_retention_score(node, node_up, i, "gain", node_ancestral_state_dict, protein_localisation_dict))
			continue
		if float(node_ancestral_state_dict[node.name][i]) <= 0.5 and float(node_ancestral_state_dict[node_up.name][i]) > 0.5:  # if there is a loss
			ancestral_state_changes.extend(get_retention_score(node, node_up, i, "loss", node_ancestral_state_dict, protein_localisation_dict))
			continue
		else:
			ancestral_state_changes.extend(["none", "na", "na"])
	return ancestral_state_changes
	
def get_retention_score(node, parent_node, subcellular_location, direction, node_ancestral_state_dict, protein_localisation_dict):
	child_node_0 = re.findall("'([^']*)'", str(parent_node.children))[0]
	child_node_1 = re.findall("'([^']*)'", str(parent_node.children))[1]
	outgroup_node = parent_node.children[0] if child_node_0 != node.name else parent_node.children[1]
	node_location = (node_ancestral_state_dict[node.name])[subcellular_location]
	outgroup_node_location = (node_ancestral_state_dict[parent_node.name])[subcellular_location]   # Assume the location of the parent node for the outgroup (prevents against changes in location due to subsequent node changes)
	leaf_retention_node = leaf_retention(node.name, node, node_location, subcellular_location, protein_localisation_dict)
	leaf_retention_outgroup = leaf_retention(outgroup_node.name, outgroup_node, outgroup_node_location, subcellular_location, protein_localisation_dict)
	return [direction, leaf_retention_node, leaf_retention_outgroup]

def leaf_retention(node_name, node_lineage, node_localisation, subcellular_location, protein_localisation_dict):
	retained_state = 0.0 if float(node_localisation) <= 0.5 else 1.0
	total_leaves, retention_in_leaves = 0.0, 0.0
	for leaf in node_lineage.iter_leaves():
		total_leaves += 1.0
		if float(protein_localisation_dict[(str(leaf))[3:]][subcellular_location]) == retained_state:
			retention_in_leaves += 1.0
	return(retention_in_leaves/total_leaves)
	
def get_duplication_score(node, species_tree_node_name, species_tree, orthogroup_tree):
	species_tree_descendants = get_descending_species(species_tree_node_name, species_tree)
	child_nodes = re.findall("'([^']*)'", (str((node).children))) 
	child_node1_descendants = get_descending_species(child_nodes[0], orthogroup_tree)
	child_node2_descendants = get_descending_species(child_nodes[1], orthogroup_tree)
	if any(i in child_node1_descendants for i in child_node2_descendants) == True:
		type = "Duplication" if len(species_tree_descendants) > 1 else "Terminal_Duplication"
		retention_child1 = percentage_retention(child_node1_descendants, species_tree_descendants)
		retention_child2 = percentage_retention(child_node2_descendants, species_tree_descendants)		
		retention_score = min(retention_child1, retention_child2)		
	else:
		retention_score = percentage_retention((child_node1_descendants.union(child_node2_descendants)), species_tree_descendants)
		type = "Speciation" if len(species_tree_descendants) > 1 else "Terminal_Speciation"
	return [type, float(retention_score)]
	
def get_descending_species(node, tree):
	tree_node = tree.search_nodes(name = node)[0]
	species_IDs = []
	for leaf in tree_node.iter_leaves():
		species = ((leaf.name).split('_'))[0] if "_" in leaf.name else leaf.name
		species_IDs.append(species)
	return(set(species_IDs))
	
def percentage_retention(sample_set, total_set):                               # Find the number of species from the total_set (species tree) that are retained in the sample_set (orthogroup tree)
	overlap = sample_set.intersection(total_set)
	retention = (float(len(overlap)) / float(len(total_set)))*100
	return(retention)

def main():	
	os.chdir(species_tree_location)
	species_tree = Tree("Phytozome10_constrainedTree_rooted_labelled.tree", format=1)
	number_of_orthogroups = 0
	all_results = []

	os.chdir(orthogroup_tree_recon_location)
	for filename in glob.glob("OG0*.locus.tree"):
		number_of_orthogroups += 1
		orthogroup = (filename[:-11])                 # Pull out the orthogroup number
		print orthogroup

		orthogroup_tree = Tree(filename, format=1)
		orthogroup_recon_file = (orthogroup + ".locus.recon")
		node2node_dict = link_nodes(orthogroup_recon_file)                           # Make a dictionary linking the nodes in the gene tree to nodes in the species tree
		
		os.chdir(location_of_localisation_predictions)
		protein_localisation_dict = get_protein_locations(orthogroup)
		
		os.chdir(location_of_ACE_likelihoods)
		node_ancestral_state_dict = get_ancestral_likelihoods(orthogroup)            # Fill dict in format - node: likelihood chloroplast, likelihood_mitochondria, likelihood_secretory, likelihood_peroxisome
		
		for node in orthogroup_tree.iter_descendants():
			if "_" in node.name:      # Ignore if it is a species (i.e. leaf or terminal branch) and not a node...may change this later to incorporate 
				continue
			species_tree_node = node2node_dict[node.name]
			relocalisations = get_relocalisations(node, node.up, node_ancestral_state_dict, protein_localisation_dict)
			node_type = get_duplication_score(node, node2node_dict[node.name], species_tree, orthogroup_tree)
			node_result = tuple([orthogroup, node.name, species_tree_node] + node_type + relocalisations)
			all_results.append(node_result)

		os.chdir(orthogroup_tree_recon_location)

	os.chdir(output_location)		
	with open(result_file, "wb") as outfile:
		w = csv.writer(outfile, delimiter = ',', quotechar = ' ', quoting = csv.QUOTE_MINIMAL)
		w.writerow(["orthogroup", "orthogroup_tree_node", "species_tree_node", "orthogroup_tree_node_event", "retention_score", "C_relocalisation", "C_ingroup_retention", "C_outgroup_retention", "M_relocalisation", "M_ingroup_retention", "M_outgroup_retention", "S_relocalisation", "S_ingroup_retention", "S_outgroup_retention", "P_relocalisation", "P_ingroup_retention", "P_outgroup_retention"])    
		for item in all_results:
			w.writerow(item)
		
		
		
	print number_of_orthogroups	
	
if __name__ == '__main__':
	main()
	
##### QUESTION #### For the outgroup node should I look for the localisation of the node up from the one with the change or the next child node??????????????????????????
	
	

	
