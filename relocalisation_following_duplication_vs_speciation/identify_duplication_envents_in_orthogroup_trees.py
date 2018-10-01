#!/usr/bin/env python

# This is a program which identifies duplication nodes in the orthogroup trees produced by Phyldog. 

import os
import glob
from ete3 import Tree
import re
import csv

def link_nodes(recon_file):
	d = {}
	with open(recon_file) as file:
		for line in file:
			token = line.rstrip().split()
			if "n" in token[0]:	              # If the line corresponds to a node and not to a gene id
				d[token[0]] = token[1]
	return d
	
def get_descending_species(node, tree):
	tree_node = tree.search_nodes(name = node)[0]
	species_IDs = []
	for leaf in tree_node.iter_leaves():
		if "_" in leaf.name:
			species = ((leaf.name).split('_'))[0]
		else:
			species = leaf.name
		species_IDs.append(species)
	return(set(species_IDs))
	
def percentage_retention(sample_set, total_set):                               # Find the number of species from the total_set (species tree) that are retained in the sample_set (orthogroup tree)
	overlap = sample_set.intersection(total_set)
	retention = (float(len(overlap)) / float(len(total_set)))*100
	return(retention)
	

def main():	
	os.chdir("/cellar/rona/Phytozome10/Phyldog")
	species_tree = Tree("Phytozome10_constrainedTree_rooted_labelled.tree", format=1)
	number_of_orthogroups = 0
	os.chdir("/home/rona/Phytozome10/Phyldog/output_full_dataset")
	for filename in glob.glob("OG0*.locus.tree"):
		number_of_orthogroups += 1
		orthogroup = (filename[:-11])                 # Pull out the orthogroup number
		print orthogroup
		orthogroup_tree = Tree(filename, format=1)
		orthogroup_recon_file = (orthogroup + ".locus.recon")
		
		# Make a dictionary linking the nodes in the gene tree to nodes in the species tree
		node2node_dict = link_nodes(orthogroup_recon_file)
		
		os.chdir("/home/rona/ancestral_state_estimation/orthogroup_duplications")
		with open((orthogroup + ".duplications.csv"), 'wb') as outfile:
			w = csv.writer(outfile, delimiter = ',', quotechar = ' ', quoting = csv.QUOTE_MINIMAL)
			w.writerow(["orthogroup_node", "species_tree_node", "type", "duplication_score"])         # Duplication score is the lowest of the two child nodes. Type is dup if there is overlap in the two child nodes, terminal duplication if only one descendant species, else it is a speciation 
			for node in orthogroup_tree.iter_descendants():
				if "_" in node.name:      # Ignore if it is a species (i.e. leaf or terminal branch) and not a node...may change this later to incorporate 
					continue
				species_tree_node = node2node_dict[node.name]
				
				##
				species_tree_descendants = get_descending_species(species_tree_node, species_tree)
				child_nodes = re.findall("'([^']*)'", (str((node).children))) 
				child_node1_descendants = get_descending_species(child_nodes[0], orthogroup_tree)
				child_node2_descendants = get_descending_species(child_nodes[1], orthogroup_tree)
				
				if any(i in child_node1_descendants for i in child_node2_descendants) == True:
					retention_child1 = percentage_retention(child_node1_descendants, species_tree_descendants)
					retention_child2 = percentage_retention(child_node2_descendants, species_tree_descendants)
					duplication_score = min(retention_child1, retention_child2)
					if len(species_tree_descendants) > 1:
						type = "Duplication"
					else:
						type = "Terminal_Duplication"
				else:
					type = "Speciation"
					duplication_score = 0.0
							
				w.writerow([node.name, species_tree_node, type, duplication_score])
				
		os.chdir("/home/rona/Phytozome10/Phyldog/output_full_dataset")	
			
	print("\n")		
	print("number of orthogroups: " + str(number_of_orthogroups))

			
if __name__ == '__main__':
	main()	

	
	


	
	
	
