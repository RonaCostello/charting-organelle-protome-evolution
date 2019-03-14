#!/usr/bin/python

# A program which identifies changes in ancestral state along an orthogroup tree and applies a filtration method to identify only those changes which have been retained in the majority of descendant species.

import glob
import csv
from ete3 import Tree
import os
import re
import time
import sys

if len(sys.argv) < 3:
	print("Please specify a relocalisation retention score between 0.0 and 1.0 and name of the output file")
	print("e.g.")
	print("script 0.75 output_0.75")
	quit()

result_file = (sys.argv[2] + ".csv")
species_tree_tally_results = (result_file + ".tally.csv")
species_tree_location = ("/home/rona/Phytozome10/Phyldog")
location_of_orthogroup_trees = ("/home/rona/Phytozome10/Phyldog/output_full_dataset")
location_of_localisation_predictions = ("/home/rona/chapter-1/subcellular_localisation_prediction/TargetP_PredAlgo_PTS1_PTS2/results_by_orthogroup/RC_5")
location_of_ACE_likelihoods =("/home/rona/chapter-1/ancestral_state_estimation/ancestral_state_likelihoods/TargetP_PredAlgo_PTS1_PTS2_likelihoods/RC_5")
output_location = ("/home/rona/chapter-1/ancestral_state_estimation/git-repos/charting-organelle-protome-evolution/ancestral_state_reconstruction")

relocalisation_retention_score = float(sys.argv[1])

# Function that takes an orthogroup gene tree and uses the associated recon file to build a dictionary with each node of the gene tree as a key with its corresponding node on the species tree as a value
def map_nodes_to_species_tree(orthogroup):
	gene_tree_species_tree_node_match = {}
	with open("%s.locus.recon" %orthogroup) as file:
		for line in file:
			token = line.rstrip().split()
			if "n" in token[0]:	              # If the line corresponds to a node and not to a gene id
				gene_tree_species_tree_node_match[token[0]] = token[1]
	return gene_tree_species_tree_node_match


def get_gene_locations(orthogroup):
	subcellular_location_of_gene = {}
	with open("%s.location.data" %(orthogroup)) as file:
		for line in file:
			token = line.rstrip().split()
			subcellular_location_of_gene[token[0]] = (token[1], token[2], token[3], token[4])  # Dict: geneID = chloro(1/0), mito(1/0), signal(1/0), perox(1/0)
	return subcellular_location_of_gene

def get_ancestral_states(orthogroup):
	ancestral_state_at_orthogroup_node = {}
	with open("%s.TargetP_PredAlgo_peroxisomal_likelihoods.csv" %(orthogroup)) as csvfile:
		reader = csv.DictReader(csvfile)
		for row in reader:
			# Fill dict in format - node: likelihood chloroplast, likelihood_mitochondria, likelihood_secretory, likelihood_PTS1
			ancestral_state_at_orthogroup_node[row['node']] = [row['likelihoods_chloroplast'], row['likelihoods_mitochondria'], row['likelihoods_secretory'], row['likelihoods_peroxisomal']]
	for key in ancestral_state_at_orthogroup_node:
		for i in range(4):
			if ancestral_state_at_orthogroup_node[key][i] == "na" or ancestral_state_at_orthogroup_node[key][i] == "NA":
				ancestral_state_at_orthogroup_node[key][i] = 'na'
				continue
			if "-" in str(ancestral_state_at_orthogroup_node[key][i]):
				ancestral_state_at_orthogroup_node[key][i] = 0
				continue
			if float(ancestral_state_at_orthogroup_node[key][i]) >= 0.5:
				ancestral_state_at_orthogroup_node[key][i] = 1
				continue
			else:
				ancestral_state_at_orthogroup_node[key][i] = 0
				continue
	return ancestral_state_at_orthogroup_node

def identify_changes(orthogroup, orthogroup_tree, number_of_genes, ancestral_state_at_orthogroup_node, subcellular_location_of_gene, gain_and_loss_on_species_tree, gene_tree_species_tree_node_match, outwriter):
	for i in range(4):                #### MUST BE 4 IF YOU WANT TO CONSIDER CHLORO, MITO, SIGNALLING, AND PEROX!! ####
		if i == 0:
			L = "C"
		if i == 1:
			L = "M"
		if i == 2:
			L = "S"
		if i == 3:
			L = "P"
		for node in orthogroup_tree.iter_descendants():
			if "_" in node.name:      # Ignore if it is a species (i.e. leaf or terminal branch) and not a node...may change this later to incorporate
				continue
			if ancestral_state_at_orthogroup_node[node.name][i] > ancestral_state_at_orthogroup_node[(node.up).name][i]:  # if there is a gain
				child_nodes = re.findall("'([^']*)'", (str((node.up).children)))
				if child_nodes[0] == node.name:                          # Find the node that is the ancestor of the outgroup, i.e the unchanged lineage.
					outgroup_node_name = child_nodes[1]
				else:
					outgroup_node_name = child_nodes[0]
				outgroup_node = orthogroup_tree.search_nodes(name = outgroup_node_name)[0]
				# Start tallys for each of the groups
				total_leaves = 0.0
				total_leaves_outgroup = 0.0
				no_retained_in_descendants = 0.0   # that is the leaves under the changed nodes
				no_retained_in_outgroup = 0.0
				for leaf in node.iter_leaves():
					total_leaves += 1.0
					gene = (str(leaf))[3:]
					if int(subcellular_location_of_gene[gene][i]) == 1:
						no_retained_in_descendants += 1.0
				for leaf in outgroup_node.iter_leaves():
					total_leaves_outgroup += 1.0
					gene = (str(leaf))[3:]
					if int(subcellular_location_of_gene[gene][i]) == 0:
						no_retained_in_outgroup += 1.0
				# Now check if for this node, the change passes the 2-way filtration. Can adjust the cutoff for % retention.
				# Also can make it a one-way filtration by blocking out the second condition.
				if (no_retained_in_descendants/total_leaves) > relocalisation_retention_score and (no_retained_in_outgroup/total_leaves_outgroup) > relocalisation_retention_score:
					gain_and_loss_on_species_tree[gene_tree_species_tree_node_match[node.name]][2*(i)] += 1
					outwriter.writerow([orthogroup, L, "G", node.name, (node.up).name, gene_tree_species_tree_node_match[node.name], number_of_genes])

			# Repeat if there is a loss:
			if ancestral_state_at_orthogroup_node[node.name][i] < ancestral_state_at_orthogroup_node[(node.up).name][i]:  # if there is a loss
				child_nodes = re.findall("'([^']*)'", (str((node.up).children)))
				if child_nodes[0] == node.name:
					outgroup_node_name = child_nodes[1]
				else:
					outgroup_node_name = child_nodes[0]
				outgroup_node = orthogroup_tree.search_nodes(name = outgroup_node_name)[0]
				total_leaves = 0.0
				total_leaves_outgroup = 0.0
				no_retained_in_descendants = 0.0
				no_retained_in_outgroup = 0.0
				for leaf in node.iter_leaves():
					total_leaves += 1.0
					gene = (str(leaf))[3:]
					if int(subcellular_location_of_gene[gene][i]) == 0:
						no_retained_in_descendants += 1.0
				for leaf in outgroup_node.iter_leaves():
					total_leaves_outgroup += 1.0
					gene = (str(leaf))[3:]
					if int(subcellular_location_of_gene[gene][i]) == 1:
						no_retained_in_outgroup += 1.0
				if (no_retained_in_descendants/total_leaves) > relocalisation_retention_score and (no_retained_in_outgroup/total_leaves_outgroup) > relocalisation_retention_score:
					gain_and_loss_on_species_tree[gene_tree_species_tree_node_match[node.name]][(2*(i)) + 1] += 1
					outwriter.writerow([orthogroup, L, "L", node.name, (node.up).name, gene_tree_species_tree_node_match[node.name], number_of_genes])

def main():
	os.chdir(output_location)
	with open("2-way_filtered_results." + timestr + ".csv", 'wb') as csvfile:
		outwriter = csv.writer(csvfile, delimiter = ',', quotechar = ' ', quoting = csv.QUOTE_MINIMAL)
		outwriter.writerow(["orthogroup", "Location", "Change", "node", "node_up", "species_tree_node", "number of genes"])
		os.chdir(species_tree_location)
		species_tree = Tree("Phytozome10_constrainedTree_rooted_labelled.tree", format=1)
		# Initiate a dictionary which will hold the numbers of gains and losses at each node in the species tree (all start with a value of zero).
		gain_and_loss_on_species_tree = {} # node: Cgain, Closs, Mgain, Mloss, Sgain, Sloss, Pgain, Ploss
		for node in species_tree.traverse():
			gain_and_loss_on_species_tree[node.name] = [0, 0, 0, 0, 0, 0, 0, 0]  # Format chloro gain, chloro loss, mito gain, mito loss, secretory gain, secretory loss, PTS1 gain, PTS1 loss

		number_of_orthogroups = 0
		os.chdir(location_of_orthogroup_trees)
		for filename in glob.glob("OG*.locus.tree"):      # Iterate through the orthogroup (gene) tree files
			number_of_orthogroups += 1
			orthogroup = (filename[:-11])                 # Pull out the orthogroup number
			orthogroup_tree = Tree(filename, format=1)
			number_of_genes = 0
			for leaf in orthogroup_tree:                  # The number of leaves in the tree = the number of genes in the orthogroup
				number_of_genes += 1
			if int(number_of_genes) < 50:                 # Skip orthogroups of less than 50 genes
				continue
			print(orthogroup)

			gene_tree_species_tree_node_match = map_nodes_to_species_tree(orthogroup)     # Fill dict [gene tree node] = gene tree node, duplication or speciation
			os.chdir(location_of_localisation_predictions)
			subcellular_location_of_gene = get_gene_locations(orthogroup)            # For each gene in the orthogroup, initiate a dict to store its predicited subcellular localisation
			os.chdir(location_of_ACE_likelihoods)
			ancestral_state_at_orthogroup_node = get_ancestral_states(orthogroup)        # For each node in the orthogroup tree, initiate a dict to store its ancestral state estimation

			## Identify changes in ancestral state on the orthogroup tree
			identify_changes(orthogroup, orthogroup_tree, number_of_genes, ancestral_state_at_orthogroup_node, subcellular_location_of_gene, gain_and_loss_on_species_tree, gene_tree_species_tree_node_match, outwriter)
			os.chdir(location_of_orthogroup_trees)

	for key in gain_and_loss_on_species_tree:
		if key.startswith("N"):
			print(key + " " + str(gain_and_loss_on_species_tree[key]))
		else:
			continue

	os.chdir(output_location)
	with open("tally_of_filtered_results." + timestr + ".csv", 'wb') as f:
		w = csv.writer(f, delimiter = ',', quotechar = ' ', quoting = csv.QUOTE_MINIMAL)
		w.writerow(["node", "chloroplast gain", "chloroplast lost", "mitochondria gain", "mitochondria lost", "secretory gain", "secretory lost", "peroxisome gain", "peroxisome lost"])
		for key in gain_and_loss_on_species_tree:
			if key.startswith("N"):
				w.writerow([key, (gain_and_loss_on_species_tree[key])[0], (gain_and_loss_on_species_tree[key])[1], (gain_and_loss_on_species_tree[key])[2], (gain_and_loss_on_species_tree[key])[3], (gain_and_loss_on_species_tree[key])[4], (gain_and_loss_on_species_tree[key])[5], (gain_and_loss_on_species_tree[key])[6],(gain_and_loss_on_species_tree[key])[7]])
			else:
				continue

	print("Number of orthogroups: " + str(number_of_orthogroups))

	# Script to turn the csv output of the tallied changes into an array for illustrator
	# NOT TESTED YET #
	with open("tally_of_filtered_results." + timestr + ".csv") as file:
		content = file.readlines()[1:]
	content = [x.strip() for x in content]
	with open("tally_of_filtered_results.tally." + timestr + ".txt", "w") as f:
		f.write(str(content))


if __name__ == '__main__':
	main()
