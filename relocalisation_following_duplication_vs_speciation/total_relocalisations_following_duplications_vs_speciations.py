#!/usr/bin/env python

## Pull out numbers for a hypergeometric test to analyse whether relocalisations are more common following duplication events ## 

import os
import glob
from ete3 import Tree
import csv


location_of_dup_reloc_files = ("/home/rona/ancestral_state_estimation/relocalisations_duplications/TargetP_PredAlgo_PTS1_PTS2_likelihoods/RC_5")
location_of_trees =  ("/home/rona/Phytozome10/Phyldog/output_full_dataset")


class NodeSet(object):
	""" A set is a node and it's two child nodes. Sets have the following properties:
	
	Attributes:	
		Orthogroup: A string that represents the orthogroup the set is from
		Parent_node: The node on the orthogroup tree which is the parent node of the set. This is either a duplication or speciation node
		Child_node_1 and Child_node_2: Either of these nodes might have a relocalisation. If one does then the set is said to be a relocalisation set
		Dict: a dictionary which contains the information duplication and relocalisation info for the set
	"""
	
	def __init__(self, orthogroup, parent_node, child_node_1, child_node_2, dict):
		""" Return a set object """
		self.orthogroup = orthogroup
		self.parent_node = parent_node
		self.child_node_1 = child_node_1
		self.child_node_2 = child_node_2
		self.dict = dict
		
	def get_set_event(self):
		""" Return whether the parent node is a duplication, terminal duplication or speciation node """
		return (self.dict[self.parent_node])[0:2]
		
	def get_set_relocalisations(self, location):
		""" Return whether either child node has a relocalisation, depending on the location or cutoff. 
		A Location can be C, M, S, P or a combination. The location alters the index of the dictionary looked at """
		if location == 'C': s = slice(2,5)
		if location == 'M': s = slice(5,8)
		if location == 'S': s = slice(8,11)
		if location == 'P': s = slice(11,14)
		if location == 'A' : s = slice(2,14)
		
		if "_" in self.child_node_1:
			return self.dict[self.child_node_2][s]
		if "_" in self.child_node_2:
			return self.dict[self.child_node_1][s]
		return self.dict[self.child_node_1][s] + self.dict[self.child_node_2][s]
		
	def is_duplication(self, cutoff):
		""" Is the parent node of the set a duplication with a score that matches the given cutoff? """
		set_event = self.get_set_event()
		return (set_event[0] == 'Duplication' and float(set_event[1]) >= float(cutoff))
		
	def is_relocalisation(self, cutoff):
		relocalisations = self.get_set_relocalisations('A')
		for i in range(len(relocalisations)):
			item = relocalisations[i]
			if item == "gain" or item == "loss":
				if min(float(relocalisations[i + 1]), float(relocalisations[i + 2])) >= cutoff:
					return True
			i += 3
		return False

	def is_speciation(self, cutoff):
		set_event = self.get_set_event()
		return (set_event[0] == 'Speciation' and float(set_event[1]) >= float(cutoff))

def make_node_dict(data_file):
	d = {}
	with open(data_file) as f:
		for line in f:
			token = line.rstrip().split(",")
			d[token[1]] = token[3:]
	return d

def make_node_2_node_dict(orthogroup):
	d = {}
	with open(orthogroup + ".locus.recon") as f:
		for line in f:
			token = line.rstrip().rsplit()
			d[token[0]] = token[1]
	return d
	
number_of_sets = 0
number_of_sets_with_duplication = 0
number_of_sets_with_speciation = 0
number_of_sets_with_relocalisation = 0
number_of_sets_with_duplication_and_relocalisation = 0
number_of_sets_with_speciation_and_relocalisation = 0

# Load Species Tree
os.chdir("/cellar/rona/Phytozome10/Phyldog")
species_tree = Tree("Phytozome10_constrainedTree_rooted_labelled.tree", format=1)
# Create dictionary for each location, with each node on the species tree as a key with a starting value of 0. This value will be updated if there is a loss/gain at that node.
number_of_relocalisations_on_species_tree_node = {}
for node in species_tree.traverse():
	number_of_relocalisations_on_species_tree_node[node.name] = [0, 0, 0, 0]  #format - node.name: relocalisations following dup, total_dups, relocalisations from spec, total_specs


os.chdir(location_of_trees)
for filename in glob.glob("OG0*.locus.tree"):
	orthogroup = (filename[:-11])
	node_2_node_dict = make_node_2_node_dict(orthogroup)    # ortho_node = species_node

	orthogroup_tree = Tree(filename, format=1)
	os.chdir(location_of_dup_reloc_files)
	duplication_relocalisation_file = ("%s.relocalisations_and_duplication.csv" %orthogroup)
	node_dict = make_node_dict(duplication_relocalisation_file)

	for node in orthogroup_tree.iter_descendants():
		if "_" in node.name:      # Ignore if it is a species (i.e. leaf or terminal branch) and not a node...may change this later to incorporate
			continue
		if "_" in (node.children[0]).name and "_" in (node.children[1]).name:
			continue

		node_set = NodeSet(orthogroup, node.name, (node.children[0]).name, (node.children[1]).name, node_dict)
		if node_set.get_set_event()[0] == 'Terminal_Duplication':
			continue
		number_of_sets += 1
		boolean1 = node_set.is_duplication(100.0)
		boolean2 = node_set.is_relocalisation(0.75)
		boolean3 = node_set.is_speciation(100.0)

		if boolean1: number_of_sets_with_duplication += 1
		if boolean2: number_of_sets_with_relocalisation += 1
		if boolean3: number_of_sets_with_speciation += 1
		if boolean1 and boolean2: number_of_sets_with_duplication_and_relocalisation += 1
		if boolean3 and boolean2: number_of_sets_with_speciation_and_relocalisation += 1

		if boolean1:
			number_of_relocalisations_on_species_tree_node[node_2_node_dict[node.name]][1] += 1
			if boolean2:
				number_of_relocalisations_on_species_tree_node[node_2_node_dict[node.name]][0] += 1

		if boolean3:
			number_of_relocalisations_on_species_tree_node[node_2_node_dict[node.name]][3] += 1
			if boolean2:
				number_of_relocalisations_on_species_tree_node[node_2_node_dict[node.name]][2] += 1

	os.chdir(location_of_trees)


print "Number of sets in whole population = " + str(number_of_sets)
print "Number of sets in sample population (with a duplication) = " + str(number_of_sets_with_duplication)
print "Number of successes in whole population (sets with a relocalisation) = " + str(number_of_sets_with_relocalisation)
print "Number of successes in the sample pop (sets with a duplication and relocalisation) = " + str(number_of_sets_with_duplication_and_relocalisation)

print "Number of sets with a speciation = " + str(number_of_sets_with_speciation)
print "Number of sets with a speciation and relocalisation = " + str(number_of_sets_with_speciation_and_relocalisation)


for key in number_of_relocalisations_on_species_tree_node:
	if key.startswith("N"):
		print str(key) + ": " + str(number_of_relocalisations_on_species_tree_node[key])

# os.chdir("/cellar/rona/Scripts/relocalisations_and_duplications")
# with open("tally_of_filtered_results_all_orthogroups_18_04_18.csv", 'wb') as f:
# 	w = csv.writer(f, delimiter = ',', quotechar = ' ', quoting = csv.QUOTE_MINIMAL)
# 	w.writerow(["Species_tree_node", "reloc_following_dup", "number_of_dups", "reloc_following_spec", "number_of_specs"])
# 	for key in number_of_relocalisations_on_species_tree_node:
# 		if key.startswith("N"):
# 			w.writerow([key, (number_of_relocalisations_on_species_tree_node[key])[0], (number_of_relocalisations_on_species_tree_node[key])[1], (number_of_relocalisations_on_species_tree_node[key])[2], (number_of_relocalisations_on_species_tree_node[key])[3]])
# 		else:
# 			continue