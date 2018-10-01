#21/08/2017

# R script which uses the package Phytools and the function ACE to conduct ancestral state estiamtion
# Uses the subcellular localisation binary (0s and 1s) data table for each orthogroup as input
 
# Creates a list of orthogroups from the Phyldog directory and then conducts ACE on each one iteratively saving the likelihoods at each node to a data frame and then producing a file for each orthogroup

# Update 04/09/2017: Conduct ACE for chloroplast, mitochondria, secretory and peroxisomal states

library(phytools)

#location of the phyldog trees and recon files:
directory_of_phyldog_output = '/home/rona/Phytozome10/Phyldog/output_full_dataset'


file.names <- dir(directory_of_phyldog_output, pattern=".recon")
# Create a list of all the orthogroups, which can be later iterated over:
orthogroups <- gsub(".locus.recon", "", file.names)


for(i in 1:length(orthogroups)){
	likelihoods <- data.frame("node" = character(), "likelihoods_not_chloroplast" = numeric(), "likelihoods_chloroplast" = numeric(), "likelihoods_not_mitochondria" = numeric(), "likelihoods_mitochondria" = numeric(), "likelihoods_not_secretory" = numeric(), "likelihoods_secretory" = numeric(), "likelihoods_not_peroxisomal" = numeric(), "likelihoods_peroxisomal" = numeric(), stringsAsFactors = FALSE) 
	node_list <- 0
	setwd(directory_of_phyldog_output)
	print(orthogroups[i])
	tree <- read.newick(sprintf("%s.locus.tree", orthogroups[i]))	
	# print(tree$node.label)
	# print(length(tree$node.label))
	setwd('/home/rona/subcellular_localisation_prediction/TargetP_PredAlgo_PTS1_PTS2/results_by_orthogroup/RC_5')
	subcellular_location_data <- read.table(sprintf("%s.location.data", orthogroups[i]), header=FALSE)
	rownames(subcellular_location_data) = subcellular_location_data[,1]
	subcellular_location_data <- subcellular_location_data[match(tree$tip.label, subcellular_location_data$V1),]
	#print(length(row.names(subcellular_location_data))) #prints the number of genes in the orthogroup
	
	if(any(subcellular_location_data$V2 == 1) && any(subcellular_location_data$V2 == 0)){
		ARD_reconstruction_chloroplast <- ace(subcellular_location_data$V2, tree, type = "discrete", model = "ARD")
		likelihoods_of_chloroplast_recon <- as.data.frame(ARD_reconstruction_chloroplast$lik.anc)
		row.names(likelihoods_of_chloroplast_recon) <- tree$node.label
		#print(length(rownames(likelihoods_of_chloroplast_recon)))
		likelihoods_not_chloroplast <- likelihoods_of_chloroplast_recon[,1]
		likelihoods_chloroplast <- likelihoods_of_chloroplast_recon[,2]
		node_list <- (rownames(likelihoods_of_chloroplast_recon))
	} else {
		likelihoods_not_chloroplast <- as.list(rep("na", length(tree$node.label)))
		likelihoods_chloroplast <- as.list(rep("na", length(tree$node.label)))

	}
	if(any(subcellular_location_data$V3 == 1) && any(subcellular_location_data$V3 == 0)){
		ARD_reconstruction_mitochondria <- ace(subcellular_location_data$V3, tree, type = "discrete", model = "ARD")
		likelihoods_of_mitochondria_recon <- as.data.frame(ARD_reconstruction_mitochondria$lik.anc)
		row.names(likelihoods_of_mitochondria_recon) <- tree$node.label
		likelihoods_not_mitochondria <- likelihoods_of_mitochondria_recon[,1]
		likelihoods_mitochondria <- likelihoods_of_mitochondria_recon[,2]
		node_list <- (rownames(likelihoods_of_mitochondria_recon))
	} else {
		likelihoods_not_mitochondria <- as.list(rep("na", length(tree$node.label)))
		likelihoods_mitochondria <- as.list(rep("na", length(tree$node.label)))
	}
	if(any(subcellular_location_data$V4 == 1) && any(subcellular_location_data$V4 == 0)){
		ARD_reconstruction_secretory <- ace(subcellular_location_data$V4, tree, type = "discrete", model = "ARD")
		likelihoods_of_secretory_recon <- as.data.frame(ARD_reconstruction_secretory$lik.anc)
		row.names(likelihoods_of_secretory_recon) <- tree$node.label
		likelihoods_not_secretory <- likelihoods_of_secretory_recon[,1]
		likelihoods_secretory <- likelihoods_of_secretory_recon[,2]
		node_list <- (rownames(likelihoods_of_secretory_recon))
	} else {
		likelihoods_not_secretory <- as.list(rep("na", length(tree$node.label)))
		likelihoods_secretory <- as.list(rep("na", length(tree$node.label)))
	}
	if(any(subcellular_location_data$V5 == 1) && any(subcellular_location_data$V5 == 0)){
		ARD_reconstruction_peroxisomal <- ace(subcellular_location_data$V5, tree, type = "discrete", model = "ARD")
		likelihoods_of_peroxisomal_recon <- as.data.frame(ARD_reconstruction_peroxisomal$lik.anc)
		row.names(likelihoods_of_peroxisomal_recon) <- tree$node.label
		likelihoods_not_peroxisomal <- likelihoods_of_peroxisomal_recon[,1]
		likelihoods_peroxisomal <- likelihoods_of_peroxisomal_recon[,2]
		node_list <- (rownames(likelihoods_of_peroxisomal_recon))
	} else {
		likelihoods_not_peroxisomal <- as.list(rep("na", length(tree$node.label)))
		likelihoods_peroxisomal <- as.list(rep("na", length(tree$node.label)))
	}
	# print(length(tree$node.label))
	if(node_list == 0){
		for(x in 1:(length(tree$node.label))){
			likelihoods[nrow(likelihoods) + 1,] <- c( (tree$node.label)[x], likelihoods_not_chloroplast[x], likelihoods_chloroplast[x], likelihoods_not_mitochondria[x], likelihoods_mitochondria[x], likelihoods_not_secretory[x], likelihoods_secretory[x], likelihoods_not_peroxisomal[x], likelihoods_peroxisomal[x])
	}} else {
		for(x in 1:(length(tree$node.label))){
		likelihoods[nrow(likelihoods) + 1,] <- c( (node_list)[x], likelihoods_not_chloroplast[x], likelihoods_chloroplast[x], likelihoods_not_mitochondria[x], likelihoods_mitochondria[x], likelihoods_not_secretory[x], likelihoods_secretory[x], likelihoods_not_peroxisomal[x], likelihoods_peroxisomal[x])
	}}	
	
	setwd("/home/rona/ancestral_state_estimation/ancestral_state_likelihoods/TargetP_PredAlgo_PTS1_PTS2_likelihoods/RC_5")
	write.csv(likelihoods, file=(sprintf("%s.TargetP_PredAlgo_peroxisomal_likelihoods.csv", orthogroups[i])), row.names=FALSE)
}

print("Number of orthogroups " + length(orthogroups))

###### NA may be returned due to underflow and should only affect v small orthogroups - I checked this and only occured in 3 v. small (<15 genes) orthogroups #######
	
