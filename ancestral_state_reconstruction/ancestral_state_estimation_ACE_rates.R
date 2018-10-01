#21/08/2017

# R script which uses the package Phytools and the function ACE to conduct ancestral state estiamtion
# Uses the subcellular localisation binary (0s and 1s) data table for each orthogroup as input
 
# Produces a list of orthogroups from the Phyldog files and then conducts ACE on the genes in each one iteratively and produces a csv file containing the rates of gain and loss 

# 04/09/2017 - updated to read data tables of the format: gene/chloroplast result (1 or 0) / mitochondria result (1 or 0) / secretory result (1 or 0) / Peroxsiome result (1 or 0) / Other location (1 or 0) 

library(phytools)

path = "/home/rona/Phytozome10/Phyldog/output_full_dataset/"
 
# Extract the file names of the orthogroups and create a list of them:
file.names <- dir(path, pattern=".recon")                        
orthogroups <- gsub(".locus.recon", "", file.names)
print(length(orthogroups))

# Initiate an empty data frame which will store the results 
rates <- data.frame("orthogroup" = character(), "size" = integer(), "ratio_chloroplast" = numeric(), "freq_chloroplast" = integer(), "ratio_mitochondria" = numeric(), "freq_mitochondria" = integer(), "ratio_secretory" = numeric(), "freq_secretory" = integer(), "ratio_peroxisome" = numeric(), "freq_peroxisome" = integer(), "ratio_other" = numeric(), "freq_other" = integer(), stringsAsFactors = FALSE) 

for(i in 1:length(orthogroups)){
	setwd('/home/rona/Phytozome10/Phyldog/output_full_dataset/')
	print(orthogroups[i])
	tree <- read.newick(sprintf("%s.locus.tree", orthogroups[i]))
	setwd('/home/rona/subcellular_localisation_prediction/TargetP_PredAlgo_PTS1_PTS2/results_by_orthogroup/RC_5')
	subcellular_location_data <- read.table(sprintf("%s.location.data", orthogroups[i]), header=FALSE)
	rownames(subcellular_location_data) = subcellular_location_data[,1]
	subcellular_location_data <- subcellular_location_data[match(tree$tip.label, subcellular_location_data$V1),]
	if(any(subcellular_location_data$V2 == 1) && any(subcellular_location_data$V2 == 0)){
		ARD_reconstruction_chloroplast <- ace(subcellular_location_data$V2, tree, type = "discrete", model = "ARD")
		rates_chloroplast <- (ARD_reconstruction_chloroplast$rates)
		ratio_chloroplast <- (rates_chloroplast[2])/(rates_chloroplast[1])   #rate of gain over rate of loss
		freq_chloroplast <- (as.data.frame(table(subcellular_location_data$V2)))[2,2]
	} else {
		ratio_chloroplast <- "na"
		freq_chloroplast <- "na"
		}		
	if(any(subcellular_location_data$V3 == 1) && any(subcellular_location_data$V3 == 0)){
		ARD_reconstruction_mitochondria <- ace(subcellular_location_data$V3, tree, type = "discrete", model = "ARD")
		rates_mitochondria <- (ARD_reconstruction_mitochondria$rates)
		ratio_mitochondria <- (rates_mitochondria[2])/(rates_mitochondria[1])   #rate of gain over rate of loss
		freq_mitochondria <- (as.data.frame(table(subcellular_location_data$V3)))[2,2]
	} else {
		ratio_mitochondria <- "na"
		freq_mitochondria <- "na"
		}		
	if(any(subcellular_location_data$V4 == 1) && any(subcellular_location_data$V4 == 0)){
		ARD_reconstruction_secretory <- ace(subcellular_location_data$V4, tree, type = "discrete", model = "ARD")
		rates_secretory <- (ARD_reconstruction_secretory$rates)
		ratio_secretory <- (rates_secretory[2])/(rates_secretory[1])   #rate of gain over rate of loss
		freq_secretory <- (as.data.frame(table(subcellular_location_data$V4)))[2,2]
	} else {
		ratio_secretory <- "na"
		freq_secretory <- "na"
		}		
	if(any(subcellular_location_data$V5 == 1) && any(subcellular_location_data$V5 == 0)){
		ARD_reconstruction_peroxisome <- ace(subcellular_location_data$V5, tree, type = "discrete", model = "ARD")
		rates_peroxisome <- (ARD_reconstruction_peroxisome$rates)
		ratio_peroxisome <- (rates_peroxisome[2])/(rates_peroxisome[1])   #rate of gain over rate of loss
		freq_peroxisome <- (as.data.frame(table(subcellular_location_data$V5)))[2,2]
	} else {
		ratio_peroxisome <- "na"
		freq_peroxisome <- "na"
		}		
	if(any(subcellular_location_data$V6 == 1) && any(subcellular_location_data$V6 == 0)){
		ARD_reconstruction_other <- ace(subcellular_location_data$V6, tree, type = "discrete", model = "ARD")
		rates_other <- (ARD_reconstruction_other$rates)
		ratio_other <- (rates_other[2])/(rates_other[1])   #rate of gain over rate of loss
		freq_other <- (as.data.frame(table(subcellular_location_data$V6)))[2,2]
	} else {
		ratio_other <- "na"
		freq_other <- "na"
		}		
	
	rates[nrow(rates) + 1,] <- c(orthogroups[i], nrow(subcellular_location_data), ratio_chloroplast, freq_chloroplast, ratio_mitochondria, freq_mitochondria, ratio_secretory, freq_secretory, ratio_peroxisome, freq_peroxisome, ratio_other, freq_other)
	
	
}
setwd("/home/rona/ancestral_state_estimation/ancestral_state_rates/RC_5")
write.csv(rates, file="rate_TargetP_PredAlgo_PTS1_raw.csv", row.names=FALSE)














