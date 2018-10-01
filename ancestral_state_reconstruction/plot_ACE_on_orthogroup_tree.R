library(phytools)
#library(ape)
#library(geiger)

## Things hashed out are needed to compute recon trees directly from orthoFinder :) ## 

tree_file <- "OG0002885.locus.tree"
# tree_with_names_file <- "OG0002508.tree_for_ACE_resolved_names.txt"
my_data_file <- "OG0002885.location.data"          
output_file_name <- "OG0002885_ACE_chloro.svg"
colour = 'green'

tree <- (read.newick(tree_file))
#tree_names <- (read.newick(tree_with_names_file))
my_data <- read.table(my_data_file, header = FALSE)
rownames(my_data) = my_data[,1]
my_data <- my_data[match(tree$tip.label, my_data$V1),]

#t2 <- tree
#t2$edge.length[t2$edge.length==0] <- max(nodeHeights(tree))*0.000001
# DO ACE on t2

ARD_reconstruction <- ace(my_data$V2, tree, type = "discrete", model = "ARD")   # 2 = chloroplast, 3 = mitochondria, 4 = secretory, 5 = peroxisome

svg(filename = output_file_name, width = 10, height = 10, pointsize = 9, bg = "white")
plotTree(tree, fsize=0.8, ftype="i") #or tree_names
nodelabels(node=1:tree$Nnode+Ntip(tree), pie=ARD_reconstruction$lik.anc, piecol = c("white", colour), cex=0.3)
tiplabels(pch = 22, bg = ifelse(my_data[tree$tip.label, ][,2],colour,"white"), col="black",adj = c(0.51, 0.5), cex = 0.6) # 2 = chloroplast, 3 = mitochondria, 4 = secretory, 5 = peroxisome

dev.off()