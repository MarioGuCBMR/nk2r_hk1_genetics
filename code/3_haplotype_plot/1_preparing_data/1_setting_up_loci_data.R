##############
#INTRODUCTION#
##############

#This code sets up the input for fine-mapping of several loci of interest.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

path_2_files <-  "/projects/kilpelainen-AUDIT/people/zlc436/answer_reviewer_5" #change it with your own path.

setwd(path_2_files)

#########################################
#Let's get the loci and see what happens#
#########################################

variant_info <- c("10", "70929740", "71367422", "chr10:71112240")
variant_info <- as.data.frame(t(variant_info))
colnames(variant_info) <- c("chromosome", "start_", "end_", "chr_pos")

################################################################
#Alright, let's go and set up the data needed for running CARMA#
################################################################

fwrite(variant_info, "/projects/kilpelainen-AUDIT/people/zlc436/answer_reviewer_5/output/2_haplotype_plot/3_fine_mapped_data/data_loci.txt", quote=FALSE, sep = " ", col.names = FALSE, row.names = FALSE)
