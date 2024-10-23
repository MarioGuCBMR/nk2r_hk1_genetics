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

path_2_input <- "/projects/kilpelainen-AUDIT/people/zlc436/answer_reviewer_5"

setwd(path_2_input)

annot_data <- fread("../../../../kilpelainen-AUDIT/data/tools/cS2G/files/baselineLD.txt")
annot_data$chr_pos <- paste("chr", annot_data$Chr, ":", annot_data$POS, sep="")

#########################################
#Let's get the loci and see what happens#
#########################################

loci_ <- paste(10, "_", "70929740", "_", "71367422", sep="")
loci_ <- as.vector(loci_)

###################################################
#Let's get the SNPs with high PIP for each loci!!!#
###################################################

for(loci in loci_){
  
  path_2_res <- paste("output/2_haplotype_plot/3_fine_mapped_data/carma_res/", loci, "_annotated.RDS", sep="")
  path_2_ss <- paste("output/2_haplotype_plot/3_fine_mapped_data/loci_ss_aligned/", loci, ".txt", sep="")
  
  ss_df <- fread(path_2_ss)
  
  #Important to know which variants are in the annotations:
  
  ss_df <- ss_df[which(ss_df$chr_pos%in%annot_data$chr_pos),]
  
  carma_res <- readRDS(path_2_res)
  
  index_pips <- which(as.numeric(unlist(carma_res[[1]][1])) > 0.1)
  
  pips <- as.numeric(unlist(carma_res[[1]][1]))[index_pips]
  causal_df <- ss_df[index_pips,]
  
  causal_df$pip <- pips 
  
  output_path <- paste("output/2_haplotype_plot/3_fine_mapped_data/carma_res_clean/", loci, "_01.txt", sep="")
  
  fwrite(causal_df, output_path)
  
}
