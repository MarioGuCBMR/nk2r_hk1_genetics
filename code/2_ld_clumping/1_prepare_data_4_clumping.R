##############
#INTRODUCTION#
##############

#This code gets the significant variants for each trait.

#################
#Loading library#
#################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

formatting_data_4_clumping <- function(ss_gw){
  #This function formats the data for the clumping.
  
  #STEP 2: now make the SNP column which is SNP:A1:A2
  
  rsid <- ifelse(str_detect(ss_gw$variant, "rs") == TRUE, ss_gw$variant, ss_gw$chr_pos)
  SNP <- paste(rsid, ":", ss_gw$other_allele, ":", ss_gw$effect_allele, sep = "")
  
  #Now we select the columns that we want:
  
  ss_gw$SNP <- SNP
  ss_gw$rsid <- rsid
  
  final_df <- ss_gw %>%
    select(SNP, p_value, chr_pos, effect_allele, other_allele, rsid)
  
  colnames(final_df) <- c("SNP", "pval", "chr_pos", "effect_allele", "other_allele", "rsid")
  
  return(final_df)
  
}

##############
#Loading data#
##############

setwd("H:/From_SUND/TACR2_project/answer_reviewer_5")

input_data <- fread("output/2_haplotype_plot/1_ss_of_interest/hba1c_curated.txt")

########################################
#Get the data ready in the right format#
########################################

data_4_clumping <- formatting_data_4_clumping(input_data)

dir.create("output/2_haplotype_plot/2_clumped_data/1_data_4_clumping/")

fwrite(data_4_clumping, "output/2_haplotype_plot/2_clumped_data/1_data_4_clumping/data_4_clumping.txt")

#Also let's create the other outputs:

dir.create("output/2_haplotype_plot/2_clumped_data/clumped_data/")
dir.create("output/2_haplotype_plot/2_clumped_data/clumped_data/clumped_data")
dir.create("output/2_haplotype_plot/2_clumped_data/clumped_data/qsub_files")
