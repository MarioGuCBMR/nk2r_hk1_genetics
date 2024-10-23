##############
#INTRODUCTION#
##############

#This is a code to curate hba1c data from Che et al 2021 (European data)
#Data was downloaded 03/08/2022

#Let's go:

###########
#libraries#
###########

library(data.table)
library(tidyverse)

#######################
#Loading hba1c data#
#######################

#We are gonna load the hba1c from 2024. 

project_path <- "H:/From_SUND/TACR2_project/answer_reviewer_5" #change it with your own path.

setwd(project_path)

hba1c <- fread("raw_data/GWAS_sumstats_EUR__invnorm_glycatedhaemoglobinhba1c__TOTALsample.tsv.gz")

######################
#STEP 1: GENERATE IDs#
######################

#We are gonna set the chr_pos:

hba1c$chr_pos <- paste("chr", hba1c$CHR, ":", hba1c$BP, sep = "")

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

summary(as.numeric(hba1c$CHR)) #no problems

###########################
#STEP 3: REMOVE MAF < 0.01#
###########################

summary(as.numeric(hba1c$A1FREQ)) #We need to remove the NAs, and those with EAF > 0.99 or EAF < 0.01

hba1c_eaf_OK <- hba1c[which(hba1c$A1FREQ > 0.01),]
hba1c_eaf_OK <- hba1c_eaf_OK[which(hba1c_eaf_OK$A1FREQ < 0.99),]

summary(hba1c_eaf_OK$A1FREQ) #worked like a charm.

hba1c_corrected_eaf <- hba1c_eaf_OK
  
##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

hba1c_mhc <- hba1c_corrected_eaf[which(as.numeric(hba1c_corrected_eaf$CHR) == 6 & as.numeric(hba1c_corrected_eaf$BP) >= 26000000 & as.numeric(hba1c_corrected_eaf$BP) <= 34000000),]

summary(as.numeric(hba1c_mhc$CHR)) #perfect!!
summary(as.numeric(hba1c_mhc$BP)) #perfect!!

hba1c_end <- hba1c_corrected_eaf[which(!(hba1c_corrected_eaf$chr_pos%in%hba1c_mhc$chr_pos)),]

#################
#Change columns!#
#################

hba1c_end <- hba1c_end %>%
  select(SNP, CHR, BP, ALLELE1, ALLELE0, A1FREQ, BETA, SE, P, chr_pos)

colnames(hba1c_end) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "chr_pos")

hba1c_end$sample_size <- 460000

#########################
#We can save this data!!#
#########################

fwrite(hba1c_end, "output/1_curated_data/hba1c_curated.txt")

##########################
#Let's save the region!!!#
##########################

hba1c_region <- hba1c_end[which(hba1c_end$chromosome == 10 & hba1c_end$base_pair_location >= 70929740 & hba1c_end$base_pair_location <= 71367422),] #1944

#How many are genome-wide significant?

length(which(hba1c_region$p_value < 5e-08)) #1012

fwrite(hba1c_region, "output/2_haplotype_plot/1_ss_of_interest/hba1c_curated.txt")

