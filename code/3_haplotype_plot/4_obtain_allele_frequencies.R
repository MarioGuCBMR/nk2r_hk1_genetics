##############
#INTRODUCTION#
##############

#Let's make the haplotype plot.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(ggrepel)
library(otargen)

##############
#Loading data#
##############

path_2_input <- "H:/From_SUND/TACR2_project/answer_reviewer_5/"

setwd(path_2_input)

hba1c <- fread("output/2_haplotype_plot/1_ss_of_interest/hba1c_curated.txt")

clumped_data <- fread("output/2_haplotype_plot/2_clumped_data/clumped_data/clumped_data/strict/hba1c_strict_clumped.txt") #15!!
fm_data <- fread("output/2_haplotype_plot/3_fine_mapped_data/carma_res_clean/10_70929740_71367422_01.txt") #14!!

##################################
#STEP 1: let's assess the overlap#
##################################

#1 get only genome-wide significant variants:

fm_gw <- fm_data[which(as.numeric(fm_data$p_value) < 5e-08),] #9

length(which(fm_gw$chr_pos%in%clumped_data$chr_pos)) #just one overlap - goes to show that ld-clumping is a bit shitty

##########################################################################
#Let's do two tables, then, one with independent and one with fine-mapped#
##########################################################################

clumped_hba1c <- hba1c[which(hba1c$chr_pos%in%clumped_data$chr_pos),]

fm_hba1c <- hba1c[which(hba1c$chr_pos%in%fm_gw$chr_pos),]

#Let's align to the positive allele:

clumped_hba1c_pos <- clumped_hba1c

new_a1 <- ifelse(as.numeric(clumped_hba1c$beta) < 0, clumped_hba1c$other_allele, clumped_hba1c$effect_allele)
new_a2 <- ifelse(as.numeric(clumped_hba1c$beta) < 0, clumped_hba1c$effect_allele, clumped_hba1c$other_allele)
new_beta <- ifelse(as.numeric(clumped_hba1c$beta) < 0, as.numeric(clumped_hba1c$beta)*(-1), as.numeric(clumped_hba1c$beta))

clumped_hba1c_pos$effect_allele <- new_a1
clumped_hba1c_pos$other_allele <- new_a2
clumped_hba1c_pos$beta <- new_beta

#Let's do the same for fine-mapped variants:

fm_hba1c_pos <- fm_hba1c

new_a1 <- ifelse(as.numeric(fm_hba1c$beta) < 0, fm_hba1c$other_allele, fm_hba1c$effect_allele)
new_a2 <- ifelse(as.numeric(fm_hba1c$beta) < 0, fm_hba1c$effect_allele, fm_hba1c$other_allele)
new_beta <- ifelse(as.numeric(fm_hba1c$beta) < 0, as.numeric(fm_hba1c$beta)*(-1), as.numeric(fm_hba1c$beta))

fm_hba1c_pos$effect_allele <- new_a1
fm_hba1c_pos$other_allele <- new_a2
fm_hba1c_pos$beta <- new_beta

##################################
#Let's add the allele frequencies#
##################################

clumped_hba1c_pos$afr_maf <- NA
clumped_hba1c_pos$amr_maf <- NA
clumped_hba1c_pos$eas_maf <- NA
clumped_hba1c_pos$eur_maf <- NA
clumped_hba1c_pos$consequence <- NA
clumped_hba1c_pos$gene <- NA

skip_to_next <- FALSE

for(index_variant in seq(1, length(clumped_hba1c_pos$variant))){
  
  # Let's make sure we retrieve the build37 data...
  
  variant_info <- tryCatch(variantInfo(clumped_hba1c_pos$variant[index_variant]), error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next){
    
    skip_to_next <- FALSE
    next()
    
  }  
  
  clumped_hba1c_pos$afr_maf[index_variant] <- variant_info$gnomadAFR
  clumped_hba1c_pos$amr_maf[index_variant] <- variant_info$gnomadAMR
  clumped_hba1c_pos$eas_maf[index_variant] <- variant_info$gnomadEAS
  clumped_hba1c_pos$eur_maf[index_variant] <- variant_info$gnomadNFE
  clumped_hba1c_pos$gene[index_variant] <- variant_info$nearestCodingGene.symbo
  clumped_hba1c_pos$consequence[index_variant] <- variant_info$mostSevereConsequence
  
}

#We will get the gene data later - let's clean the data now:

clean_clumped_hba1c_pos <- clumped_hba1c_pos %>%
  select(variant, chr_pos, effect_allele, other_allele, afr_maf, amr_maf, eas_maf, eur_maf, beta, standard_error, p_value)

write.csv(clean_clumped_hba1c_pos, "tables/clumped_variants_w_all_ancestries_strict.csv", quote=FALSE)

################################################
#Let's do the same for the fine-mapped variants#
################################################

fm_hba1c_pos$afr_maf <- NA
fm_hba1c_pos$amr_maf <- NA
fm_hba1c_pos$eas_maf <- NA
fm_hba1c_pos$eur_maf <- NA
fm_hba1c_pos$consequence <- NA
fm_hba1c_pos$gene <- NA

skip_to_next <- FALSE

for(index_variant in seq(1, length(fm_hba1c_pos$variant))){
  
  # Let's make sure we retrieve the build37 data...
  
  variant_info <- tryCatch(variantInfo(fm_hba1c_pos$variant[index_variant]), error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next){
    
    skip_to_next <- FALSE
    next()
    
  }  
  
  fm_hba1c_pos$afr_maf[index_variant] <- variant_info$gnomadAFR
  fm_hba1c_pos$amr_maf[index_variant] <- variant_info$gnomadAMR
  fm_hba1c_pos$eas_maf[index_variant] <- variant_info$gnomadEAS
  fm_hba1c_pos$eur_maf[index_variant] <- variant_info$gnomadNFE
  fm_hba1c_pos$gene[index_variant] <- variant_info$gnomadNFE
  fm_hba1c_pos$consequence[index_variant] <- variant_info$mostSevereConsequence
  
}

#We will get the gene data later - let's clean the data now:

clean_fm_hba1c_pos <- fm_hba1c_pos %>%
  select(variant, chr_pos, effect_allele, other_allele, afr_maf, amr_maf, eas_maf, eur_maf, beta, standard_error, p_value)

write.csv(clean_fm_hba1c_pos, "tables/causal_variants_w_all_ancestries.csv")
