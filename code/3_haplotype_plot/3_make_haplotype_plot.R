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

###################
#Loading functions#
###################

mgu_ggLD <- function(data, custom_breaks, custom_labels, hkdc1_breaks,  hk1_breaks, nk2r_breaks, tspan15_breaks, point_labels = point_labels, point_breaks = point_breaks){
  data <- tidyr::as_tibble(data)
  colnames(data) <- c(1:ncol(data))
  n <- length(data)
  
  # Tidy data, only taking unique pairs of data
  values <- data %>%
    dplyr::mutate(idx1 = c(1:nrow(data))) %>%
    tidyr::pivot_longer(!.data$idx1, names_to = "idx2", values_to = "LD") %>%
    dplyr::mutate(dplyr::across(idx2, as.double)) %>%
    dplyr::filter(!duplicated(paste(pmax(.data$idx1, .data$idx2), pmin(.data$idx1, .data$idx2), sep = "_"))) %>%
    tidyr::unite("id", .data$idx1:.data$idx2, remove = FALSE) %>%
    dplyr::mutate(diff = abs(idx2 - idx1))
  
  # Calculate coordinates for geom_polygon
  positions <- dplyr::bind_rows(values, values, values, values) %>%
    dplyr::group_by(diff, idx1) %>%
    dplyr::mutate(add_index1 = c(0, 1, 0, 1),
                  add_index2 = c(0, -1, 0, 1),
                  minus1_index = c(1, 1, 0, 1)) %>%
    dplyr::mutate(x = diff * 5 / n + 10 / n * (idx1 - minus1_index) + 5 / n * add_index1,
                  y = 5 - diff * 5 / n + 5 / n * add_index2) %>%
    dplyr::ungroup()
  
  #Let's make a dataframe for positions of clumped variants:
  
  positions_breaks <- positions[which(positions$idx1%in%custom_breaks),]
  positions_breaks <- positions_breaks[which(positions_breaks$idx1 == positions_breaks$idx2),]
  positions_breaks <- positions_breaks[which(duplicated(positions_breaks$id) == FALSE),]

  positions_breaks$label <- custom_labels
  
  #Let's change some of them manually so that they do not overlap.
  #This can be a quite tough trick to pull off...
  
  positions_breaks$x[3] <- positions_breaks$x[3]+0.07
  positions_breaks$x[5] <- positions_breaks$x[5]-0.05
  positions_breaks$x[7] <- positions_breaks$x[6]+0.08
  
  #Now let's do a similar dataframe for the other breaks:
  
  positions_hkdc1 <- positions[which(positions$idx1%in%hkdc1_breaks),]
  positions_hkdc1 <- positions_hkdc1[which(positions_hkdc1$idx1 == positions_hkdc1$idx2),]
  positions_hkdc1 <- positions_hkdc1[which(duplicated(positions_hkdc1$id) == FALSE),]
  
  positions_hk1 <- positions[which(positions$idx1%in%hk1_breaks),]
  positions_hk1 <- positions_hk1[which(positions_hk1$idx1 == positions_hk1$idx2),]
  positions_hk1 <- positions_hk1[which(duplicated(positions_hk1$id) == FALSE),]
  
  positions_nk2r <- positions[which(positions$idx1%in%nk2r_breaks),]
  positions_nk2r <- positions_nk2r[which(positions_nk2r$idx1 == positions_nk2r$idx2),]
  positions_nk2r <- positions_nk2r[which(duplicated(positions_nk2r$id) == FALSE),]
  
  positions_tspan15 <- positions[which(positions$idx1%in%tspan15_breaks),]
  positions_tspan15 <- positions_tspan15[which(positions_tspan15$idx1 == positions_tspan15$idx2),]
  positions_tspan15 <- positions_tspan15[which(duplicated(positions_tspan15$id) == FALSE),]
  
  #Finally, let's add the text for the genes in another set of positions:
  
  positions_hkdc1_name <- positions_hkdc1[3,]
  positions_hk1_name <- positions_hk1[3,]
  positions_nk2r_name <- positions_nk2r[2,]
  positions_tspan15_name <- positions_tspan15[3,]
  
  positions_hkdc1_name$label <- "HKDC1"
  positions_hk1_name$label <- "HK1"
  positions_nk2r_name$label <- "NK2R"
  positions_tspan15_name$label <- "TSPAN15"
  
  #Now, finally for real, for realsies, we are going to get the positions for the location
  
  positions_location <- positions[which(positions$idx1%in%point_breaks),]
  positions_location <- positions_location[which(positions_location$idx1 == positions_location$idx2),]
  positions_location <- positions_location[which(duplicated(positions_location$id) == FALSE),]
  
  positions_location$label <- point_labels
  
  #custom_palette <- colorRampPalette(c("darkblue", "lightblue", "red", "darkred"))(n = 10)
  
  # ggplot2
  p <- positions %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_polygon(ggplot2::aes(fill = .data$LD, group = .data$id)) +
    ggplot2::theme_void() +
    #ggplot2::scale_fill_distiller(type = "seq", palette = "YlOrBr", direction = 1)
    ggplot2::scale_fill_gradient(low = "white", high = "red")
  
  
  p_2 <- p + #ggplot2::geom_segment(data= positions_breaks, aes(x = x, xend= x, y=y+0.13, yend = y+0.16), linewidth=0.05) + 
    geom_text(data = positions_breaks, aes(x = x+0.01, y = y+0.42, label = label, angle = 90), vjust = -0.15, size=0.9) 
  
  #p_2 <- p + ggplot2::geom_segment(data= positions_breaks, aes(x = x, xend= x, y=y, yend = y+0.055), linewidth=0.05) +   geom_label(data = positions_breaks, aes(x = x + 0.01, y = y + 0.15, label = label),
  #                                                                                                                                  angle = 90, vjust = -0.15, size = 0.2, fontface = "bold")
  
  p_3 <- p_2 +  ggplot2::geom_segment(aes(x = positions_hkdc1$x[1], xend= positions_hkdc1$x[length(positions_hkdc1$x)], y=5.7, yend = 5.7), linewidth=0.25) + #HKDC1 track line 
    ggplot2::geom_segment(aes(x = positions_hk1$x[1], xend= positions_hk1$x[length(positions_hk1$x)], y=5.7, yend = 5.7), linewidth=0.25) + #HK1 track line 
    ggplot2::geom_segment(aes(x = positions_nk2r$x[1], xend= positions_nk2r$x[length(positions_nk2r$x)], y=5.7, yend = 5.7), linewidth=0.25) + #HK1 track line 
    ggplot2::geom_segment(aes(x = positions_tspan15$x[1], xend= positions_tspan15$x[length(positions_tspan15$x)], y=5.7, yend = 5.7), linewidth=0.25)  #HK1 track line 
    
  p_4 <- p_3 +     
    
    geom_text(data = positions_hkdc1_name, aes(x = x+0.1, y = 5.80, label = label, angle=35), vjust = -0.15, size=1) + 
    geom_text(data = positions_hk1_name, aes(x = x+0.1, y = 5.80, label = label, angle=35), vjust = -0.15, size=1) + 
    geom_text(data = positions_nk2r_name, aes(x = x+0.12, y = 5.80, label = label, angle=35), vjust = -0.15, size=1) +
    geom_text(data = positions_tspan15_name, aes(x = x+0.12, y = 5.80, label = label, angle=35), vjust = -0.15, size=1)
  
  p_5 <- p_4 + ggplot2::geom_segment(aes(x = x[1]-0.01, xend=x[length(x)], y=5.17, yend = 5.17), linewidth=0.15) +
    #ggplot2::geom_segment(aes(x = x, xend=x, y=5.1207, yend = 5.13), linewidth=0.10) +
    ggplot2::geom_text(data = positions_location, aes(x = x+0.01, y = 5.06, label = label, angle = 90), vjust = -0.15, size=0.6) 

  
  return(p_5)
  
}


mgu_ggLD_text <- function(data, custom_breaks, custom_labels, hkdc1_breaks,  hk1_breaks, nk2r_breaks, tspan15_breaks, point_labels = point_labels, point_breaks = point_breaks){
  data <- tidyr::as_tibble(data)
  colnames(data) <- c(1:ncol(data))
  n <- length(data)
  
  # Tidy data, only taking unique pairs of data
  values <- data %>%
    dplyr::mutate(idx1 = c(1:nrow(data))) %>%
    tidyr::pivot_longer(!.data$idx1, names_to = "idx2", values_to = "LD") %>%
    dplyr::mutate(dplyr::across(idx2, as.double)) %>%
    dplyr::filter(!duplicated(paste(pmax(.data$idx1, .data$idx2), pmin(.data$idx1, .data$idx2), sep = "_"))) %>%
    tidyr::unite("id", .data$idx1:.data$idx2, remove = FALSE) %>%
    dplyr::mutate(diff = abs(idx2 - idx1))
  
  # Calculate coordinates for geom_polygon
  positions <- dplyr::bind_rows(values, values, values, values) %>%
    dplyr::group_by(diff, idx1) %>%
    dplyr::mutate(add_index1 = c(0, 1, 0, 1),
                  add_index2 = c(0, -1, 0, 1),
                  minus1_index = c(1, 1, 0, 1)) %>%
    dplyr::mutate(x = diff * 5 / n + 10 / n * (idx1 - minus1_index) + 5 / n * add_index1,
                  y = 5 - diff * 5 / n + 5 / n * add_index2) %>%
    dplyr::ungroup()
  
  #Let's make a dataframe for positions of clumped variants:
  
  positions_breaks <- positions[which(positions$idx1%in%custom_breaks),]
  positions_breaks <- positions_breaks[which(positions_breaks$idx1 == positions_breaks$idx2),]
  positions_breaks <- positions_breaks[which(duplicated(positions_breaks$id) == FALSE),]
  
  positions_breaks$label <- custom_labels
  
  #Let's change some of them manually so that they do not overlap.
  #This can be a quite tough trick to pull off...
  
  positions_breaks$x[3] <- positions_breaks$x[3]+0.05
  positions_breaks$x[5] <- positions_breaks$x[5]-0.05
  positions_breaks$x[7] <- positions_breaks$x[6]+0.07
  
  #Now let's do a similar dataframe for the other breaks:
  
  positions_hkdc1 <- positions[which(positions$idx1%in%hkdc1_breaks),]
  positions_hkdc1 <- positions_hkdc1[which(positions_hkdc1$idx1 == positions_hkdc1$idx2),]
  positions_hkdc1 <- positions_hkdc1[which(duplicated(positions_hkdc1$id) == FALSE),]
  
  positions_hk1 <- positions[which(positions$idx1%in%hk1_breaks),]
  positions_hk1 <- positions_hk1[which(positions_hk1$idx1 == positions_hk1$idx2),]
  positions_hk1 <- positions_hk1[which(duplicated(positions_hk1$id) == FALSE),]
  
  positions_nk2r <- positions[which(positions$idx1%in%nk2r_breaks),]
  positions_nk2r <- positions_nk2r[which(positions_nk2r$idx1 == positions_nk2r$idx2),]
  positions_nk2r <- positions_nk2r[which(duplicated(positions_nk2r$id) == FALSE),]
  
  positions_tspan15 <- positions[which(positions$idx1%in%tspan15_breaks),]
  positions_tspan15 <- positions_tspan15[which(positions_tspan15$idx1 == positions_tspan15$idx2),]
  positions_tspan15 <- positions_tspan15[which(duplicated(positions_tspan15$id) == FALSE),]
  
  #Finally, let's add the text for the genes in another set of positions:
  
  positions_hkdc1_name <- positions_hkdc1[3,]
  positions_hk1_name <- positions_hk1[3,]
  positions_nk2r_name <- positions_nk2r[2,]
  positions_tspan15_name <- positions_tspan15[3,]
  
  positions_hkdc1_name$label <- "HKDC1"
  positions_hk1_name$label <- "HK1"
  positions_nk2r_name$label <- "NK2R"
  positions_tspan15_name$label <- "TSPAN15"
  
  #Now, finally for real, for realsies, we are going to get the positions for the location
  
  positions_location <- positions[which(positions$idx1%in%point_breaks),]
  positions_location <- positions_location[which(positions_location$idx1 == positions_location$idx2),]
  positions_location <- positions_location[which(duplicated(positions_location$id) == FALSE),]
  
  positions_location$label <- point_labels
  
  #custom_palette <- colorRampPalette(c("darkblue", "lightblue", "red", "darkred"))(n = 10)
  
  # ggplot2
  p <- positions %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$y)) +
    #ggplot2::geom_polygon(ggplot2::aes(fill = .data$LD, group = .data$id)) +
    #ggplot2::geom_polygon(data = data.frame(), aes()) +  # Empty geom_polygon
    ggplot2::annotate("rect", xmin = -1, xmax = 10, ymin = -2, ymax = 7, alpha = 0) +
    ggplot2::theme_void() +
    #ggplot2::scale_fill_distiller(type = "seq", palette = "YlOrBr", direction = 1)
    ggplot2::scale_fill_gradient(low = "white", high = "red") 
    #ggplot2::xlim(-2, 7) +
    #ggplot2::ylim(-2, 7)
  
  
  p_2 <- p + ggplot2::geom_segment(data= positions_breaks, aes(x = x, xend= x, y=y+0.13, yend = y+0.16), linewidth=0.05) + 
    geom_text(data = positions_breaks, aes(x = x+0.01, y = y+0.28, label = label, angle = 90), vjust = -0.15, size=0.6) 
  
  #p_2 <- p + ggplot2::geom_segment(data= positions_breaks, aes(x = x, xend= x, y=y, yend = y+0.055), linewidth=0.05) +   geom_label(data = positions_breaks, aes(x = x + 0.01, y = y + 0.15, label = label),
  #                                                                                                                                  angle = 90, vjust = -0.15, size = 0.2, fontface = "bold")
  
  p_3 <- p_2 +  ggplot2::geom_segment(aes(x = positions_hkdc1$x[1], xend= positions_hkdc1$x[length(positions_hkdc1$x)], y=5.5, yend = 5.5), linewidth=0.25) + #HKDC1 track line 
    ggplot2::geom_segment(aes(x = positions_hk1$x[1], xend= positions_hk1$x[length(positions_hk1$x)], y=5.5, yend = 5.5), linewidth=0.25) + #HK1 track line 
    ggplot2::geom_segment(aes(x = positions_nk2r$x[1], xend= positions_nk2r$x[length(positions_nk2r$x)], y=5.5, yend = 5.5), linewidth=0.25) + #HK1 track line 
    ggplot2::geom_segment(aes(x = positions_tspan15$x[1], xend= positions_tspan15$x[length(positions_tspan15$x)], y=5.5, yend = 5.5), linewidth=0.25)  #HK1 track line 
  
  p_4 <- p_3 +     
    
    geom_text(data = positions_hkdc1_name, aes(x = x+0.1, y = 5.60, label = label, angle=35), vjust = -0.15, size=1) + 
    geom_text(data = positions_hk1_name, aes(x = x+0.1, y = 5.60, label = label, angle=35), vjust = -0.15, size=1) + 
    geom_text(data = positions_nk2r_name, aes(x = x+0.12, y = 5.60, label = label, angle=35), vjust = -0.15, size=1) +
    geom_text(data = positions_tspan15_name, aes(x = x+0.12, y = 5.65, label = label, angle=35), vjust = -0.15, size=1)
  
  p_5 <- p_4 + 
    ggplot2::geom_segment(aes(x = -0.01, xend=10, y=5.12, yend = 5.12), linewidth=0.15) +
    #ggplot2::geom_segment(aes(x = x, xend=x, y=5.1207, yend = 5.13), linewidth=0.10) +
    ggplot2::geom_text(data = positions_location, aes(x = x+0.01, y = 5.05, label = label, angle = 90), vjust = -0.15, size=0.5) 
  
  p_6 <- p_5 + ylim(4.75, 6)
  
  return(p_6)
  
}

function.plot.save <- function(p, filename.prefix="MISSING-FILENAME-PREFIX", filename.suffix="MISSING-FILENAME-SUFFIX") {
  ### USE: UNIVERSIAL
  ### Function exports a ggplot object to pdf file
  ### INPUT: 
  # p: a ggplot object
  ### OUTPUT
  # NONE
  # ======================== Save plot =============================== # 
  #p.filename <- sprintf("tissue_plot_%s.pdf", filename.suffix)
  p.filename <- sprintf("tissue_plot_%s_%s.pdf", filename.prefix, filename.suffix)
  suppressWarnings(ggsave(file=p.filename, width=6, height=6, units="in", dpi=100))
  print(sprintf("Saved plot %s", p.filename))
}

function.plot.save_2 <- function(p, filename.prefix="MISSING-FILENAME-PREFIX", filename.suffix="MISSING-FILENAME-SUFFIX") {
  ### USE: UNIVERSIAL
  ### Function exports a ggplot object to pdf file
  ### INPUT: 
  # p: a ggplot object
  ### OUTPUT
  # NONE
  # ======================== Save plot =============================== # 
  #p.filename <- sprintf("tissue_plot_%s.pdf", filename.suffix)
  p.filename <- sprintf("tissue_plot_%s_%s.pdf", filename.prefix, filename.suffix)
  suppressWarnings(ggsave(file=p.filename, width=1800, height=600, units="px", dpi=25))
  print(sprintf("Saved plot %s", p.filename))
}


##############
#Loading data#
##############

path_2_input <- "H:/From_SUND/TACR2_project/answer_reviewer_5/"

setwd(path_2_input)

hba1c <- fread("output/2_haplotype_plot/1_ss_of_interest/hba1c_curated.txt")

clumped_data <- fread("output/2_haplotype_plot/2_clumped_data/clumped_data/clumped_data/hba1c_clumped.txt") #55!!
fm_data <- fread("output/2_haplotype_plot/3_fine_mapped_data/carma_res_clean/10_70929740_71367422_01.txt") #14!!

fm_data <- fm_data[which(as.numeric(fm_data$p_value) < 5e-08),] #trick to avoid issues with fine-mapping

ss_df <- fread("output/2_haplotype_plot/3_fine_mapped_data/loci_ss_aligned/10_70929740_71367422.txt")
ld_df <- readRDS("output/2_haplotype_plot/3_fine_mapped_data/ld_matrices/10_70929740_71367422.RDS")

#################################################################
#Let's filter the data cuz including all SNPs just does not work#
#################################################################

index_ <- which(as.numeric(ss_df$p_value) < 5e-08)

ss_df <- ss_df[index_,]

#################################
#Let's try to generate a heatmap#
#################################

ld_df <- as.data.frame(ld_df)

ld_df <- ld_df[,index_][index_,]

colnames(ld_df) <- ss_df$variant
rownames(ld_df) <- ss_df$variant

ld_df <- ld_df^2

#############################################
#Let's make a dataframe with the proper data#
#############################################

library(ggLD)

n <- dim(ld_df)[1]

#Let's add the clumped variants!

custom_breaks <-   which(ss_df$chr_pos%in%fm_data$chr_pos)#the indexes of where the clumped data is located
custom_labels <- ss_df$variant[custom_breaks]  # Specify the corresponding labels
custom_labels <- unlist(str_split(custom_labels, "_"))
custom_labels <- custom_labels[which(!(str_detect(custom_labels, "A") | str_detect(custom_labels, "C") | str_detect(custom_labels, "G") | str_detect(custom_labels, "T")))]

#Let's get the hk1_breaks to make the gene tracks:

hkdc1_breaks <- which(ss_df$base_pair_location >= 70980088 & ss_df$base_pair_location <= 71027308)
hk1_breaks <- which(ss_df$base_pair_location >= 71029740 & ss_df$base_pair_location <= 71161640)
nk2r_breaks <- which(ss_df$base_pair_location >= 71163659 & ss_df$base_pair_location <= 71176674)
tspan15_breaks <- which(ss_df$base_pair_location >= 71211221 & ss_df$base_pair_location <= 71267422)

#Let's finally add the breaks of the variants when the finish each region:

bp_df <- ss_df %>%
  select(variant, chr_pos, base_pair_location)

#options(digits=5)

bp_df$point <- round(as.numeric(bp_df$base_pair_location)/(1e06), 5)

bp_df <- bp_df[which(duplicated(bp_df$point) == FALSE),]

ok_vect <- c(70.930, 71.000, 71.010, 71.020, 71.030, 71.040, 71.050, 71.060, 71.070, 71.080, 71.090, 
             71.100, 71.110, 71.120, 71.130, 71.140, 71.150, 71.160, 71.170, 71.180, 71.190, 
             71.200, 71.210, 71.220, 71.230, 71.240, 71.250, 71.260, 71.270, 71.280, 71.290, 
             71.300, 71.350, 71.360, 71.365)

point_breaks_df <- bp_df[which(round(bp_df$point,3)%in%round(ok_vect, 3)),]
point_breaks_df <- point_breaks_df[which(duplicated(round(point_breaks_df$point,2)) == FALSE),]

#This did not work perfectly, but works well for us:

point_breaks_df[length(point_breaks_df$variant)] <- bp_df[length(bp_df$variant)]
point_breaks_df <- rbind(point_breaks_df, bp_df[which(bp_df$variant == "rs7902486"),]) #140 missing 
point_breaks_df <- rbind(point_breaks_df, bp_df[which(bp_df$variant == "rs118012763"),]) #210 missing
point_breaks_df <- rbind(point_breaks_df, bp_df[which(bp_df$variant == "rs1665585"),]) #280 missing
point_breaks_df <- rbind(point_breaks_df, bp_df[which(bp_df$variant == "rs2616125"),]) #350 missing

point_breaks_df <- point_breaks_df[order(point_breaks_df$base_pair_location),] #there is no 71.350 so we will have to put the closest one available

point_breaks <- which(ss_df$variant%in%point_breaks_df$variant)
point_labels <- c("70.93", "71.00", "71.01", "71.02", "71.03", "71.04", "71.05", "71.06", "71.07", "71.08", "71.09", 
                  "71.10", "71.11", "71.12", "71.13", "71.14", "71.15", "71.16", "71.17", "71.18", "71.19", 
                  "71.20", "71.21", "71.22", "71.23", "71.24", "71.25", "71.26", "71.27", "71.28", "71.29", 
                  "71.30", "71.35", "71.36", "71.37")

# Let's add the text!

p <- mgu_ggLD_text(ld_df, custom_breaks, custom_labels, hkdc1_breaks,  hk1_breaks, nk2r_breaks, tspan15_breaks, point_labels = point_labels, point_breaks = point_breaks)
  
#function.plot.save(p, filename.prefix="haplotype_guideline_", filename.suffix="plot")

ggsave(filename = "your_plot.pdf", plot = p, device = "pdf", width = 4, height = 2, dpi = 25)

#Let's try saving the data with as much quality as possible:

p <- mgu_ggLD(ld_df, custom_breaks, custom_labels, hkdc1_breaks,  hk1_breaks, nk2r_breaks, tspan15_breaks, point_labels = point_labels, point_breaks = point_breaks)

ggsave(filename = "hq_ld_plot.pdf", plot = p, device = "pdf", width = 4, height = 2, dpi = 25)

p_full <- p

tiff("H:/From_SUND/TACR2_project/answer_reviewer_5/output/2_haplotype_plot/4_plots/1_haplotype_test_90_degrees.tiff", res = 350, width=2400, height=1250)
p_full
dev.off()

