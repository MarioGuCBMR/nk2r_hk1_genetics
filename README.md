# NK2R-HK1 genetics

A repository for the genetic analysis performed for HK1-NK2R-TSPAN15 locus in the paper: "NK2R control of energy expenditure and feeding to treat metabolic diseases"

In particular, this code performs the analysis described below:

***SNP-level associations in the HK1-NK2R-TSPAN15 region***

Utilizing the latest and largest European GWAS of HbA1c available in T2D Knowledge Portal (query on 10/05/2024) [1,2] we retrieved summary statistics from Jurgens et al [3] and took forward 1944 common variants (minimum allele frequency > 1%) located 100 kilobases (kb) upstream from HK1 transcription starting site (TSS) and 100kb downstream from TSPAN15 TSS (10:70929740-71367422). 

We then assessed the number of lead, independent, genome-wide significant variants by performed LD-clumping with with plink1.9 [4] and utilizing the European 1000 Genomes reference panel Phase 3 Version 5 [5]. We set a r2 threshold of 0.01 and distance threshold of 1000 kilobases (kb).

We performed fine-mapping of HbA1c associations in 10:70929740-71367422 locus with CARMA [6], a software designed to correct for differences in linkage disequilibrium (LD) between summary statistics and LD reference panels. We performed CARMA with the default settings, utilizing European 1000 Genomes reference panel Phase 3 Version 5 with annotations [5] and took forward those variants that presented a posterior inclusion probability (PIP) > 0.1. Finally, we retrieved in dbSNP [7] minimum allele frequencies in five genetic ancestries: African American (AFR), Admixed American (AMR), East Asian (EAS), European (EUR) and South Asian (SAS).

Finally, we utilized Open Target Genetics (OTG) [8] to query the variant-to-gene (V2G) scores utilized for gene prioritization and the associations of causal variants with expression quantitative trait loci (eQTLs). 

## Dependencies:

All analyses were performed in Rstudio (2022.07.2+576) with R (4.1.3). Data were loaded and manipulated using data.table (1.14.2) and tidyverse (1.3.1). Linkage disequilibrium (LD) operations were performed using plink1.9 [4], ggLD [9] and LDLink [10]. 

## Figures:

Raw figures can be found in the output folder, which were then beautified manually. We utilized these raw figures to have an easy represantation of the locus which we used to represent the signals in the final version of the manuscript.

## References:

1) HbA1c phenotype page. 2024 May 21; https://t2d.hugeamp.org/phenotype.html?phenotype=HBA1C (RRID:SCR_003743).
2) Costanzo, M. C. et al.  The Type 2 Diabetes Knowledge Portal: An open access genetic resource dedicated to type 2 diabetes and related traits. Cell metabolism, 35, 695–710 (2023).
3) Jurgens, S. J. et al. Analysis of rare genetic variation underlying cardiometabolic diseases and traits among 200,000 individuals in the UK Biobank. Nature genetics, 54, 240–250 (2022).
4) Purcell S et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. Am J Hum Genet, 81, 559-575 (2007).
5) 1000 Genomes Project Consortium et al. A global reference for human genetic variation. Nature, 526, 68–74. (2015).
6) Yang, Z. et al.CARMA is a new Bayesian model for fine-mapping in genome-wide association meta-analyses. Nature genetics, 55, 1057–1065 (2023).
7) Sherry,S.T. et al. dbSNP—Database for Single Nucleotide Polymorphisms and Other Classes of Minor Genetic Variation. Genome Res., 9, 677–679 (1999).
8) Ghoussaini, M. et al. Open Targets Genetics: systematic identification of trait-associated genes using large-scale genetics and functional genomics. Nucleic acids research, 49, 1311–1320 (2022)
9) https://github.com/mmkim1210/ggLD
10) Machiela MJ et al LDLink: a web-based application for exploring population-specific haplotype structure and linking correlated alleles of possible functional variants. Bioinformatics. 1, 31 (2015).
