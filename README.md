# NK2R-HK1 genetics

A repository for the genetic analysis performed for NK2R and HK1.

ADD FIGURE here.

## Source of data:

To reproduce the data of these analysis the following data needs to be saved in the /raw_data folder from this repository once it is downloaded.

MAGIC HbA1c trans-ancestry genome-wide associations [Ji Che et al]: https://magicinvestigators.org/downloads/MAGIC1000G_HbA1c_TA.tsv.gz
MAGIC HbA1c trans-ancestry exome-wide associations [Sarah Willems et al]: https://magicinvestigators.org/downloads/MAGIC1000G_HbA1c_TA.tsv.gz

UKBB HbA1c exome associations: Data was downloaded by querying TACR2 in Genebass, selecting HbA1c as phenotype and downloading the CSV for SNP-level associations. You can find the file used in our analysis in /raw_data as genebass_file.csv_2023_02_26_19_24_47. Importantly, we did not use any filter to avoid excluding the variants  that were excluded in the burden analysis since we were interested in all SNP-level associations available. 

HUGEAMP T2D KP's HbA1c associations for HK1-NK2R-TSPAN15 loci. These can be found in /raw_data as hugeamp_associations.csv. To obtain the data go to the following link and download the CSV in the section "Variants in region (Ancestry: All)" for all ancestries: https://t2d.hugeamp.org/region.html?chr=10&end=71226674&phenotype=HBA1C&start=71113659.

## Dependencies:

All analysis were performed in Rstudio () with R (). 

The following libraries were used:

X
X
X
X

## Figures:

The raw figures are here:

