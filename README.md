# SG17-135
This is a repository containing scripts for the generation of figures within a publication on an XDR strain of Salmonella sourced from an Australian Gull.

DOI: **To be added post publication**

## Contents

|Directory|Contents|
|--------|------|
|analysis|Output of bioinformatic analysis tools|
|assemblies|Genome assemblies from Enterobase|
|fasta_dbs|Nucleotide databases uses in abricate|
|figures|Final figures of manuscript|
|flags|Flag images used in generation of figures|
|logs|Logs generated during bioinformatic analyses|
|metadata|Metadata on samples under analysis|
|scripts|Scripts used for post processing of data and data visualisation|
|supplemental_material|Supplementary material for the manuscripts|

## analysis
### abricate
Contains abricate output, concatenated from each resulting output file:
* abricate.txt - Contains abricate data on plasmid, virulence and AMR gene carriage
* abricate_PAIs.txt - Contains abricate data on Salmonella Pathogenicity island (SPI) and pathogenicity associated island (PAI) carriage

### snp_outputs
Contains output of snplord pipeline
* Agona195 - Contains snp matrix and phylogenetic tree for large subset of Agona genomes (n=195)
* Agona80 - Contains snp matrix and phylogenetic tree for small subset of Agona genomes closely related to SG17-135 (HC5:4181 strains) (n=80)

### pointfinder
Contains output of pointfinder analysis
* pointfinder_results.txt - lists AMR associated SNPs associated with particular genomes from the Agona195 subset

## scripts
### Figure 1.R
This script is used to generate Figure 1, as well as Supplementary Table 1 which combines Metadata and Genotypic data of samples

### Figure_2.R
This script is used to generate Figure 2, as well as Supplementary Table 4 which combines lists the co-association of AMR genes and IncX scaffolds

### Random_subset_IncX.R
This script is used to (pseudo)randomly select 10 IncX-positive strains from the Agona80 subset for analysis using BRIG

### SPI-analysis.R
This script is used to classify strains as SPI/PAI positive or negative based on >=60% (discontiguous) coverage and >=95% nucleotide identity for a given genetic element









## Majestic Seagull



<p align="center">
  <img width="600" height="600" src="https://raw.githubusercontent.com/maxlcummins/SG17-135/master/majestic_gull_dot_jpeg.png">
</p>





## R Session info
```
R version 3.6.0 (2019-04-26)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

Random number generation:
 RNG:     Mersenne-Twister 
 Normal:  Inversion 
 Sample:  Rounding 
 
locale:
[1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] reshape2_1.4.3  readr_1.3.1     dplyr_0.8.3     magrittr_1.5    pheatmap_1.0.12 ggtree_1.16.6  

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2         pillar_1.4.2       compiler_3.6.0     RColorBrewer_1.1-2 plyr_1.8.4         tools_3.6.0        zeallot_0.1.0      jsonlite_1.6      
 [9] tidytree_0.2.6     tibble_2.1.3       gtable_0.3.0       nlme_3.1-141       lattice_0.20-38    pkgconfig_2.0.2    rlang_0.4.0        rstudioapi_0.10   
[17] rvcheck_0.1.3      yaml_2.2.0         parallel_3.6.0     treeio_1.8.2       stringr_1.4.0      vctrs_0.2.0        hms_0.5.1          grid_3.6.0        
[25] tidyselect_0.2.5   glue_1.3.1         R6_2.4.0           tidyr_0.8.3        ggplot2_3.2.1      purrr_0.3.2        scales_1.0.0       backports_1.1.4   
[33] assertthat_0.2.1   ape_5.3            colorspace_1.4-1   labeling_0.3       stringi_1.4.3      lazyeval_0.2.2     munsell_0.5.0      crayon_1.3.4 
```
