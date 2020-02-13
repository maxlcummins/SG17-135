library(ggtree)
library(pheatmap)
library(magrittr)
library(dplyr)
library(readr)
library(reshape2)

#Replace SG17-135 on the line below to the github repo containing this script, then unhash the line below
#setwd("SG17-135")

#Read in the abricate genotype data sheet (small number of rows for colname reassignment)
df <- read_delim("abricate_HC5-4181.txt", 
                 "\t", escape_double = FALSE, trim_ws = TRUE, n_max = 10)

#Colname reassignment
colnames(df)[c(1,9:10)] <- c("name","perc_coverage", "perc_identity")
df_colnames <- colnames(df)

#Read in full abricate genotype data sheet
df <- read_delim("abricate_HC5-4181.txt", 
                 "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE, skip = 1)

#Remove cases where there are multiple headers from concatenation of abricate reports
df <- df %>% filter(X2 != "SEQUENCE")

#Colname reassignment
colnames(df) <- df_colnames

#Convert percent coverage and identity to numeric type to allow filtering
df$perc_coverage <- as.numeric(df$perc_coverage)
df$perc_identity <- as.numeric(df$perc_identity)

#Filter to perc_coverage and perc_identity > 90%
df <- df %>% filter(perc_coverage > 90) %>% filter( perc_identity > 90)

#Filter out non-virulence-associated genes
df <- df %>% filter(DATABASE %in% c("vfdb"))

#Select the columns "name" and "GENE"
df2 <- df %>%select(contains("name"), contains("GENE"))

#Cast the dataframe into a more usable format
df3 <- dcast(data = df2, name ~ GENE, length, drop = FALSE)

#Replace NAs (indicating a lack of a gene hit) with a zero
df3[is.na(df3)] <- 0

#Trim excess characters the assembly names and reassign this to rownames
rownames(df3) <- gsub("\\..*","",df3$name)

#Remove the column that has assembly names
df3 <- df3[,2:ncol(df3)]

#Change multiple hits to a binary hit (gene considered present or absent, regardless of copy number)
df3[df3 > 1] <- 1

#Sum the number of hits per gene
colSums(df3) -> sums

#reassign Colnames to include the number the format "gene_A 50% (n=41/82)"
colnames(df3) <- paste0(colnames(df3), " ", sums, "/", nrow(df3), " (", round(sums/nrow(df3)*100), "%)")

#Read in the tree file
tree <- read.tree(file = "SG17-135.clean.fullcore.tree")

#trim the names of the assemblies in the tree tip labels
tree$tip.label <- gsub("\\..*","",tree$tip.label)

#reassign strain SG17-135's assembly barcode to its strain name (tree tip label)
tree$tip.label <- gsub("SAL_PB9279AA_AS","SG17-135",tree$tip.label)

#reassign strain SG17-135's assembly barcode to its strain name (data frame with genotypic data)
rownames(df3) <- gsub("SAL_PB9279AA_AS", "SG17-135", rownames(df3))

#reassign column name
colnames(Metadata_HC5_4181)[3] <- "name"

#Generate the tree
p <- ggtree(tree) %<+% 
  df3 +
  geom_tiplab(size = 2, align = TRUE) 

#Generate the heatmap
a <- gheatmap(p = p, data = df3,
              colnames_offset_y = -0.1,
              low = "white",
              high = "red",
              font.size = 1.8,
              hjust = 0,
              colnames_position = "top",
              colnames_angle = 90,
              offset = 0.10,
              width = 6) + ggplot2::ylim(NA, 88) + theme(legend.position = "none")

#Display the heatmap
a




