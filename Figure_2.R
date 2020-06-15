library(ggtree)
library(pheatmap)
library(magrittr)
library(dplyr)
library(readr)
library(reshape2)
library(ComplexHeatmap)
library(ggplot2)

### Check lines where there are hardcoded changes, such as line 21, if you plan to modify this script's inputs.

#Replace the variable below with the path to your SG17-135 repo
path_to_repo <- "/Users/maxcummins/Dropbox/Doctorate/Manuscripts/Salmonella_AMR/Submission_2-mSphere/New_analysis/SG17-135"

#Changes working directory
setwd(path_to_repo)

#Read in the abricate genotype data sheet (small number of rows for colname reassignment)
df <- read_delim("analysis/abricate/abricate.txt", 
                 "\t", escape_double = FALSE, trim_ws = TRUE, n_max = 10)

#Colname reassignment
colnames(df)[c(1,10:11)] <- c("name","perc_coverage", "perc_identity")
df_colnames <- colnames(df)

#Read in full abricate genotype data sheet
df <- read_delim("analysis/abricate/abricate.txt", 
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
df <- df %>% filter(DATABASE %in% c("vfdb", "card", "plasmidfinder"))

df$GENE <- paste0(df$DATABASE,"_",df$GENE)

df$GENE <- gsub("^vfdb","v", df$GENE)
df$GENE <- gsub("^card","r", df$GENE)
df$GENE <- gsub("^plasmidfinder","p", df$GENE)

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
tree <- read.tree(file = "analysis/snp_outputs/fasttree/SG17-135.clean.fullcore.tree")

#trim the names of the assemblies in the tree tip labels
tree$tip.label <- gsub("\\..*","",tree$tip.label)

#reassign strain SG17-135's assembly barcode to its strain name (tree tip label)
tree$tip.label <- gsub("Reference","SG17-135",tree$tip.label)

#reassign strain SG17-135's assembly barcode to its strain name (data frame with genotypic data)
rownames(df3) <- gsub("SAL_HC4750AA_AS", "SG17-135", rownames(df3))

res <- df3 %>% select(starts_with("r"))
plas <- df3 %>% select(starts_with("p"))
vir <- df3 %>% select(starts_with("v"))

res[res == 1] <- "R"
plas[plas == 1] <- "P"
vir[vir == 1] <- "V"

df4 <- cbind(res,plas,vir)

rownames(df4) <- rownames(df3)

df4[df4 == 0] <- "N"

#Read in metadata
Metadata <- read_delim("S.agona-3-6-20.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)

cgMLST <- read_delim("cgMLST.txt", 
           "\t", escape_double = FALSE, trim_ws = TRUE)

cgMLST <- cgMLST %>% select(Uberstrain, starts_with("HC"))

Metadata <- left_join(Metadata, cgMLST)

Metadata$working_name <- Metadata$`Assembly Barcode`

Metadata$working_name <- gsub("SAL_HC4750AA_AS", "SG17-135", Metadata$working_name)

Metadata <- Metadata %>% filter(working_name %in% tree$tip.label)

#reassign column name
#colnames(Metadata)[3] <- "name"

#ComplexHeatmap::Heatmap(as.matrix(df3), show_row_names = FALSE, column_names_gp = gpar(fontsize = 3))

#Define colors for gene-type dependent coloring of gene hits
colorgenotype <- c("N" = "white", "R" = "#bebada", "V" = "#fb8072", "P" = "#80b1d3")

df4 <- df4 %>% select(starts_with("r_"),starts_with("p_"))

#Generate the tree
p <- ggtree(tree) %<+% 
  df3 +
  geom_tiplab(size = 2, align = TRUE)

#Generate the heatmap
a <- gheatmap(p = p, data = df4,
              colnames_offset_y = -0.1,
              font.size = 1.8,
              hjust = 0,
              colnames_position = "top",
              colnames_angle = 90,
              offset = 0.10,
              width = 1) + theme(legend.position = "none") +
  scale_fill_manual(values = colorgenotype, na.value = 'grey')
#+ ggplot2::ylim(NA, 88)

#Display the heatmap
a




