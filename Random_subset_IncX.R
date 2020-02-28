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

df <- df %>% filter(X2 != "SAL_BA7906AA_AS.scaffold.fasta" & X2 != "SAL_CA4730AA_AS.scaffold.fasta")

#Colname reassignment
colnames(df) <- df_colnames

#Convert percent coverage and identity to numeric type to allow filtering
df$perc_coverage <- as.numeric(df$perc_coverage)
df$perc_identity <- as.numeric(df$perc_identity)

#Filter to perc_coverage and perc_identity > 90%
df <- df %>% filter(perc_coverage > 90) %>% filter( perc_identity > 90)

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

#Change column names to more widely used versions
colnames(df3) <- gsub("APH\\(3''\\)-Ib","strA",colnames(df3))
colnames(df3) <- gsub("APH\\(6\\)-Id","strB",colnames(df3))
colnames(df3) <- gsub("Qnrs1","qnrs1",colnames(df3))
colnames(df3) <- gsub("TEM-1","blaTEM-1",colnames(df3))
colnames(df3) <- gsub("^CMY","blaCMY",colnames(df3))
colnames(df3) <- gsub("^CTX","blaCTX",colnames(df3))

df3$Assembly_barcode <- rownames(df3)

df3 <- df3 %>% select(Assembly_barcode, everything())

write.table(df3, "Supplementary Table 2.txt", sep = "\t", row.names = FALSE)

IncX <- df3[df3$`IncX1_1 74/80 (92%)` == 1,] 

IncX_names <- rownames(IncX)

set.seed(1)

test3 <- base::sample(IncX_names, size = 10) %>% sort()



