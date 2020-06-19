library(tidytree)
library(ggtree)
library(pheatmap)
library(magrittr)
library(dplyr)
library(readr)
library(reshape2)
library(ComplexHeatmap)
library(ggplot2)

### IMPORTANT ###
### Check lines where there are hardcoded changes, such as line 21,
### if you plan to modify this script's inputs.

#Replace the variable below with the path to your SG17-135 repo
path_to_repo <-
        "/Users/maxcummins/Dropbox/Doctorate/Manuscripts/Salmonella_AMR/Submission_2-mSphere/SG17-135"

#Changes working directory
setwd(path_to_repo)

path_to_tree <- "analysis/snp_outputs/"

path_to_abricate <- "analysis/abricate/abricate.txt"

path_to_metadata <- "Metadata/S.agona-3-6-20.txt"

path_to_cgMLST <- "Metadata/cgMLST.txt"

path_to_pais <- "/Users/maxcummins/Dropbox/Doctorate/Manuscripts/Salmonella_AMR/SG17-135/analysis/abricate/abricate_PAIs.txt"

#Changes working directory
setwd(path_to_repo)

#Read in the abricate genotype data sheet (small number of rows for colname reassignment)
df <- read_delim(
        path_to_abricate,
        "\t",
        escape_double = FALSE,
        trim_ws = TRUE,
        n_max = 10
)

#Colname reassignment
colnames(df)[c(1, 10:11)] <-
        c("name", "perc_coverage", "perc_identity")
df_colnames <- colnames(df)

#Read in full abricate genotype data sheet
df <- read_delim(
        path_to_abricate,
        "\t",
        escape_double = FALSE,
        trim_ws = TRUE,
        col_names = FALSE,
        skip = 1
)

#Remove cases where there are multiple headers from concatenation of abricate reports
df <- df %>% filter(X2 != "SEQUENCE")

#Colname reassignment
colnames(df) <- df_colnames

#Convert percent coverage and identity to numeric type to allow filtering
df$perc_coverage <- as.numeric(df$perc_coverage)
df$perc_identity <- as.numeric(df$perc_identity)

#Filter to perc_coverage and perc_identity > 90%
df <-
        df %>% filter(perc_coverage > 90) %>% filter(perc_identity > 90)

#Filter out non-virulence-associated genes
#df <-
#  df %>% filter(DATABASE %in% c("vfdb", "card", "plasmidfinder"))

#Remove SPI associated VFs
df <- df[!grepl("SPI", df$PRODUCT), ]

################################################################################
####Remove columns that we dont need
################################################################################

#Change csg genes to be operons where they co-occur
df <- df[!grepl("csg(B|C|D|E|F|G)", df$PRODUCT), ]
df$GENE <- gsub("csgA", "csgABCDEFG", df$GENE)

#Change fim genes to be operons where they co-occur
df <- df[!grepl("fim(D|F|H|I)", df$PRODUCT), ]
df$GENE <- gsub("fimC", "fimCDFHI", df$GENE)

#Change pap genes to be operons where they co-occur
df <- df[!grepl("pap(D|E|F|G|H|J|K)", df$PRODUCT), ]
df$GENE <- gsub("papC", "papCDEFGHJK", df$GENE)


################################################################################
####Clean gene names
################################################################################

df$GENE <- gsub("AAC(3)-IId", "aac(3)-IId", df$GENE)
df$GENE <- gsub("AAC(3)-IV", "aac(3)-IV", df$GENE)
df$GENE <- gsub("AAC(3)-VIa", "aac(3)-VIa", df$GENE)
df$GENE <- gsub("AAC(6')-Iy", "aac(6')-Iy", df$GENE)
df$GENE <- gsub("ANT(3'')-IIa", "ant(3'')-IIa", df$GENE)
df$GENE <- gsub("APH(3'')-Ib", "aph(3'')-Ib", df$GENE)
df$GENE <- gsub("APH(3')-Ia", "aph(3')-Ia", df$GENE)
df$GENE <- gsub("APH(3')-IIa", "aph(3')-IIa", df$GENE)
df$GENE <- gsub("APH(4)-Ia", "aph(4)-Ia", df$GENE)
df$GENE <- gsub("APH(6)-Id", "aph(6)-Id", df$GENE)
df$GENE <- gsub("CMY-59", "blaCMY-59", df$GENE)
df$GENE <- gsub("CTX-M-55", "blaCTX-M-55", df$GENE)
df$GENE <- gsub("TEM-1", "blaTEM-1", df$GENE)
df$GENE <- gsub("TEM-150", "blaTEM-150", df$GENE)
df$GENE <- gsub("FosA7", "fosA7", df$GENE)
df$GENE <- gsub("QnrS1", "qnrS1", df$GENE)

################################################################################
####Remove unwanted AMR genes (not considered actual AMR genes by most)
################################################################################

df <- df %>% filter(GENE != "golS")
df <- df %>% filter(GENE != "sdiA")


################################################################################
################################################################################






#Prepend the database name to the gene name so we can identify specific databases later
df$GENE <- paste0(df$DATABASE, "_", df$GENE)

#trim db name to something smaller
df$GENE <- gsub("^vfdb", "v", df$GENE)
df$GENE <- gsub("^card", "r", df$GENE)
df$GENE <- gsub("^plasmidfinder", "p", df$GENE)

#Select the columns "name" and "GENE"
df2 <- df %>% select(contains("name"), contains("GENE"))

#Cast the dataframe into a more usable format
df3 <- dcast(data = df2, name ~ GENE, length, drop = FALSE)

#Replace NAs (indicating a lack of a gene hit) with a zero
df3[is.na(df3)] <- 0

#Trim excess characters the assembly names and reassign this to rownames
rownames(df3) <- gsub("\\..*", "", df3$name)

#Remove the column that has assembly names
df3 <- df3[, 2:ncol(df3)]

#Change multiple hits to a binary hit (gene considered present or absent, regardless of copy number)
df3[df3 > 1] <- 1

#Sum the number of hits per gene
colSums(df3) -> sums

#reassign Colnames to include the number the format "gene_A 50% (n=41/82)"
colnames(df3) <-
        paste0(colnames(df3),
               " ",
               sums,
               "/",
               nrow(df3),
               " (",
               round(sums / nrow(df3) * 100),
               "%)")

#Read in the tree file
tree <-
        read.tree(file = path_to_tree)

#trim the names of the assemblies in the tree tip labels
tree$tip.label <- gsub("\\..*", "", tree$tip.label)

#reassign strain SG17-135's assembly barcode to its strain name (tree tip label)
tree$tip.label <- gsub("Reference", "SG17-135", tree$tip.label)

#reassign strain SG17-135's assembly barcode to its strain name (data frame with genotypic data)
rownames(df3) <- gsub("SAL_HC4750AA_AS", "SG17-135", rownames(df3))

#Read in the abricate genotype data sheet (small number of rows for colname reassignment)
pais_df <-
        read_delim(
                path_to_pais,
                "\t",
                escape_double = FALSE,
                trim_ws = TRUE,
                n_max = 10
        )

#Colname reassignment
colnames(pais_df)[c(1, 10:11)] <-
        c("name", "perc_coverage", "perc_identity")
pais_df_colnames <- colnames(pais_df)

#Re-read in PAI abricate genotype data sheet
pais_df <-
        read_delim(
                path_to_pais,
                "\t",
                escape_double = FALSE,
                trim_ws = TRUE,
                col_names = FALSE,
                skip = 1
        )

#Remove cases where there are multiple headers from concatenation of abricate reports
pais_df <- pais_df %>% filter(X2 != "SEQUENCE")

#Colname reassignment
colnames(pais_df) <- pais_df_colnames

#Convert percent coverage and identity to numeric type to allow filtering
pais_df$perc_coverage <- as.numeric(pais_df$perc_coverage)
pais_df$perc_identity <- as.numeric(pais_df$perc_identity)

#Filter to perc_coverage to > 50% and perc_identity > 90%
pais_df <-
        pais_df %>% filter(perc_coverage > 50) %>% filter(perc_identity > 90)

#add DB identifier for later
pais_df$GENE <- paste0(pais_df$DATABASE, "_", pais_df$GENE)

#Trims DB identifier
pais_df$GENE <- gsub("^Salmonella_PIAs", "v2", pais_df$GENE)

#Select the columns "name" and "GENE"
pais_df2 <- pais_df %>% select(contains("name"), contains("GENE"))

#Cast the dataframe into a more usable format
pais_df3 <-
        dcast(data = pais_df2, name ~ GENE, length, drop = FALSE)

#Replace NAs (indicating a lack of a gene hit) with a zero
pais_df3[is.na(pais_df3)] <- 0

#Trim excess characters the assembly names and reassign this to rownames
rownames(pais_df3) <- gsub("\\..*", "", pais_df3$name)

#Remove the column that has assembly names
pais_df3 <- pais_df3[, 2:ncol(pais_df3)]

#Change multiple hits to a binary hit (gene considered present or absent, regardless of copy number)
pais_df3[pais_df3 > 1] <- 1

#Sum the number of hits per gene
colSums(pais_df3) -> sums

#reassign Colnames to include the number the format "gene_A n=41/82 (50%)"
colnames(pais_df3) <-
        paste0(colnames(pais_df3),
               " ",
               sums,
               "/",
               nrow(pais_df3),
               " (",
               round(sums / nrow(pais_df3) * 100),
               "%)")

#reassign strain SG17-135's assembly barcode to its strain name (data frame with genotypic data)
rownames(pais_df3) <-
        gsub("SAL_HC4750AA_AS", "SG17-135", rownames(pais_df3))

#join the dfs of abricate data and PAI data
df3 <- cbind(df3, pais_df3)

#isolate the individual DBs so we can edit their contents for colorisation in the heatmap later
res <- df3 %>% select(starts_with("r"))
plas <- df3 %>% select(starts_with("p"))
vir <- df3 %>% select(starts_with("v_"))
vir2 <- df3 %>% select(starts_with("v2"))

#DB specific cell content changing for colorisation in the heatmap later
res[res == 1] <- "R"
plas[plas == 1] <- "P"
vir[vir == 1] <- "V"
vir2[vir2 == 1] <- "V2"

#sort columnanmes alphabetically
res <- res %>% select(sort(names(.)))
plas <- plas %>% select(sort(names(.)))
vir <- vir %>% select(sort(names(.)))
vir2 <- vir2 %>% select(sort(names(.)))

#bind all the DBs again together
df4 <- cbind(res, plas, vir, vir2)

#asign rownames again
rownames(df4) <- rownames(df3)

#replace 0's with Ns
df4[df4 == 0] <- "N"

#Read in metadata
Metadata <- read_delim(path_to_metadata,
                       "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE)

#Read in cgMLST data
cgMLST <- read_delim(path_to_cgMLST,
                     "\t",
                     escape_double = FALSE,
                     trim_ws = TRUE)

#Select HC columns from cgMLST
cgMLST <- cgMLST %>% select(Uberstrain, starts_with("HC"))

#bind Metadata and cgMLST tables
Metadata <- left_join(Metadata, cgMLST)

#Add new column working_name
Metadata$working_name <- Metadata$`Assembly Barcode`

#replace the assembly barcode with strain name
Metadata$working_name <-
        gsub("SAL_HC4750AA_AS", "SG17-135", Metadata$working_name)

#Filter metadata table to only contain strains from the tree
Metadata <- Metadata %>% filter(working_name %in% tree$tip.label)

#Change the column order to put working_name first
Metadata <- Metadata %>% select(working_name, everything())

#set rowname to be equal to working name
rownames(Metadata) <- Metadata$working_name

#Create a list of colnames for HCC columns from cgMLST table
cols <- Metadata %>% select(starts_with("HC")) %>% colnames()

#reassign the above columns as factors
Metadata[cols] <- lapply(Metadata[cols], factor)

#Designate strains as HC5-4181 or Other
Metadata$HC5_or_other <- gsub("^4181$*", "HC5-4181", Metadata$HC5)
Metadata$HC5_or_other <-
        gsub("^[0-9].*", "Other", Metadata$HC5_or_other)

#Replace ND's in Source Niche with NA
Metadata$`Source Niche` <-
        gsub("ND", NA, Metadata$`Source Niche`)

#Original colour scheme for heatmap 
#Define colors for gene-type dependent coloring of gene hits
#colorgenotype <-
#  c(
#    "N" = "white",
#    "R" = "#bebada",
#    "V2" = "black",
#    "V" = "#fb8072",
#    "P" = "#80b1d3"
#  )

df4 <- df4 %>%
        #select(starts_with("v_"),starts_with("v2_"))
        select(starts_with("r_"),
               starts_with("p_"),
               starts_with("v_"),
               starts_with("v2_"))

colnames(Metadata) <- gsub(" ", "_", colnames(Metadata))

tree2 <- groupClade(tree, c(207, 206))

#Generate the tree
p <- ggtree(tree2, aes(color = group)) %<+%
        Metadata +
        geom_tiplab(size = 1.5, align = TRUE, aes(color = HC5_or_other)) +
        geom_tippoint(size = 1,
                      aes(color = Source_Niche),
                      show.legend = TRUE) #+ geom_text2(size = 2, aes(subset=!isTip, label=node), hjust=-.5)
#get node labels
#

colval <- c(
        'steelblue', #  'clade 1'
        #'steelblue', #  'middle clade?'
        'firebrick', #  'clade 2'
        '#8d68ca', #  'Source Niche - Environment'
        '#b6933f', #  'Source Niche - Food'
        'black', #  'tip labels - HC4181'
        '#609bce', #  'Source Niche - Human'
        '#64a758', #  'Source Niche - Livestock'
        'white', #  'Gene absent - N'
        #'green', #  'Source Niche - ND' - removed as ND was replaced with NA
        'gray47', #  'tip labels 2 - Other'    
        '#80b1d3', #  'Plasmid genes - P'    
        '#c85994', #  'Source Niche - Poultry'    
        '#bebada', #  'AMR genes - R'    
        '#fb8072', #  'vfdb genes - V'
        '#cd2a18', #  'SPI genes - V2'
        '#8d68ca' #  "Source Niche - Wild Animal"
)

colobject <- c(
        'clade 1',
        #'middle clade?',
        'clade 2',
        'Source Niche - Environment',
        'Source Niche - Food',
        'tip labels - HC4181',
        'Source Niche - Human',
        'Source Niche - Livestock',
        'Gene absent - N',
        #'Source Niche - ND', - removed as ND was replaced with NA
        'tip labels 2 - Other',
        'Plasmid genes - P',
        'Source Niche - Poultry',
        'AMR genes - R',
        'vfdb genes - V',
        'SPI genes - V2',
        "Source details - Wild Animal"
)

colsdf <- data.frame(colobject, colval)

#Generate the heatmap
gheatmap(
        p = p,
        data = df4,
        colnames_offset_y = -0.1,
        font.size = 1.8,
        hjust = 0,
        colnames_position = "top",
        #colnames = FALSE,
        colnames_angle = 90,
        offset = 0.0125,
        width = 2,
        color = 'grey'
) +
        theme(legend.position = "none") +
        scale_fill_manual(
                aesthetics = c("colour", "fill"),
                values = colval,
                na.value = 'grey'
        ) +
        ggplot2::ylim(NA, 230)

#Display the heatmap
a

#p <- ggtree(tree2, aes(color = group)) %<+%
#  Metadata #+
##scale_fill_manual(values=c("firebrick","steelblue",rainbow(n = 10)))
#
#p +
#  geom_tiplab(size = 1, align = TRUE, aes(color = HC5_or_other)) +
#  geom_tippoint(size = 0.5,
#                aes(color = Source_Niche),
#                show.legend = TRUE) +
#  geom_text2(aes(subset = !isTip, label = node),
#             hjust = -.3,
#             size = 1) +
#  scale_fill_manual(values = c(rainbow(n = 10), "firebrick", "steelblue"))
#