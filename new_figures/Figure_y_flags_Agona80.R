library(tidytree)
library(ggtree)
library(pheatmap)
library(magrittr)
library(dplyr)
library(readr)
library(reshape2)
library(ComplexHeatmap)
library(ggplot2)
library(ggimage)

### IMPORTANT ###
### Check lines where there are hardcoded changes, such as line 21,
### if you plan to modify this script's inputs.

#Replace the variable below with the path to your SG17-135 repo
path_to_repo <-
        "/Users/maxcummins/Dropbox/Doctorate/Manuscripts/Salmonella_AMR/Submission_2-mSphere/SG17-135"

#Changes working directory
setwd(path_to_repo)

path_to_tree <- "analysis/snp_outputs/HC5-4181/HC5-4181_SG17-135.clean.fullcore.tree"

path_to_abricate <- "analysis/abricate/abricate.txt"

path_to_metadata <- "Metadata/S.agona-3-6-20.txt"

path_to_cgMLST <- "Metadata/cgMLST.txt"

path_to_pais <- "analysis/PAIs_present_absent.csv"

#Changes working directory
setwd(path_to_repo)


#Read in the tree file
tree <-
        read.tree(file = path_to_tree)



#trim the names of the assemblies in the tree tip labels
tree$tip.label <- gsub("\\..*", "", tree$tip.label)

#reassign strain SG17-135's assembly barcode to its strain name (tree tip label)
tree$tip.label <- gsub("Reference", "SG17-135", tree$tip.label)

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

#Trim excess characters the assembly names and reassign this to rownames
df$name <- gsub("\\..*", "", df$name)

#Replace "SAL_HC4750AA_AS" with SG17-135
df$name <- gsub("SAL_HC4750AA_AS", "SG17-135", df$name)

df <- df %>% filter(name %in% tree$tip.label)

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

#AMR genes
df$GENE <- gsub("^AAC", "aac", df$GENE)
df$GENE <- gsub("^ANT", "ant", df$GENE)
df$GENE <- gsub("^APH","aph", df$GENE)
df$GENE <- gsub("CMY-59", "blaCMY-59", df$GENE)
df$GENE <- gsub("CTX-M-55", "blaCTX-M-55", df$GENE)
df$GENE <- gsub("TEM-1", "blaTEM-1", df$GENE)
df$GENE <- gsub("FosA7", "fosA7", df$GENE)
df$GENE <- gsub("QnrS1", "qnrS1", df$GENE)

#Plasmid genes
df$GENE <- gsub("_Gamma_1", "(Gamma)", df$GENE)
df$GENE <- gsub("_1_Alpha", "(Alpha)", df$GENE)
df$GENE <- gsub("_1.*","",df$GENE)


################################################################################
####Remove unwanted AMR genes (not considered actual AMR genes by most)
################################################################################

df <- df %>% filter(GENE != "golS")
df <- df %>% filter(GENE != "sdiA")


################################################################################
####Replace strA/B names (more commonly used names for these genes)
################################################################################

df$GENE <- gsub("aph\\(3''\\)-Ib", "strA", df$GENE)
df$GENE <- gsub("aph\\(6\\)-Id", "strB", df$GENE)


#Prepend the database name to the gene name so we can identify specific databases later
df$GENE <- paste0(df$DATABASE, "_", df$GENE)
df$GENE <- gsub(".*EC_custom_Feb_2020_","", df$GENE)

#trim db name to something smaller
df$GENE <- gsub("^vfdb", "v", df$GENE)
df$GENE <- gsub("^card", "r", df$GENE)
df$GENE <- gsub("^plasmidfinder", "p", df$GENE)

#Read in pointfinder data on resistance snps
res_snps <- read_delim("analysis/pointfinder/pointfinder_results.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

#fix colnames of res_snps for consistency
colnames(res_snps) <- gsub("name_file","name", colnames(res_snps))

#fix colnames of res_snps for consistency
colnames(res_snps) <- gsub("Mutation","GENE", colnames(res_snps))

#fix colnames of res_snps for consistency
res_snps$GENE <- gsub("^","r_", res_snps$GENE)

#fix sample names of res_snps for consistency
res_snps$name <- gsub("/.*.tsv","",res_snps$name)

#Cast the dataframe into a more usable format
res_snps2 <- dcast(data = res_snps, name ~ GENE, length, drop = FALSE)

#Replace "SAL_HC4750AA_AS" with SG17-135
res_snps2$name <- gsub("SAL_HC4750AA_AS", "SG17-135", res_snps2$name)

#Select the columns "name" and "GENE"
df2 <- df %>% select(contains("name"), contains("GENE"))

#Cast the dataframe into a more usable format
df3 <- dcast(data = df2, name ~ GENE, length, drop = FALSE)

#Bind Res_SNP data to abricate data
df3 <- left_join(df3, res_snps2)

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

#reassign strain SG17-135's assembly barcode to its strain name (data frame with genotypic data)
rownames(df3) <- gsub("SAL_HC4750AA_AS", "SG17-135", rownames(df3))

#Read in the abricate genotype data sheet (small number of rows for colname reassignment)
pais_df <- read_csv("analysis/PAIs_present_absent.csv")

colnames(pais_df)[1] <- "name"


#Replace "SAL_HC4750AA_AS" with SG17-135
pais_df$name <- gsub("SAL_HC4750AA_AS", "SG17-135", pais_df$name)

pais_df <- pais_df %>% filter(name %in% tree$tip.label)

################################################################################
####Clean column names
################################################################################

#PAI_genes
colnames(pais_df) <- gsub("_NC_003198","", colnames(pais_df))
colnames(pais_df) <- gsub("_sitABCD_AF128999","", colnames(pais_df))

################################################################################
################################################################################

#add in colname change prepend with v2_

#add DB identifier for later
colnames(pais_df)[2:ncol(pais_df)] <- paste0("v2_", colnames(pais_df)[2:ncol(pais_df)])

#Set rownames to be equal to the name column
rownames(pais_df) <- pais_df$name

#Remove the column that has assembly names
pais_df2 <- pais_df[, 2:ncol(pais_df)]

#Change multiple hits to a binary hit (gene considered present or absent, regardless of copy number)
pais_df2[pais_df2 > 1] <- 1

#Sum the number of hits per gene
colSums(pais_df2) -> sums

#reassign strain SG17-135's assembly barcode to its strain name (data frame with genotypic data)
rownames(pais_df2) <-
        gsub("SAL_HC4750AA_AS", "SG17-135", rownames(pais_df))

pais_df2$name <- pais_df$name


#reassign Colnames to include the number the format "gene_A n=41/82 (50%)"
colnames(pais_df2) <-
        paste0(colnames(pais_df2),
               "  ",
               sums,
               "/",
               nrow(pais_df2),
               " (",
               round(sums / nrow(pais_df2) * 100),
               "%)")

#reassign strain SG17-135's assembly barcode to its strain name (data frame with genotypic data)
rownames(pais_df3) <-
        gsub("SAL_HC4750AA_AS", "SG17-135", rownames(pais_df3))

df3$name <- rownames(df3)


pais_df2$name <- pais_df$name

#join the dfs of abricate data and PAI data
df3 <- left_join(df3, pais_df2)

#isolate the individual DBs so we can edit their contents for colorisation in the heatmap later
res <- df3 %>% select(starts_with("r"))
plas <- df3 %>% select(starts_with("p"))
vir <- df3 %>% select(starts_with("v_"))
vir2 <- df3 %>% select(starts_with("v2"))

#remove any remaining columns with no hits
res <- res[,colSums(res) > 0]
plas <- plas[,colSums(plas) > 0]
vir <- vir[,colSums(vir) > 0]
vir2 <- vir2[,colSums(vir2) > 0]

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
rownames(df4) <- df3$name

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

#Remove the DB prepend from df4 column names
colnames(df4) <- gsub("^[^_]+_", "", colnames(df4))

#Remove spaces from column names for metadata - this usually causes issues
colnames(Metadata) <- gsub(" ", "_", colnames(Metadata))

#Provides clade numbers to colour by.
#If you dont know what the node labels are use the line below under "get node labels"
#You may need to change the sizing..
tree2 <- groupClade(tree, c(87, 86))

Metadata$Flag <- Metadata$Country
Metadata$Flag <- gsub("Australia", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Australia.jpg", Metadata$Flag)
Metadata$Flag <- gsub("Denmark", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Denmark.jpg", Metadata$Flag)
Metadata$Flag <- gsub("United States", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/USA.png", Metadata$Flag)
Metadata$Flag <- gsub("^Ireland$", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Ireland.png", Metadata$Flag)
Metadata$Flag <- gsub("Vietnam", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Vietnam.png", Metadata$Flag)
Metadata$Flag <- gsub("United Kingdom", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/United_Kindgom.png", Metadata$Flag)
Metadata$Flag <- gsub("Mexico", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Mexico.png", Metadata$Flag)
Metadata$Flag <- gsub("Ghana", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Ghana.png", Metadata$Flag)
Metadata$Flag <- gsub("United Arab Emirates", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/UAE.png", Metadata$Flag)
Metadata$Flag <- gsub("Scotland", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Scotland.png", Metadata$Flag)
Metadata$Flag <- gsub("Germany", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Germany.png", Metadata$Flag)
Metadata$Flag <- gsub("Northern Ireland", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Northern_Ireland.jpg", Metadata$Flag)


#Generate the tree
p <- ggtree(tree2) %<+%
        Metadata +
        geom_tiplab(size = 1.25,
                    align = TRUE,
                    offset = 0.007,
                    aes(color = HC5_or_other)
        ) +
        geom_tippoint(size = 1,
                      aes(color = Source_Niche),
                      show.legend = TRUE) +
        geom_tiplab(aes(image = Flag), geom="image", size = 0.016, align = FALSE, offset = 0.0007
        ) + geom_treescale(x = 0.1, y = -2.25, offset = -2.4, fontsize = 2)#+
#get node labels
#geom_text2(size = 2, aes(subset=!isTip, label=node), hjust=-.5)


colval <- c(
        #'steelblue', #  'clade 1'                       √
        #'firebrick', # 'middle clade?' 
        #'steelblue', #  'clade 2'                       √
        '#b6933f', #      'Source Niche - Food'           √
        'black', #      'HC5_or_other (HC5-4181)        √
        '#609bce', #      'Source Niche - Human'          √
        'white', #      'Genes - 0 (N)                  √
        '#80b1d3', #      'Plasmid Genes - 1 (P)          √       
        '#bebada', #      'AMR genes - 1 - (R)            √ 
        '#fb8072', #        'VFDB genes -1 - (V)'         √
        '#cd2a18', #  'PAI genes - 1 (V2)                 √    
        '#8d68ca' #  'Source Niche - Wild Animal'        √
)

pdf(file = "new_figures/Flag_Figure_x_Agona80.pdf", paper = "a4r")

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
        offset = 0.037,
        width = 2,
        color = 'grey'
) +
        theme(legend.position = "none") +
        scale_fill_manual(
                aesthetics = c("colour", "fill"),
                values = colval,
                na.value = 'grey'
        ) +
        ggplot2::ylim(NA, 100)


dev.off()

df_by_scaff <- df %>% group_by(name, SEQUENCE) %>% dplyr::summarise(same_scaff=paste(unique(GENE), collapse=" "))

inc_by_scaff <- df_by_scaff %>% filter(grepl("Inc", same_scaff)) 

inc_by_scaff$Inc <- gsub(".* p_", "p_", inc_by_scaff$same_scaff)

spl <-strsplit(as.character(inc_by_scaff$same_scaff), " ")
genes <- data.frame(name = inc_by_scaff$name,
                        SEQUENCE = inc_by_scaff$SEQUENCE,
                        inc = inc_by_scaff$Inc,
                        gene1 = sapply(spl, "[", 1),
                        gene2 = sapply(spl, "[", 2),
                        gene3 = sapply(spl, "[", 3),
                        gene4= sapply(spl, "[", 4),
                        gene5 = sapply(spl, "[", 5),
                        gene6 = sapply(spl, "[", 6),
                        gene7 = sapply(spl, "[", 7),
                        gene8= sapply(spl, "[", 8),
                        gene9 = sapply(spl, "[", 9),
                        gene10 = sapply(spl, "[", 10),
                        gene11 = sapply(spl, "[", 11),
                        gene12= sapply(spl, "[", 12),
                        gene13 = sapply(spl, "[", 13),
                        gene14 = sapply(spl, "[", 14),
                        gene15 = sapply(spl, "[", 15),
                        gene16 = sapply(spl, "[", 16),
                        gene17 = sapply(spl, "[", 17),
                        gene18 = sapply(spl, "[", 18),
                        gene19 = sapply(spl, "[", 19),
                        gene20 = sapply(spl, "[", 20),
                        gene21 = sapply(spl, "[", 21))

genes <- melt(genes, id=1:3, value.name = "gene", na.rm = TRUE)

genes <- genes %>% select(-starts_with("variable"))

genes <- dcast(data = genes, name + inc ~ gene, length, drop = TRUE)

IncX1 <- genes %>% filter(inc == "p_IncX1", p_IncX1 == 1)

IncX1_df <- IncX1[,3:ncol(IncX1)]

IncX1_df <- IncX1_df[,colSums(IncX1_df) > 0]

rownames(IncX1_df) <- IncX1$name

pheatmap(IncX1_df)

