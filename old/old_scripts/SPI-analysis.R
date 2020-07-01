library(pheatmap)
library(ggplot)
library(magrittr)
library(dplyr)
library(reshape2)
library(readr)
library(OneR)
library(microbenchmark)

#path_to_pais <- "/Users/maxcummins/Dropbox/Doctorate/Manuscripts/Salmonella_AMR/SG17-135/analysis/abricate/abricate_PAIs.txt"
path_to_pais <- "/Users/maxcummins/Dropbox/Doctorate/Manuscripts/Salmonella_AMR/Submission_2-mSphere/SG17-135/analysis/abricate/abricate_PAIs_CT18.txt"

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

#Filter to perc_identity > 95%
#pais_df <-
pais_df <- pais_df %>% filter(perc_identity > 95)

#Trim excess characters the assembly names and reassign this to rownames
pais_df$name <- gsub("\\..*", "", pais_df$name)

#Replace "SAL_HC4750AA_AS" with SG17-135
pais_df$name <- gsub("SAL_HC4750AA_AS", "SG17-135", pais_df$name)

pais_df$newcov <- gsub("\\/.*","", pais_df$COVERAGE)
pais_df$length_gene <- gsub(".*\\/","", pais_df$COVERAGE)

pais_df$newcov <- gsub("-",":", pais_df$newcov)

new_df <- pais_df %>% group_by(name, GENE, length_gene) %>% filter(perc_coverage > 5) %>% summarise(start =paste(sort(unique(newcov)), collapse=","), end = paste(sort(unique(newcov)), collapse=",")) #%>% filter(grepl("SPI-1_", GENE))
#new_df <- pais_df %>% group_by(name, GENE, length_gene) %>% summarise(start =paste(sort(unique(newcov)), collapse=","), end = paste(sort(unique(newcov)), collapse=",")) #%>% filter(grepl("SPI-1_", GENE))

new_df$end <- gsub("[0-9]+:","", new_df$end)
new_df$start <- gsub(":[0-9]+","", new_df$start)


spl <-strsplit(as.character(new_df$start), ",")
start_coord <- data.frame(name = new_df$name, gene = new_df$GENE,
           length_gene = new_df$length_gene,
           chunk1 = sapply(spl, "[", 1),
           chunk2 = sapply(spl, "[", 2),
           chunk3 = sapply(spl, "[", 3),
           chunk4= sapply(spl, "[", 4),
           chunk5 = sapply(spl, "[", 5),
           chunk6 = sapply(spl, "[", 6),
           chunk7 = sapply(spl, "[", 7),
           chunk8 = sapply(spl, "[", 8))

start_coord <- melt(start_coord, id=1:3, value.name = "start")

start_coord <- start_coord %>% select(-starts_with("variable")) 

spl <-strsplit(as.character(new_df$end), ",")
end_coord <- data.frame(name = new_df$name, gene = new_df$GENE,
                          length_gene = new_df$length_gene,
                          chunk1 = sapply(spl, "[", 1),
                          chunk2 = sapply(spl, "[", 2),
                          chunk3 = sapply(spl, "[", 3),
                          chunk4= sapply(spl, "[", 4),
                          chunk5 = sapply(spl, "[", 5),
                          chunk6 = sapply(spl, "[", 6),
                          chunk7 = sapply(spl, "[", 7),
                          chunk8 = sapply(spl, "[", 8))

end_coord <- melt(end_coord, id=1:3, value.name = "end")

end_coord <- end_coord %>% select(-starts_with("variable")) 

coords <- start_coord

coords$end <- end_coord$end

coords <- coords[complete.cases(coords),]

unique(coords$length_gene)

coords$start <- as.numeric(coords$start)
coords$end <- as.numeric(coords$end)
coords$length_gene <- as.numeric(levels(coords$length_gene))[coords$length_gene]


coords$percentage <-  (((coords$end-coords$start)+1)/coords$length_gene)*100

test <- coords# %>% filter(name == "SAL_AB7542AA_AS", gene == "SPI-12_NC_006905_P4") %>% arrange(desc(end))

list_ <- vector(mode = "list", length = 0)

for(sample in unique(test$name)){
        test2 <- test %>% filter(name == sample)
        for(gene_ in unique(test$gene)){
                test3 <- test2 %>% filter(gene == gene_)
                length_of_gene <- test3$length_gene[1]
                if(is.na(length_of_gene) == FALSE){
                        range_matrix <- rep(0, times = length_of_gene)
                        for(hit in 1:nrow(test3)){
                                start_ <- test3[hit,4]
                                end_ <- test3[hit,5]
                                range_matrix[start_:end_] <- 1
                                range_matrix[range_matrix > 1] <- 1
                                }
                        }
                newline <- c(sample, gene_, round((sum(range_matrix)/length_of_gene)*100, digits = 3))
                list_ <- append(list_,newline)
        }

}

abc <- length(list_)/3

df <- data.frame(matrix(unlist(list_), nrow = abc, byrow=T), stringsAsFactors = F)

colnames(df) <- c("name","GENE","Coverage_percentage")

df$Coverage_percentage[is.na(df$Coverage_percentage)] <- 0

df$Coverage_percentage <- as.numeric(df$Coverage_percentage)

final_table <- dcast(df, name ~ GENE)

final_final_table <- final_table[1:nrow(final_table),2:ncol(final_table)]

final_final_table_2 <- final_final_table

final_final_table[final_final_table < 60] <- 0
final_final_table[final_final_table >= 60] <- 1

rownames(final_final_table) <- final_table$name

final_table <- final_final_table

write.csv(final_table, "analysis/PAIs_present_absent.csv")

pheatmap(final_final_table, fontsize_row = 2)

