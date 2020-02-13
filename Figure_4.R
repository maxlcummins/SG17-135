library(ggtree)
library(pheatmap)
library(magrittr)
library(dplyr)
library(readr)
library(reshape2)

df <- read_delim("~/Dropbox/Doctorate/Manuscripts/Salmonella_AMR/Enterobase/Assemblies_HC5-4181/abricate_HC5-4181.txt", 
                 "\t", escape_double = FALSE, trim_ws = TRUE)

colnames(df)[c(1,9:10)] <- c("name","perc_coverage", "perc_identity")

df_colnames <- colnames(df)

df <- read_delim("~/Dropbox/Doctorate/Manuscripts/Salmonella_AMR/Enterobase/Assemblies_HC5-4181/abricate_HC5-4181.txt", 
                 "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE, skip = 1)

df <- df %>% filter(X2 != "SEQUENCE")

colnames(df) <- df_colnames

df$perc_coverage <- as.numeric(df$perc_coverage)

df$perc_identity <- as.numeric(df$perc_identity)

df <- df %>% filter(perc_coverage > 90) %>% filter( perc_identity > 90)

df <- df %>% filter(DATABASE %in% c("card"))

df2 <- df %>%select(contains("name"), contains("GENE"))

df3 <- dcast(data = df2, name ~ GENE, length, drop = FALSE)

df3[is.na(df3)] <- 0

rownames(df3) <- gsub("\\..*","",df3$name)

df3 <- df3[,2:ncol(df3)]

df3$FQR <- 1

df3[df3 > 1] <- 1

colSums(df3) -> sums

colnames(df3) <- paste0(colnames(df3), " ", sums, "/", nrow(df3), " (", round(sums/nrow(df3)*100), "%)")

tree <- read.tree(file = "~/Dropbox/Doctorate/Manuscripts/Salmonella_AMR/Enterobase/HC5-4181_analysis/SG17-135.clean.fullcore.tree")

tree$tip.label <- gsub("\\..*","",tree$tip.label)
tree$tip.label <- gsub("SAL_PB9279AA_AS","SG17-135",tree$tip.label)

rownames(df3) <- gsub("SAL_PB9279AA_AS", "SG17-135", rownames(df3))

colnames(Metadata_HC5_4181)[3] <- "name"

tips <- as.data.frame(tree$tip.label)

colnames(tips) <- "name"


p <- ggtree(tree) %<+% 
  df3 +
  geom_tiplab(size = 2, align = TRUE) 

p

a <- gheatmap(p = p, data = df3,
              colnames_offset_y = -0.1,
              low = "white",
              high = "blue",
              font.size = 2,
              hjust = 0,
              colnames_position = "top",
              colnames_angle = 90,
              offset = 0.04,
              width = 2) + ggplot2::ylim(NA, 88) + theme(legend.position = "none")

a




