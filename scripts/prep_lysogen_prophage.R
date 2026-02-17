library(dplyr)

setwd("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/input_data/")

# Genomad -----------------------------------------------------------------
gen_vib <- read.delim("votu_clusters.tsv", header = F)
colnames(gen_vib) <- c("Virus", "members")
gen_vib <- gen_vib %>% mutate(Genomad = case_when(grepl("provirus", gen_vib$members) ~ "provirus", 
                                                  .default = "lytic"),
                              VIBRANT = case_when(grepl("fragment", gen_vib$members) ~ "provirus", 
                                                  .default = "lytic")) %>%
  select(-members)

gen_vib[gen_vib$Virus == "MAH_S4_1189", "Virus"] <- "MAH_S4_1189_1_checkv_1-6616/21217"
gen_vib[gen_vib$Virus == "MAH_S4_267", "Virus"] <- "MAH_S4_267_1_checkv_1-17879/48598"
gen_vib[gen_vib$Virus == "POI_B4_117", "Virus"] <- "POI_B4_117_1_checkv_1-4335/63936"
add <- data_frame(Virus = c("MAH_S4_1189_2_checkv_16880-21217/21217", 
                            "MAH_S4_267_2_checkv_34651-48598/48598", 
                            "POI_B4_117_2_checkv_30486-63936/63936"),
                  Genomad = rep("lytic", 3),
                  VIBRANT = rep("lytic", 3))
gen_vib <- rbind(gen_vib, add)

vibrant <- read.delim("vibrant_votus_genome_quality.tsv", header = F, col.names = c("Virus", "life", "quality"))
vibrant <- unique(vibrant[ ,1:2])
gen_vib <- merge(gen_vib, vibrant, by = "Virus", all.x = TRUE)
gen_vib <- gen_vib %>% mutate(VIBRANT = case_when(VIBRANT == "lytic" & life == "lysogenic" ~ "lysogenic",
                                       .default = VIBRANT)) %>%
  select(-life)

# CheckV ------------------------------------------------------------------
checkv <- read.delim("gv_votus_checkv/checkv_provirus_names.txt", header = F, col.names = c("name"))
checkv$Virus <- gsub("_1_checkv.*", "", checkv$name)

gen_vib <- merge(gen_vib, checkv, by = "Virus", all.x = TRUE)

gen_vib$Virus <- ifelse(is.na(gen_vib$name), 
                        gen_vib$Virus, 
                        gen_vib$name)
gen_vib <- gen_vib %>% select(-name) %>% 
  mutate(CheckV = case_when(grepl("checkv", Virus) ~ "provirus", .default = "lytic"))

# Metacerberus ------------------------------------------------------------
metacerb <- read.delim("annotations-prodigal_viruses_2p5kb_1vg/final_annotation_summary.tsv", quote = "")
integ <- metacerb[grepl("integrase", metacerb$product, ignore.case = T), c("target")]
recomb <- metacerb[grepl("recombinase", metacerb$product, ignore.case = T), c("target")]
transpo <- metacerb[grepl("transposase", metacerb$product, ignore.case = T), c("target")]
excis <- metacerb[grepl("excisionase", metacerb$product, ignore.case = T), c("target")]
ci_cro <- metacerb[grepl(" cro[/ ]", metacerb$product, ignore.case = T), c("target")]
rep <- metacerb[grepl("repressor", metacerb$product, ignore.case = T), ]
rep <- rep[grepl("phage", rep$product, ignore.case = T), c("target")]
par <- metacerb[grepl("ParB|ParA", metacerb$product), ]
par <- par[par$HMM %in% c("PHROG", "VOG"), c("target")]
cii <- metacerb[grepl("Phage regulatory protein CII", metacerb$product, ignore.case = T), c("target")]

prophage_genes <- unique(c(integ, recomb, transpo, excis, ci_cro, rep, par, cii))
prophages <- unique(gsub("_[^_]*$", "", prophage_genes))
prophage_df <- as.data.frame(prophages)
prophage_df <- prophage_df %>% 
  mutate(fixed = case_when(grepl("_[0-9]*_[0-9]", prophages) & !grepl("provirus|metator", prophages) ~ gsub("_[^_]*$", "", prophages),
                           .default = prophages), 
         .keep = "unused")
prophage_df$`Lysogenic Genes` <- "lysogenic"
colnames(prophage_df)[1] <- "Virus"

# BLAST -------------------------------------------------------------------
blast <- read.delim("prophage_quality_bins_no_derep_filtered_blastn_95ID_100cov_500bpflank.tsv", header = F)
blast <- unique(data.frame(Virus = unique(blast[ , 1])))
blast$`BLAST MAGs` <- "provirus"

# Merged Table ------------------------------------------------------------

final_table <- merge(gen_vib, prophage_df, by = "Virus", all.x = T) %>%
  merge(blast, by = "Virus", all.x = T) 


final_table[is.na(final_table)] <- "lytic"

final_table <- final_table %>% rowwise() %>%
  mutate(n_lysogen_ids = sum(c_across(Genomad:`BLAST MAGs`) == "lysogenic"),
         n_provirus_ids = sum(c_across(Genomad:`BLAST MAGs`) == "provirus"),
         .after = 1)

write.table(final_table, "~/Documents/SPACE_2023/viruses/figures/data/viral_lifestyle_summary.tsv",
           sep = "\t", row.names = F, quote = F)                       
