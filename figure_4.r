library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(plotly)
library(vegan)
library(pairwiseAdonis)
library(ggpubr)
library(rstatix)
library(scales)
library(circlize)
library(ggpattern)
library(ggnewscale)
library(ggsci)
library(paletteer)
library(plotrix)

pal <- c("Carbon fixation by Calvin cycle" = "#120e99", 
         "Other carbon fixation pathways" = "#377EB8",
         "Carotenoid biosynthesis" = "#ffb15e",
         "Photosynthesis" = "#FC8D62", 
         "Photosynthesis - antenna proteins" = "#f0cdb1", 
         "Porphyrin metabolism" = "#009476",
         "Glycolysis / Gluconeogenesis" = "#bf5328",
         "Sulfur metabolism" = "#F781BF",
         "Sulfur relay system" = "#dba1f7", 
         "Pentose phosphate pathway" = "#8DA0CB",
         "Citrate cycle (TCA cycle)" = "#63e0dc")
# Initial Data Importing and Wrangling ------------------------------------
setwd("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/input_data/")
#Import metacerberus annotations
metacerb <- read.table("annotations-prodigal_viruses_2p5kb_1vg/final_annotation_summary.tsv",
                       header = T, sep = "\t", quote = "")
#Create column for Viral contig association
metacerb$Virus <- gsub("_[^_]*$", "", metacerb$target)
#Create column for extracting gene orders
metacerb$gene_num <- gsub("^.*_", "", metacerb$target)
#Read in the kegg annotations from metacerberus
kegg <- read.table("annotations-prodigal_viruses_2p5kb_1vg/annotation_summary_KOFam_all_KEGG.tsv",
                    header = T, sep = "\t", quote = "")
#Create list of viral databases for identifying viral genes
vdbs <- c("VOG", "PHROG", "GVDB", "PVOG")
#Create list of phage specific genes
phage_hit <- metacerb[grepl("phage", metacerb$product, ignore.case = T), "best_hit"] 
#I manually checked the above list and decided to exclude 2 genes
phage_hit <- phage_hit[!phage_hit %in% c("COG5340", "antiphage_deaminase")]
#Create list of virus specific genes from PFAM
pfam_vir <- metacerb[grepl("vir", metacerb$product, ignore.case = T),]
pfam_vir <- pfam_vir[pfam_vir$HMM == "PFAM", "best_hit"]
#Load in KEGG database connections, this was done according to the KO walkthrough in the Anvi'o walkthroughs
kegg_db <- read.delim("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/KEGG/KO_Orthology_R_formated_micro.tsv")
#Trimming whitespace to make sure things merge properly
kegg_db <- apply(kegg_db, 2, trimws) %>% as.data.frame()
#Identify KEGG genes associated with viruses
kegg_vir <- kegg_db[grepl("vir", kegg_db$pathway, ignore.case = T), ]
#Include genes in the Viral proteins metabolism
kegg_vir <- rbind(kegg_vir, kegg_db[kegg_db$metabolism == "Viral proteins", ]) %>% unique()
#Make a list of all of the viral identified genes
vhits <- c(phage_hit, pfam_vir, kegg_vir$ko)
#Assign a flag for if hits are to one of our identified viral genes
metacerb <- metacerb %>% mutate(vgene = case_when(HMM %in% vdbs | best_hit %in% vhits ~ TRUE, .default = FALSE))
#Assign the max and min gene numbers for viral genes in each virus
metacerb <- metacerb %>% group_by(Virus) %>% mutate(vmax = max(gene_num[vgene == TRUE]), vmin =min(gene_num[vgene == TRUE]))
#Keep only non-viral genes in between min and max
metacerb_amgs <- metacerb %>% filter(vmin < gene_num & vmax > gene_num) %>%
  filter(vgene != TRUE)
#Reduce the kegg data set to relevant columns and remove empty hit rows
meta_kegg <- kegg %>% select(target, best_hit, product) %>% filter(best_hit != "")
#Merge with the KEGG database for additional KEGG annotation information
meta_kegg <- merge(meta_kegg, kegg_db, by.x = "best_hit", by.y = "ko", all.x = T)
#Add the virus column
meta_kegg$Virus <- gsub("_[^_]*$", "", meta_kegg$target)
#Make sure whitespaces are removed
meta_kegg <- apply(meta_kegg, 2, trimws) %>% as.data.frame()
#Filter KEGG hits to only potential AMGs in between viral genes
meta_kegg <- meta_kegg[meta_kegg$target %in% metacerb_amgs$target, ]
#write.table(meta_kegg,
#            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/amg_table_stringent.tsv",
#            sep = "\t", quote = FALSE, row.names = FALSE)
setwd("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/")
#Load in RNA data
rna_per <- read.delim("data/all_votu_relative_fractional_transcription.tsv", sep = "\t")
rna_per <- rna_per %>%  select(-c(MAH_C3, POI_C2))
#Load in lifestyle data
life_summ <- read.delim("data/viral_lifestyle_summary.tsv") 
#Add a summary column and keep only relevant information
life_summ <- life_summ %>% mutate(temperate = case_when(n_lysogen_ids > 0 | n_provirus_ids > 0 ~ TRUE, .default = FALSE)) %>% 
  select(Virus, temperate)
#Merge the AMG data with the lifeestyle data and keep relevant information, remove redundant rows
metacerb_life <- merge(meta_kegg, life_summ, by = "Virus") %>% 
  select(Virus, category, pathway, metabolism, gene, temperate) %>% 
  unique()
#Choosing metabolisms of interest for this study
amgs_interest <- read.delim("data/amg_groups_rna.tsv")
full_kegg_filt <- meta_kegg
metabs <- c("Photosynthesis", "Photosynthesis - antenna proteins", "Carotenoid biosynthesis",
           "Porphyrin metabolism",
           "Glycolysis / Gluconeogenesis", "Pentose phosphate pathway",   
           "Sulfur metabolism", "Sulfur relay system")
paths_check <- meta_kegg %>% 
  filter(metabolism %in% metabs) %>%
  select(gene, pathway, metabolism) %>%
  unique()
  
kegg_filt <-  full_kegg_filt %>% select(Virus, gene, pathway, metabolism) %>% unique()
kegg_filt <- merge(kegg_filt, life_summ, by = "Virus") %>%
  select(Virus, pathway, metabolism, gene, temperate)

kegg_per <- merge(kegg_filt, rna_per, by = "Virus") %>% 
  select(-metabolism) %>%
  merge(amgs_interest, by = "gene") %>% 
  select(Virus, label, metabolism, gene, temperate, contains("_C")) 

kegg_per_filt <- kegg_per %>%
  group_by(gene, metabolism, label, temperate) %>%
  summarise_if(is.numeric, sum)

kegg_per_filt_long <- kegg_per_filt %>% pivot_longer(where(is.numeric), names_to = "sample", values_to = "abundance")
kegg_per_filt_long <- kegg_per_filt_long %>% mutate(Lake = case_when(grepl("LB", sample) ~ "Lime Blue",
                                                                     grepl("MAH", sample) ~ "Mahoney",
                                                                     grepl("POI", sample) ~ "Poison"),
                                                    Depth = case_when(grepl("_B", sample) ~ "Bottom",
                                                                      grepl("_S", sample) ~ "Surface",
                                                                      grepl("_C", sample) ~ "Plate")) 

add <- kegg_per_filt_long %>% filter(abundance != 0) %>% pull(abundance) %>% min()
add <- add/10

kegg_per_filt_summary <- kegg_per_filt_long %>%
  group_by(label, metabolism, temperate, Lake, Depth) %>%
  summarise(mean = mean(log10(add + abundance)), sd = std.error(log10(add + abundance))) 

all_amg_df <- kegg_per_filt_summary %>% filter(!is.na(label))
all_amg_df <- all_amg_df %>% arrange(metabolism)
ord <- read.delim("data/amg_orders.tsv")
all_amg_df$label <- factor(all_amg_df$label, levels = ord$label)
all_amg_df <- all_amg_df %>% mutate(mean = case_when(mean == log10(add) ~ NA, .default = mean))
all_amg_df <- all_amg_df %>% mutate(upper = mean + sd,
                                    lower = mean - sd)

amg_dots <- ggplot(all_amg_df, aes(x = label, y = mean, shape = Lake, color = temperate)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) + # Add points for the mean values
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.5, position = position_dodge(width = 0.5))+
  scale_y_continuous(breaks=c( -8, -6, -4, -2, 0), 
                     labels = expression(10^-8, 10^-6, 10^-4, 10^-2, 10),
                     limits = c(-8.5, 0)) +
  scale_color_manual(values = c("black", "slategray")) +
  facet_wrap(~ Lake, ncol = 1)+
  labs(x = NULL,
       y = NULL) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.99, hjust=1))+
  theme(legend.position = "none")

kegg_per_filt_long <- kegg_per_filt_long %>% ungroup()
comp_life <- compare_means(abundance ~ temperate, kegg_per_filt_long, 
                                    group.by = c("label", "Lake"),
                           p.adjust.method = "hochberg")
comp_life_sig <- comp_life %>% filter(p < 0.05)

comp_lake <- compare_means(abundance ~ Lake, kegg_per_filt_long, 
                                    group.by = c("label", "temperate"), 
                                    method = "kruskal.test")
comp_lake_tukey <-  kegg_per_filt_long %>%
  group_by(label, temperate) %>%
  tukey_hsd(abundance ~ Lake) %>%
  add_significance()
comp_lake_sig <- comp_lake_tukey %>% filter(p.adj < 0.05)

# Rank abundance curves ---------------------------------------------------
ra_rna_long <- rna_per %>% 
  pivot_longer(LB_B1:POI_S4, names_to = "sample", values_to = "abundance") %>% 
  mutate(Lake = case_when(grepl("LB", sample) ~ "Lime Blue",
                          grepl("MAH", sample) ~ "Mahoney",
                          grepl("POI", sample) ~ "Poison"),
         Depth = case_when(grepl("_B", sample) ~ "Bottom",
                           grepl("_S", sample) ~ "Surface",
                           grepl("_C", sample) ~ "Plate"))
ra_rna_means <- ra_rna_long %>%
  group_by(Virus, Lake, Depth) %>%
  summarize(mean = mean(abundance))
ra_rna_means_chemo <- ra_rna_means %>% filter(Depth == "Plate")
ra_rna_means_chemo <- ra_rna_means_chemo %>% pivot_wider(names_from = "Lake", values_from = "mean") %>%
  as.data.frame()
ra_rna_means_chemo <- ra_rna_means_chemo %>%
  mutate(LB_rank = row_number(desc(`Lime Blue`)),
         MAH_rank = row_number(desc(`Mahoney`)),
         POI_rank = row_number(desc(`Poison`)))
ra_rna_means_chemo <- merge(ra_rna_means_chemo, life_summ, by = "Virus")

LB_hl <- ra_rna_means_chemo %>%
  filter(`Lime Blue` != 0) %>%
  merge(metacerb_life, by = c("Virus", "temperate"), all.x = TRUE) %>%
  arrange(Virus) %>%
  select(Virus, `Lime Blue`, LB_rank, metabolism, temperate) %>% 
  unique()
LB_hl_light <- LB_hl %>% filter(metabolism %in% c("Photosynthesis", 
                                                      "Porphyrin metabolism", "Other carbon fixation pathways", 
                                                      "Photosynthesis - antenna proteins"))
LB_hl_sulf <- LB_hl %>% filter(metabolism %in% c("Sulfur metabolism", 
                                                 "Glycolysis / Gluconeogenesis",
                                                 "Sulfur relay system"))

LB_ra <- ra_rna_means_chemo %>%
  filter(`Lime Blue` != 0) %>%
  ggplot(aes(x = LB_rank, y = `Lime Blue`)) +
  geom_segment(data = LB_hl_light,
               aes(x = LB_rank, y = 0.000000001, xend = LB_rank, yend = `Lime Blue`, color = metabolism, linetype = temperate),
               linewidth = 0.6) +
  geom_segment(data = LB_hl_sulf,
               aes(x = LB_rank, y = 1, xend = LB_rank, yend = `Lime Blue`, color = metabolism, linetype = temperate),
               linewidth = 0.6) +
  scale_color_manual(values = pal)+
  scale_linetype_manual(values = c("solid", "dotted"))+
  new_scale_color()+
  geom_point()+
  geom_point(aes(alpha = temperate), color = "gray") +
  scale_alpha_manual(values = c(0, 1))+
  scale_y_log10(limits = c(1, 0.000000001), breaks = c(1, 0.01, 0.0001, 0.000001, 0.00000001),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + # Log-transform y-axis
  labs(title = "Lime Blue",
       x = NULL,
       y = NULL) +
  annotation_logticks(sides = "l")+
  theme_pubr()

MAH_hl <- ra_rna_means_chemo %>%
  filter(`Mahoney` != 0) %>%
  merge(metacerb_life, by = c("Virus", "temperate"), all.x = TRUE) %>%
  arrange(Virus) %>%
  select(Virus, `Mahoney`, MAH_rank, metabolism, temperate) %>% 
  unique()
MAH_hl_light <- MAH_hl %>% filter(metabolism %in% c("Photosynthesis", 
                                                    "Porphyrin metabolism", "Other carbon fixation pathways", 
                                                    "Photosynthesis - antenna proteins"))
MAH_hl_sulf <- MAH_hl %>% filter(metabolism %in% c("Sulfur metabolism", "Pentose phosphate pathway", 
                                                   "Glycolysis / Gluconeogenesis",
                                                   "Sulfur relay system"))

MAH_ra <- ra_rna_means_chemo %>%
  filter(`Mahoney` != 0) %>%
  ggplot(aes(x = MAH_rank, y = `Mahoney`)) +
  geom_segment(data = MAH_hl_light,
               aes(x = MAH_rank, y = 0.000000001, xend = MAH_rank, yend = `Mahoney`, color = metabolism, linetype = temperate),
               linewidth = 0.6) +
  geom_segment(data = MAH_hl_sulf,
               aes(x = MAH_rank, y = 1, xend = MAH_rank, yend = `Mahoney`, color = metabolism, linetype = temperate),
               linewidth = 0.6) +
  scale_color_manual(values = pal)+
  scale_linetype_manual(values = c("solid", "dotted"))+
  new_scale_color()+
  geom_point()+
  geom_point(aes(alpha = temperate), color = "gray") +
  scale_alpha_manual(values = c(0, 1))+
  scale_y_log10(limits = c(1, 0.000000001), breaks = c(1, 0.01, 0.0001, 0.000001, 0.00000001),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + # Log-transform y-axis
  labs(title = "Mahoney",
       x = NULL,
       y = NULL) +
  annotation_logticks(sides = "l")+
  theme_pubr()

POI_hl <- ra_rna_means_chemo %>%
  filter(`Poison` != 0) %>%
  merge(metacerb_life, by = c("Virus", "temperate"), all.x = TRUE) %>%
  arrange(Virus) %>%
  select(Virus, `Poison`, POI_rank, metabolism, temperate) %>% 
  unique()
POI_hl_light <- POI_hl%>% filter(metabolism %in% c("Photosynthesis", 
                                                   "Porphyrin metabolism", "Other carbon fixation pathways", 
                                                   "Photosynthesis - antenna proteins"))
POI_hl_sulf <- POI_hl %>% filter(metabolism %in% c("Sulfur metabolism", 
                                                   "Glycolysis / Gluconeogenesis",
                                                   "Sulfur relay system"))

POI_ra <- ra_rna_means_chemo %>%
  filter(`Poison` != 0) %>%
  ggplot(aes(x = POI_rank, y = `Poison`)) +
  geom_segment(data = POI_hl_light,
               aes(x = POI_rank, y = 0.000000001, xend = POI_rank, yend = `Poison`, color = metabolism, linetype = temperate),
               linewidth = 0.6) +
  geom_segment(data = POI_hl_sulf,
               aes(x = POI_rank, y = 1, xend = POI_rank, yend = `Poison`, color = metabolism, linetype = temperate),
               linewidth = 0.6) +
  scale_color_manual(values = pal)+
  scale_linetype_manual(values = c("solid", "dotted"))+
  new_scale_color()+
  geom_point()+
  geom_point(aes(alpha = temperate), color = "gray") +
  scale_alpha_manual(values = c(0, 1))+
  scale_y_log10(limits = c(1, 0.000000001), breaks = c(1, 0.01, 0.0001, 0.000001, 0.00000001),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + # Log-transform y-axis
  labs(title = "Poison",
       x = NULL,
       y = NULL) +
  annotation_logticks(sides = "l")+
  theme_pubr()

ra_curves <- ggarrange(LB_ra, MAH_ra, POI_ra, nrow = 1, common.legend = T, legend = "none")

# AMG counts --------------------------------------------------------------
plate_viruses <- ra_rna_means_chemo %>% filter((`Lime Blue` + Mahoney + Poison) != 0)

n_amgs <- metacerb_life %>% 
  group_by(Virus, metabolism, temperate) %>%
  filter(Virus %in% plate_viruses$Virus) %>%
  unique() %>%
  group_by(metabolism, temperate) %>%
  summarise(n = n())

data_pie <- n_amgs %>%
  group_by(metabolism) %>%
  arrange(metabolism, desc(temperate)) %>%
  mutate(percentage = round(n / sum(n), digits = 2),
         ypos = cumsum(percentage) - percentage / 2)
data_pie$metabolism <- factor(data_pie$metabolism, 
                              levels = c("Photosynthesis", "Photosynthesis - antenna proteins", "Other carbon fixation pathways",
                                         "Porphyrin metabolism",
                                         "Glycolysis / Gluconeogenesis", "Pentose phosphate pathway",   
                                         "Sulfur metabolism", "Sulfur relay system"))
data_pie <- data_pie %>% filter(!is.na(metabolism))
# Create the faceted pie chart with text
pies <- ggplot(data_pie, aes(x = "", y = percentage, fill = temperate)) +
  geom_bar(width = 1, stat = "identity") +
  geom_text(aes(y = ypos, label = n, color = temperate), size = 6) + # Add labels
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("black", "slategray"))+
  scale_color_manual(values = c("white", "black"))+
  facet_wrap(~ metabolism, nrow = 1) + # Facet by the 'group' variable
  theme_void() + # Remove axis and grid lines
  theme(legend.position = "none")

ggarrange(ra_curves, amg_dots, pies, ncol = 1, heights = c(1.5,4,1))

ra_leg <- get_legend(LB_ra) %>% as_ggplot()

amg_leg <- ggplot(all_amg_df, aes(x = label, y = mean, shape = Lake, color = temperate)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) + # Add points for the mean values
  scale_color_manual(values = c("black", "slategray"))
amg_leg <- get_legend(amg_leg) %>% as_ggplot()
ggarrange(ra_leg, amg_leg)

hosts <- read.delim("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/qc_votus_host_matches.tsv")

host_amgs <- meta_kegg %>% 
  filter(Virus %in% plate_viruses$Virus) %>%
  merge(hosts, by.x = "Virus", by.y = "virus")
host_amgs <- host_amgs %>% filter(metabolism %in% metabs)

#Temperate Wilcoxon Test
test_kegg <- kegg_per %>%
  group_by(gene, metabolism) %>%
  summarise_if(is.numeric, sum)

test_kegg_long <- test_kegg %>% pivot_longer(where(is.numeric), 
                                             names_to = "sample", values_to = "abundance") %>% 
  ungroup()
test_kegg_long <- test_kegg_long %>% mutate(Lake = case_when(grepl("LB", sample) ~ "Lime Blue",
                                                                     grepl("MAH", sample) ~ "Mahoney",
                                                                     grepl("POI", sample) ~ "Poison"),
                                                    Depth = case_when(grepl("_B", sample) ~ "Bottom",
                                                                      grepl("_S", sample) ~ "Surface",
                                                                      grepl("_C", sample) ~ "Plate"))
test_comp_summ <- test_kegg_long %>% group_by(metabolism, Lake) %>% summarize(mean = mean(abundance))
test_comp_lake <- compare_means(abundance ~ Lake, test_kegg_long, 
                           group.by = c("metabolism"), 
                           method = "kruskal.test")
test_comp_lake_tukey <-  test_kegg_long %>%
  group_by(metabolism) %>%
  tukey_hsd(abundance ~ Lake) %>%
  add_significance()
test_comp_lake_sig <- test_comp_lake_tukey %>% filter(p.adj < 0.05)
