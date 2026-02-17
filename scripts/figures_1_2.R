library(dplyr)
library(tidyr)
library(rstatix)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggpattern)
library(plotrix)

# MAG Information ---------------------------------------------------------
setwd("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/")
mag_func <- read.table("annotation_summary_bins.tsv", sep = "\t", quote = "",
                       col.names = c("locus_tag",	"FOAM",	"KEGG",	"COG",	"CAZy",	"PHROG",
                                     "VOG",	"Best hit",	"length_bp",	"e-value",	"score",	
                                     "EC_number",	"product"))
mag_func$contig <- gsub("_[^_]*$", "", mag_func$locus_tag)
taxa <- read.delim("mag_taxonomy.tsv")
c2b <- read.delim("contigs2bins.tsv", header = F, col.names = c("bin", "contig"))

anoxi <- mag_func %>% filter(grepl("K0894[0-3]", KEGG))
anoxi_bin <- merge(c2b, anoxi, by = "contig")
anoxi_bin <- anoxi_bin %>% filter(bin %in% anoxi_bin[duplicated(anoxi_bin$bin), "bin"])
anoxii <- mag_func %>% filter(grepl("K0892[89]", KEGG))
anoxii_bin <- merge(c2b, anoxii, by = "contig")
anox_bin <- rbind(anoxi_bin, anoxii_bin)
anox_taxa <- merge(anox_bin, taxa, by.x = "bin", by.y = "user_genome")

srb <- mag_func %>% 
  filter(KEGG %in% c("K00958", "K00394", "K00395", "K11180", "K11181", "K27196", "K27187", "K27188", "K27189", "K27190", "K27191"))
srb_bin <- merge(c2b, srb, by = "contig")
srb_taxa <- merge(srb_bin, taxa, by.x = "bin", by.y = "user_genome")
srb_taxa <- srb_taxa %>% select(bin, domain:genus) %>%
  unique()

srb_taxa_1 <- taxa %>% filter(order %in% c("o__Desulfobacterales", "o__Desulfobulbales", "o__Desulfatiglandales"))

rna_mag_per <- read.delim("all_mag_relative_transcription.tsv")
rna_mag_per <- rna_mag_per %>% select(-c(MAH_C3, POI_C2))
dna_mag_per <- read.delim("all_mag_relative_fractional_abundance.tsv")

relative_taxa <- anox_taxa %>% select(bin, domain:genus) %>%
  unique() %>%
  merge(dna_mag_per, by = "bin") %>%
  merge(rna_mag_per, by = "bin", suffixes = c("_DNA", "_RNA"))
grp_rel_taxa <- relative_taxa %>%
  mutate(anox_grp = case_match(order,
                               "o__Chlorobiales" ~ "GSB",
                               "o__Chromatiales" ~ "PSB",
                               .default = "OAP")) %>%
  group_by(anox_grp) %>%
  summarise_if(is.numeric, sum)
long <- grp_rel_taxa %>% pivot_longer(LB_B1_DNA:POI_S4_RNA, names_to = "samples", 
                                      values_to = "Amount")
long <- long %>% mutate(Lake = case_when(grepl("LB", samples) ~ "Lime Blue",
                                         grepl("MAH", samples) ~ "Mahoney",
                                         grepl("POI", samples) ~ "Poison"),
                        Depth = case_when(grepl("_B", samples) ~ "Bottom",
                                          grepl("_C", samples) ~ "Plate",
                                          grepl("_S", samples) ~ "Surface"),
                        Replicate = case_when(grepl("1", samples) ~ "1",
                                              grepl("2", samples) ~ "2",
                                              grepl("3", samples) ~ "3",
                                              grepl("4", samples) ~ "4",
                                              grepl("5", samples) ~ "5"),
                        LakeDepth = paste0(Lake, Depth),
                        Type = case_when(grepl("RNA", samples) ~ "RNA", .default = "DNA"))
summ_long <- long %>% group_by(anox_grp, Lake, Depth, Type) %>% 
  summarise(mean = mean(Amount), sd = std.error(Amount))
summ_long$Depth <- factor(summ_long$Depth, levels = c("Surface", "Plate", "Bottom"))
summ_long$anox_grp <- factor(summ_long$anox_grp, levels = c("OAP", "GSB", "PSB"))

mag_tukey <- long %>% 
  filter(LakeDepth != "MahoneyBottom") %>% 
  group_by(Lake, Depth, anox_grp) %>%
  mutate(anox_grp = factor(anox_grp, levels = c("OAP", "GSB", "PSB")),
         Depth = factor(Depth, levels = c("Surface", "Plate", "Bottom"))) %>%
  tukey_hsd(Amount ~ Type) %>%
  add_significance()

mag_tukey <- mag_tukey %>% add_xy_position(x = "anox_grp", scales = "fixed", step.increase = 0.05)
mag_tukey <- mag_tukey %>% group_by(group1, group2) %>% 
  mutate(y.position = max(y.position)) %>% 
  as.data.frame()
mag_tukey$Depth <- factor(mag_tukey$Depth, levels = c("Surface", "Plate", "Bottom"))
mag_tukey$anox_grp <- factor(mag_tukey$anox_grp, levels = c("OAP", "GSB", "PSB"))

mag_tukey_plate <- long %>% 
  filter(Depth == "Plate") %>% 
  group_by(Lake, Depth, anox_grp) %>%
  mutate(anox_grp = factor(anox_grp, levels = c("OAP", "GSB", "PSB"))) %>%
  tukey_hsd(Amount ~ Type) %>%
  add_significance()

mag_tukey_plate <- mag_tukey_plate %>% 
  add_xy_position(x = "anox_grp", scales = "fixed")
mag_tukey_plate <- mag_tukey_plate %>% 
  group_by(group1, group2) %>% 
  mutate(y.position = c(rep(0.4, 3), rep(0.85, 6))) %>% 
  as.data.frame()

micro_full <- ggplot(data = summ_long, aes(x = anox_grp, y = mean, fill = anox_grp))+
  geom_col_pattern(aes(pattern = Type), 
                   position = position_dodge(width = 0.9), 
                   colour = "black", pattern_fill = "black") +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, group = Type), 
                position = position_dodge(width = 0.9),
                width = 0.25) +
  stat_pvalue_manual(mag_tukey,  coord.flip = TRUE, hide.ns = FALSE, tip.length = 0.01)+
  coord_flip()+
  facet_grid(Depth ~ Lake) +
  scale_pattern_manual(values = c("none", "stripe")) + 
  scale_fill_manual(name = "Anox Groups", values  = c("slategray", "green4", "orchid3")) +
  labs(y = "Mean Relative Abundance/Activity", title = "Anoxic Phototroph Abundance/Activity") +
  theme_pubr()+ 
  theme(axis.title.y = element_blank(), legend.position = "bottom")
  
micro_pan <- summ_long %>%
  filter(Depth == "Plate") %>%
  ggplot(aes(x = anox_grp, y = mean, fill = anox_grp))+
  geom_col_pattern(aes(pattern = Type), 
                   position = position_dodge(width = 0.9), 
                   colour = "black", , pattern_fill = "black") +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, group = Type), 
                position = position_dodge(width = 0.9),
                width = 0.25, linewidth = 1) +
  stat_pvalue_manual(mag_tukey_plate,
                     coord.flip = TRUE, hide.ns = FALSE, tip.length = 0.01)+
  coord_flip()+
  facet_wrap(~ Lake, ncol = 1, strip.position = "left", scales = "free") +
  scale_pattern_manual(values = c("none", "circle")) + 
  scale_fill_manual(values  = c("slategray", "green4", "orchid3")) +
  labs(y = "Mean Relative Abundance/Activity",
       title = "Anoxygnic Phototrophs")+
  theme_pubr() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

# SRB Plot ----------------------------------------------------------------
srb_taxa <- taxa %>% 
  filter(order %in% c("o__Desulfobacterales", "o__Desulfobulbales", "o__Desulfatiglandales")) %>%
  rename(c(bin = "user_genome"))

srb_rel_taxa <- srb_taxa %>% select(bin, domain:genus) %>%
  unique() %>%
  merge(dna_mag_per, by = "bin") %>%
  merge(rna_mag_per, by = "bin", suffixes = c("_DNA", "_RNA"))
srb_sum_rel_taxa <- srb_rel_taxa %>%
  summarise_if(is.numeric, sum)
srb_long <- srb_sum_rel_taxa %>% pivot_longer(LB_B1_DNA:POI_S4_RNA, names_to = "samples", 
                                      values_to = "Amount")
srb_long <- srb_long %>% mutate(Lake = case_when(grepl("LB", samples) ~ "Lime Blue",
                                         grepl("MAH", samples) ~ "Mahoney",
                                         grepl("POI", samples) ~ "Poison"),
                        Depth = case_when(grepl("_B", samples) ~ "Bottom",
                                          grepl("_C", samples) ~ "Plate",
                                          grepl("_S", samples) ~ "Surface"),
                        Replicate = case_when(grepl("1", samples) ~ "1",
                                              grepl("2", samples) ~ "2",
                                              grepl("3", samples) ~ "3",
                                              grepl("4", samples) ~ "4",
                                              grepl("5", samples) ~ "5"),
                        LakeDepth = paste0(Lake, Depth),
                        Type = case_when(grepl("RNA", samples) ~ "RNA", .default = "DNA"))
srb_summ_long <- srb_long %>% group_by(Lake, Depth, Type) %>% 
  summarise(mean = mean(Amount), sd = std.error(Amount))
srb_summ_long$Depth <- factor(srb_summ_long$Depth, levels = c("Bottom", "Plate", "Surface"))

srb_tukey <- srb_long %>% filter(LakeDepth != "MahoneyBottom") %>% 
  group_by(Lake, Depth) %>%
  tukey_hsd(Amount ~ Type) %>%
  add_significance()
srb_tukey <- srb_tukey %>% 
  add_xy_position(x = "Depth", scales = "fixed", step.increase = 0.05)
srb_tukey <- srb_tukey %>% group_by(group1, group2) %>% 
  mutate(y.position = max(y.position)) %>% 
  as.data.frame()

srb_tukey_plate <- srb_long %>% 
  filter(Depth == "Plate") %>% 
  group_by(Lake, Depth) %>%
  tukey_hsd(Amount ~ Type) %>%
  add_significance()
srb_tukey_plate <- srb_tukey_plate %>% 
  add_y_position(scales = "fixed")
srb_tukey_plate <- srb_tukey_plate %>% 
  group_by(group1, group2) %>% 
  mutate(y.position = max(y.position)) %>% 
  as.data.frame()

srb_full <- ggplot(data = srb_summ_long, 
                   aes(x = Depth, y = mean, fill = Depth))+
  geom_col_pattern(aes(pattern = Type), 
                   position = position_dodge(width = 0.9), 
                   colour = "black", pattern_fill = "black") +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, group = Type), 
                position = position_dodge(width = 0.9),
                width = 0.25) +
  stat_pvalue_manual(srb_tukey,  coord.flip = TRUE, 
                     hide.ns = FALSE, tip.length = 0.01)+
  coord_flip()+
  facet_wrap(~ Lake, ncol = 1) +
  scale_pattern_manual(values = c("none", "stripe")) + 
  scale_fill_manual(values = c("darkgoldenrod", "violetred", "paleturquoise3"))+
  labs(y = "Mean Relative Abundance/Activity",
       x = "SRB MAGs",
       title = "Sulfate Reducer Abundance/Activity") +
  theme_pubr() +
  theme(legend.position = "bottom")

srb_pan <- srb_summ_long %>%
  filter(Depth == "Plate") %>%
  ggplot(aes(x = Type, y = mean, fill = Type))+
  geom_col_pattern(aes(pattern = Type), 
                   position = position_dodge(width = 0.9), 
                   colour = "black", , pattern_fill = "black") +  
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, group = Type), 
                position = position_dodge(width = 0.9),
                width = 0.25, linewidth = 1) +
  stat_pvalue_manual(srb_tukey_plate,
                     coord.flip = TRUE, hide.ns = FALSE, tip.length = 0.01)+
  coord_flip()+
  facet_wrap(~ Lake, ncol = 1) +
  scale_pattern_manual(values = c("none", "circle")) + 
  scale_fill_manual(values = c("gray", "white"))+
  labs(y = "Mean Relative Abundance/Activity",
       title = "SRB MAGs")+
  theme_pubr() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
# Read in information on Viruses and other data ------------------------------------------
setwd("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/")
life_summ <- read.delim("data/viral_lifestyle_summary.tsv") #Read in lifestyle summary
qual <- read.delim("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/input_data/gv_votus_checkv/quality_summary.tsv") #Read in checkv quality
dna_per <- read.delim("data/all_votu_relative_fractional_abundance.tsv", sep = "\t")
meta <- read.delim("data/space_2023_metadata.txt")
meta <- meta[!is.na(meta$DO), ]
samp_dep <- read.delim("data/space_2023_sampling_depth.txt")
rna_per <- read.delim("data/all_votu_relative_fractional_transcription.tsv", sep = "\t")
rna_per <- rna_per %>% select(-c(`MAH_C3`,`POI_C2`))

# Make summary tables and filter viruses ----------------------------------
life_summ <- life_summ[life_summ$Virus %in% dna_per$Virus, ] #Strictest filter 0 host genes

life_summary <- function(x) {
  life_table <- x[, 1:3] #Reduce dataset
  life_table <- life_table %>% mutate(temp_virus = case_when(n_lysogen_ids >= 1 & n_provirus_ids >= 1 ~ "Pro & Temp",
                                                             n_lysogen_ids >= 0 & n_provirus_ids >= 1 ~ "Prophage Only",
                                                             n_lysogen_ids >= 1 & n_provirus_ids >= 0 ~ "Temperate Only",
                                                             .default = "Lytic")) #Assign temperate label on criteria
  table(life_table$temp_virus) #View counts of the different types of temperate virus identifications
}

life_summary(life_summ)

# Make metadata plots -----------------------------------------------------
meta <- merge(meta, samp_dep, by = "Lake") %>% 
  pivot_wider(names_from = "sampling_type", values_from = "sampling_depth_m")

do_depth <- ggplot(meta, aes(x = Depth, y = DO))+
  geom_line(aes(color = "DO (%)"), linewidth = 1)+
  geom_line(data = meta[!is.na(meta$PAR), ], aes(y = (PAR/4.5), color = "PAR"), linewidth = 1)+
  scale_y_continuous(name = "DO (%)", sec.axis = sec_axis(~ .*4.5, name = "PAR"))+
  facet_wrap(~ Lake, ncol = 1, scales = "free_y", strip.position = "left")+
  coord_flip()+
  scale_x_reverse(name = "Depth (m)", position = "top")+
  geom_vline(aes(xintercept = Surface), colour ="paleturquoise3", lwd = 0.75, linetype = 2)+
  geom_vline(aes(xintercept = Chemocline), colour ="violetred", lwd = 0.75, linetype = 2)+
  geom_vline(aes(xintercept = Bottom), colour ="darkgoldenrod", lwd = 0.75, linetype = 2)+
  scale_color_manual(values = c("black", "yellowgreen"),
                     breaks = c("DO (%)", "PAR"),
                     name = "Measurements")+
  theme_pubr()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.title.y = element_text(angle = 90),
        legend.position = "none",
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0.5,0,2,0))


# Fractional abundance lifestyles-------------------------------------------
sum_all <- life_summ %>% 
  mutate(Lysogenic = case_when(n_lysogen_ids > 0 ~ "Temperate",
                               .default = "Lytic",),
         Provirus = case_when(n_provirus_ids > 0 ~ "Prophage",
                              .default = "Lytic")) %>%
  select(Virus, Lysogenic, Provirus)

sum_pro <- sum_all %>% mutate(Lifestyle = case_when(Lysogenic == "Temperate" ~ "Temperate", 
                                        Provirus == "Prophage" ~ "Temperate", 
                                        Lysogenic == "Lytic" & Provirus == "Lytic" ~ "Lytic"),
                              .keep = "unused")

life_dnas <- merge(sum_pro, dna_per, by = "Virus")
lifestyle_dna <- life_dnas %>% 
  group_by(Lifestyle) %>% 
  summarise_if(is.numeric, sum) %>%
  t() 
colnames(lifestyle_dna) <-lifestyle_dna[1, ]
lifestyle_dna <- lifestyle_dna[-1, ]
lifestyle_dna <- apply(lifestyle_dna, 2, as.numeric) %>%
  as.data.frame()
lifestyle_dna$sample <- colnames(life_dnas)[3:length(life_dnas)]
lifestyle_dna <- lifestyle_dna %>% mutate(Lake = case_when(grepl("LB", sample) ~ "Lime Blue",
                                      grepl("MAH", sample) ~ "Mahoney",
                                      grepl("POI", sample) ~ "Poison"),
                     Depth = case_when(grepl("_B", sample) ~ "Bottom",
                                      grepl("_S", sample) ~ "Surface",
                                      grepl("_C", sample) ~ "Plate"),
                     LakeDepth = paste0(Lake, Depth))
lifestyle_dna$Depth <- factor(lifestyle_dna$Depth, 
                              levels = c("Bottom", "Plate", "Surface"))
lifestyle_dna$Type <- "DNA"

# RPKM Figures ------------------------------------------------------------
life_rnas <- merge(sum_pro, rna_per, by = "Virus")
lifestyle_rna <- life_rnas %>% 
  group_by(Lifestyle) %>% 
  summarise_if(is.numeric, sum) %>%
  t() 
colnames(lifestyle_rna) <- lifestyle_rna[1, ]
lifestyle_rna <- lifestyle_rna[-1, ]
lifestyle_rna <- apply(lifestyle_rna, 2, as.numeric) %>%
  as.data.frame()
lifestyle_rna$sample <- colnames(life_rnas)[3:length(life_rnas)]

lifestyle_rna <- lifestyle_rna %>% mutate(Lake = case_when(grepl("LB", sample) ~ "Lime Blue",
                                                           grepl("MAH", sample) ~ "Mahoney",
                                                           grepl("POI", sample) ~ "Poison"),
                                          Depth = case_when(grepl("_B", sample) ~ "Bottom",
                                                            grepl("_S", sample) ~ "Surface",
                                                            grepl("_C", sample) ~ "Plate"),
                                          LakeDepth = paste0(Lake, Depth))
lifestyle_rna$Depth <- factor(lifestyle_rna$Depth, 
                              levels = c("Bottom", "Plate", "Surface"))
lifestyle_rna$Type <- "RNA"

lifestyle_full <- rbind(lifestyle_dna, lifestyle_rna)

life_summ <- lifestyle_full %>%
  group_by(Lake, Depth, Type) %>%
  summarise(temp_mean = mean(Temperate), temp_sd = std.error(Temperate),
            lyt_mean = mean(Lytic), lyt_sd = std.error(Lytic))

life_summ$Depth <- factor(life_summ$Depth, levels = c("Bottom", "Plate", "Surface"))

temp_tukey <- lifestyle_full %>% filter(LakeDepth != "MahoneyBottom") %>% group_by(Lake, Depth) %>%
  tukey_hsd(Temperate ~ Type) %>%
  add_significance()
temp_tukey <- temp_tukey %>% 
  add_xy_position(x = "Depth", scales = "fixed", step.increase = 0.05)
temp_tukey <- temp_tukey %>% group_by(group1, group2) %>% 
  mutate(y.position = max(y.position)) %>% 
  as.data.frame()

temp_tukey_plate <- lifestyle_full %>% 
  filter(Depth == "Plate") %>% 
  group_by(Lake, Depth) %>%
  tukey_hsd(Temperate ~ Type) %>%
  add_significance()

temp_tukey_plate <- temp_tukey_plate %>% 
  add_y_position(scales = "fixed")
temp_tukey_plate <- temp_tukey_plate %>% 
  group_by(group1, group2) %>% 
  mutate(y.position = max(y.position)) %>% 
  as.data.frame()

lyt_tukey <- lifestyle_full %>% filter(LakeDepth != "MahoneyBottom") %>% group_by(Lake, Depth) %>%
  tukey_hsd(Lytic ~ Type) %>%
  add_significance()
lyt_tukey <- lyt_tukey %>% 
  add_xy_position(x = "Depth", scales = "fixed", step.increase = 0.05)
lyt_tukey <- lyt_tukey %>% group_by(group1, group2) %>% 
  mutate(y.position = max(y.position)) %>% 
  as.data.frame()

lyt_tukey_plate <- lifestyle_full %>% 
  filter(Depth == "Plate") %>% 
  group_by(Lake, Depth) %>%
  tukey_hsd(Lytic ~ Type) %>%
  add_significance()

lyt_tukey_plate <- lyt_tukey_plate %>% 
  add_y_position(scales = "fixed")
lyt_tukey_plate <- lyt_tukey_plate %>% 
  group_by(group1, group2) %>% 
  mutate(y.position = max(y.position)) %>% 
  as.data.frame()

temp_full <- ggplot(data = life_summ, 
                    aes(x = Depth, y = temp_mean, fill = Depth))+
  geom_col_pattern(aes(pattern = Type), 
                   position = position_dodge(width = 0.9), 
                   colour = "black", pattern_fill = "black") +  
  geom_errorbar(aes(ymin = temp_mean - temp_sd, ymax = temp_mean + temp_sd, 
                    group = Type), 
                position = position_dodge(width = 0.9),
                width = 0.25) +
  stat_pvalue_manual(temp_tukey,  coord.flip = TRUE, 
                     hide.ns = FALSE, tip.length = 0.01)+
  coord_flip()+
  facet_wrap(~ Lake, ncol = 1) +
  scale_pattern_manual(values = c("none", "stripe")) + 
  scale_fill_manual(values = c("darkgoldenrod", "violetred", "paleturquoise3"))+
  labs(y = "Mean Relative Abundance/Activity",
       x = "Temperate vOTUs") +
  theme_pubr() +
  theme(axis.text.y = element_blank())

temp_pan <- life_summ %>%
  filter(Depth == "Plate") %>%
  ggplot(aes(x = Type, y = temp_mean, fill = Type))+
  geom_col_pattern(aes(pattern = Type), 
                   position = position_dodge(width = 0.9), 
                   colour = "black", , pattern_fill = "black") +  
  geom_errorbar(aes(ymin = temp_mean - temp_sd, ymax = temp_mean + temp_sd, 
                    group = Type), 
                position = position_dodge(width = 0.9),
                width = 0.25, linewidth = 1) +
  stat_pvalue_manual(temp_tukey_plate,
                     coord.flip = TRUE, hide.ns = FALSE, tip.length = 0.01)+
  coord_flip()+
  facet_wrap(~ Lake, ncol = 1) +
  scale_pattern_manual(values = c("none", "circle")) + 
  scale_fill_manual(values = c("gray", "white"))+
  labs(y = "Mean Relative Abundance/Activity",
       title = "Temperate vOTUs")+
  theme_pubr() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

lyt_full <- ggplot(data = life_summ, 
                    aes(x = Depth, y = lyt_mean, fill = Depth))+
  geom_col_pattern(aes(pattern = Type), 
                   position = position_dodge(width = 0.9), 
                   colour = "black", pattern_fill = "black") +  
  geom_errorbar(aes(ymin = lyt_mean - lyt_sd, ymax = lyt_mean + lyt_sd, 
                    group = Type), 
                position = position_dodge(width = 0.9),
                width = 0.25) +
  stat_pvalue_manual(lyt_tukey,  coord.flip = TRUE, 
                     hide.ns = FALSE, tip.length = 0.01)+
  coord_flip()+
  facet_wrap(~ Lake, ncol = 1) +
  scale_pattern_manual(values = c("none", "stripe")) + 
  scale_fill_manual(values = c("darkgoldenrod", "violetred", "paleturquoise3"))+
  labs(y = "Mean Relative Abundance/Activity",
       x = "Lytic vOTUs") +
  theme_pubr() +
  theme(axis.text.y = element_blank())

lyt_pan <- life_summ %>%
  filter(Depth == "Plate") %>%
  ggplot(aes(x = Type, y = lyt_mean, fill = Type))+
  geom_col_pattern(aes(pattern = Type), 
                   position = position_dodge(width = 0.9), 
                   colour = "black", , pattern_fill = "black") +  
  geom_errorbar(aes(ymin = lyt_mean - lyt_sd, ymax = lyt_mean + lyt_sd, 
                    group = Type), 
                position = position_dodge(width = 0.9),
                width = 0.25, linewidth = 1) +
  stat_pvalue_manual(lyt_tukey_plate,
                     coord.flip = TRUE, hide.ns = FALSE, tip.length = 0.01)+
  coord_flip()+
  facet_wrap(~ Lake, ncol = 1) +
  scale_pattern_manual(values = c("none", "circle")) + 
  scale_fill_manual(values = c("gray", "white"))+
  labs(y = "Mean Relative Abundance/Activity",
       title = "Lytic vOTUs")+
  theme_pubr() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
# Host Defenses -----------------------------------------------------------
defense <- read.delim("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/input_data/hosts/all_postderep_bins_defense_finder_systems.tsv")
defense <- defense %>% mutate(contig = gsub("_[^_]*$", "", sys_beg)) #%>% select(contig, type)
defense <- merge(c2b, defense, by = "contig") %>% select(-contig) %>% unique()
defense <- defense %>% filter(type == "RM" | type == "Cas")
def_red <- defense %>% select(bin, type) %>% rename(c(System = "type"))

def_abun_act <- def_red %>%
  unique() %>%
  merge(dna_mag_per, by = "bin") %>%
  merge(rna_mag_per, by = "bin", suffixes = c("_DNA", "_RNA"))
def_long <- def_abun_act %>% 
  pivot_longer(where(is.numeric), names_to = "sample", values_to = "Amount") %>% 
  mutate(Lake = case_when(grepl("LB", sample) ~ "Lime Blue",
                          grepl("MAH", sample) ~ "Mahoney",
                          grepl("POI", sample) ~ "Poison"),
         Depth = case_when(grepl("_B", sample) ~ "Bottom",
                           grepl("_S", sample) ~ "Surface",
                           grepl("_C", sample) ~ "Plate"),
         LakeDepth = paste0(Lake, Depth),
         Type = case_when(grepl("RNA", sample) ~ "RNA", .default = "DNA"))

def_summ <- def_long %>%
  group_by(Lake, Depth, Type, System, sample) %>%
  summarise(total = sum(Amount)) %>%
  ungroup() %>%
  pivot_wider(names_from = "System", values_from = "total") %>%
  group_by(Lake, Depth, Type) %>%
  summarise(rm_mean = mean(RM), rm_sd = std.error(RM),
            cas_mean = mean(Cas), cas_sd = std.error(Cas))

def_summ$Depth <- factor(def_summ$Depth, levels = c("Bottom", "Plate", "Surface"))

rm_tukey <- def_long %>% 
  filter(System == "RM") %>%
  filter(LakeDepth != "MahoneyBottom") %>% 
  group_by(Lake, Depth, Type, sample) %>%
  summarise(total = sum(Amount)) %>%
  ungroup() %>%
  group_by(Lake, Depth) %>%
  mutate(Depth = factor(Depth, levels = c("Surface", "Plate", "Bottom"))) %>%
  tukey_hsd(total ~ Type) %>%
  add_significance()

rm_tukey <- rm_tukey %>% add_xy_position(x = "Depth", scales = "fixed", step.increase = 0.05)
rm_tukey <- rm_tukey %>% group_by(group1, group2) %>% 
  mutate(y.position = 0.99) %>% 
  as.data.frame()
rm_tukey$Depth <- factor(rm_tukey$Depth, levels = c("Surface", "Plate", "Bottom"))

rm_tukey_plate <- def_long %>% 
  filter(System == "RM") %>%
  filter(Depth == "Plate") %>% 
  group_by(Lake, Depth, Type, sample) %>%
  summarise(total = sum(Amount)) %>%
  ungroup() %>%
  group_by(Lake, Depth) %>%
  tukey_hsd(total ~ Type) %>%
  add_significance()

rm_tukey_plate <- rm_tukey_plate %>% 
  add_xy_position(x = "Type", scales = "fixed")
rm_tukey_plate <- rm_tukey_plate %>% 
  group_by(group1, group2) %>% 
  mutate(y.position = max(y.position)) %>% 
  as.data.frame()

rm_full <- ggplot(data = def_summ, aes(x = Depth, y = rm_mean, fill = Depth))+
  geom_col_pattern(aes(pattern = Type), 
                   position = position_dodge(width = 0.9), 
                   colour = "black", pattern_fill = "black") +  
  geom_errorbar(aes(ymin = rm_mean - rm_sd, ymax = rm_mean + rm_sd, group = Type), 
                position = position_dodge(width = 0.9),
                width = 0.25) +
  stat_pvalue_manual(rm_tukey,  coord.flip = TRUE, hide.ns = FALSE, tip.length = 0.01)+
  coord_flip()+
  facet_wrap(~ Lake, ncol = 1) +
  scale_pattern_manual(values = c("none", "stripe")) + 
  scale_fill_manual(values = c("darkgoldenrod", "violetred", "paleturquoise3"))+
  labs(y = "Mean Relative Abundance/Activity",
       x = "MAGs with RM Systems") +
  theme_pubr() +
  theme(axis.text.y = element_blank())

rm_pan <- def_summ %>%
  filter(Depth == "Plate") %>%
  ggplot(aes(x = Type, y = rm_mean, fill = Type))+
  geom_col_pattern(aes(pattern = Type), 
                   position = position_dodge(width = 0.9), 
                   colour = "black", , pattern_fill = "black") +  
  geom_errorbar(aes(ymin = rm_mean - rm_sd, ymax = rm_mean + rm_sd, group = Type), 
                position = position_dodge(width = 0.9),
                width = 0.25, linewidth = 1) +
  stat_pvalue_manual(rm_tukey_plate,
                     coord.flip = TRUE, hide.ns = FALSE, tip.length = 0.01)+
  coord_flip()+
  facet_wrap(~ Lake, ncol = 1) +
  scale_pattern_manual(values = c("none", "circle")) + 
  scale_fill_manual(values = c("gray", "white"))+
  labs(y = "Mean Relative Abundance/Activity",
       title = "MAGs w/ RM Systems")+
  theme_pubr() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

cas_tukey <- def_long %>% 
  filter(System == "Cas") %>%
  filter(LakeDepth != "MahoneyBottom") %>% 
  group_by(Lake, Depth, Type, sample) %>%
  summarise(total = sum(Amount)) %>%
  ungroup() %>%
  group_by(Lake, Depth) %>%
  mutate(Depth = factor(Depth, levels = c("Surface", "Plate", "Bottom"))) %>%
  tukey_hsd(total ~ Type) %>%
  add_significance()

cas_tukey <- cas_tukey %>% add_xy_position(x = "Depth", scales = "fixed", step.increase = 0.05)
cas_tukey <- cas_tukey %>% group_by(group1, group2) %>% 
  mutate(y.position = 0.99) %>% 
  as.data.frame()
cas_tukey$Depth <- factor(cas_tukey$Depth, levels = c("Surface", "Plate", "Bottom"))

cas_tukey_plate <- def_long %>% 
  filter(System == "RM") %>%
  filter(Depth == "Plate") %>% 
  group_by(Lake, Depth, Type, sample) %>%
  summarise(total = sum(Amount)) %>%
  ungroup() %>%
  group_by(Lake, Depth) %>%
  tukey_hsd(total ~ Type) %>%
  add_significance()

cas_tukey_plate <- cas_tukey_plate %>% 
  add_xy_position(x = "Type", scales = "fixed")
cas_tukey_plate <- cas_tukey_plate %>% 
  group_by(group1, group2) %>% 
  mutate(y.position = max(y.position)) %>% 
  as.data.frame()

cas_full <- ggplot(data = def_summ, aes(x = Depth, y = cas_mean, fill = Depth))+
  geom_col_pattern(aes(pattern = Type), 
                   position = position_dodge(width = 0.9), 
                   colour = "black", pattern_fill = "black") +  
  geom_errorbar(aes(ymin = cas_mean - cas_sd, ymax = cas_mean + cas_sd, group = Type), 
                position = position_dodge(width = 0.9),
                width = 0.25) +
  stat_pvalue_manual(cas_tukey,  coord.flip = TRUE, hide.ns = FALSE, tip.length = 0.01)+
  coord_flip()+
  facet_wrap(~ Lake, ncol = 1) +
  scale_pattern_manual(values = c("none", "stripe")) + 
  scale_fill_manual(values = c("darkgoldenrod", "violetred", "paleturquoise3"))+
  labs(y = "Mean Relative Abundance/Activity",
       x = "MAGs with Cas Systems") +
  theme_pubr() +
  theme(axis.text.y = element_blank())

cas_pan <- def_summ %>%
  filter(Depth == "Plate") %>%
  ggplot(aes(x = Type, y = cas_mean, fill = Type))+
  geom_col_pattern(aes(pattern = Type), 
                   position = position_dodge(width = 0.9), 
                   colour = "black", , pattern_fill = "black") +  
  geom_errorbar(aes(ymin = cas_mean - cas_sd, ymax = cas_mean + cas_sd, group = Type), 
                position = position_dodge(width = 0.9),
                width = 0.25, linewidth = 1) +
  stat_pvalue_manual(cas_tukey_plate,
                     coord.flip = TRUE, hide.ns = FALSE, tip.length = 0.01)+
  coord_flip()+
  facet_wrap(~ Lake, ncol = 1) +
  scale_pattern_manual(values = c("none", "circle")) + 
  scale_fill_manual(values = c("gray", "white"))+
  labs(y = "Mean Relative Abundance/Activity",
       title = "MAGs w/ Cas Systems")+
  theme_pubr() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))

ggarrange(do_depth, micro_pan, srb_pan, nrow = 1, widths = c(2.5,3,3))
ggarrange(lyt_pan, temp_pan, cas_pan, rm_pan, nrow = 1,
          common.legend = T, legend = "bottom")

micro_leg <- summ_long %>%
  ggplot(aes(x = anox_grp, y = mean, fill = anox_grp))+
  geom_col() + 
  scale_fill_manual(name = "Anoxygenic phototrophs", values  = c("slategray", "green4", "orchid3"))
micro_leg <- get_legend(micro_leg) %>% as_ggplot()
abio_leg <- ggplot(meta, aes(x = Depth, y = DO))+
  geom_line(aes(color = "DO (%)"), linewidth = 1)+
  geom_line(data = meta[!is.na(meta$PAR), ], aes(y = (PAR/4.5), color = "PAR"), linewidth = 1)+
  scale_y_continuous(name = "DO (%)", sec.axis = sec_axis(~ .*4.5, name = "PAR"))+
  geom_vline(aes(xintercept = Surface), colour ="paleturquoise3", lwd = 0.75, linetype = 2)+
  geom_vline(aes(xintercept = Chemocline), colour ="violetred", lwd = 0.75, linetype = 2)+
  geom_vline(aes(xintercept = Bottom), colour ="darkgoldenrod", lwd = 0.75, linetype = 2)+
  scale_color_manual(name = "Conditions", values = c(Surface = "paleturquoise3", Plate = "violetred", Bottom = "darkgoldenrod",
                                                  `DO (%)` = "black", PAR = "yellowgreen"))

abio_leg <- get_legend(abio_leg) %>% as_ggplot()
def_leg <- def_summ %>%
  ggplot(aes(x = Type, y = rm_mean, fill = Type))+
  geom_col_pattern(aes(pattern = Type), 
                   position = position_dodge(width = 0.9), 
                   colour = "black", , pattern_fill = "black") +  
  scale_pattern_manual(values = c("none", "circle")) + 
  scale_fill_manual(values = c("gray", "white"))
def_leg <- get_legend(def_leg)  %>% as_ggplot()

ggarrange(abio_leg, micro_leg, def_leg, nrow = 1)

s <- mag_func %>% filter(grepl("K17218", KEGG))
s_bin <- merge(c2b, s, by = "contig")
s_bin <- s_bin %>% filter(bin %in% s_bin[duplicated(s_bin$bin), "bin"])
s_taxa <- merge(s_bin, taxa, by.x = "bin", by.y = "user_genome")

ggarrange(micro_full, srb_full, 
          labels = c("A", "B"), widths = c(2, 1),
          font.label = list(size = 20))

ggarrange(lyt_full, temp_full, cas_full, rm_full, 
          labels = c("A", "B", "C", "D"),
          font.label = list(size = 20),
          common.legend = TRUE, legend = "bottom")
