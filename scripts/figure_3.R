#viral transcription vs viral abundance
# Load libraries
library(readr)      # for reading TSV
library(dplyr)      # for data wrangling
library(tidyr)      # for reshaping
library(stringr)
library(dplyr)
library(ggplot2)
#install.packages("ggalluvial") 
library(ggalluvial)
library(tidyverse)
library(scales)

#color pallete
my_colors <- c(
  "#120e99", # green
  "#377EB8", # purple
  "lightblue", # orange
  "#f7edcb",
  "#ffb15e",
  "#f0cdb1", # brown
  "#bf5328",
  "#FC8D62", # salmon
  "#009476", # teal
  "#F781BF", # pink
  "#dba1f7",
  "#63e0dc",
  "#8DA0CB", # periwinkle
  "#999999",
  "black"# gray
)


#load datasheets
#phage host pair info
hoges <- read.table("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/General - silveiralab/Manuscripts/Walker & Varona et al 2025/Data/Varona_scripts/Data/qc_votus_host_matches.tsv", 
                    header=TRUE)
glimpse(hoges)      # quick look at structure
hoges2 <- hoges %>%
  mutate(
    # Standardize to uppercase for both host and virus
    host_upper  = str_to_upper(host),
    virus_upper = str_to_upper(virus),
    
    # ---- Host parsing ----
    lake_host     = str_extract(host_upper, "MAH|POI|LB"),
    lake_host = recode(lake_host, "LB" = "LIM"),
    depth_host    = str_extract(host_upper, "(?<=_)[SCB](?=\\d)"),
    replicate_host = str_extract(host_upper, "(?<=_[SCB])(\\d)"),
    
    # ---- Virus parsing ----
    lake_virus     = str_extract(virus_upper, "MAH|POI|LIM|LB"),
    lake_virus = recode(lake_virus, "LB" = "LIM"),
    depth_virus    = str_extract(virus_upper, "(?<=_)[SCB](?=\\d)"),
    replicate_virus = str_extract(virus_upper, "(?<=_[SCB])(\\d)")
  ) %>%
  # Drop the helper uppercase columns
  select(-host_upper, -virus_upper)

#modify so the crispr hits get added
# Get virus-host pairs from iphop that are CRISPR hits

#also take iphop linkages where crispr was identified and add them to the other crispr column 
iphop <- read.table("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/General - silveiralab/Manuscripts/Walker & Varona et al 2025/Data/Varona_scripts/Data/Host_prediction_to_genome_m90.csv", header=TRUE, sep=",")

iphop_crispr_hits <- iphop %>% filter(!Confidence.score<95) %>%
  filter(Main.method == "CRISPR") %>%
  transmute(virus=Virus, host = Host.genome, crispr_hit = TRUE)



#add the iphop crispr counts do the df
hoges3 <- hoges2 %>%
  left_join(iphop_crispr_hits, by = c("virus", "host")) %>%
  mutate(
    # If either existing crispr column OR iphop says crispr → mark TRUE
    crispr = crispr | if_else(is.na(crispr_hit), FALSE, crispr_hit),
    
    # Make iphop_other column: TRUE only if iphop was TRUE but not crispr
    iphop_other = iphop & !crispr
  ) %>%
  select(-crispr_hit)  # remove helper column


#vOTU tanscriptional abundace
vOTU_trxn <- read.table("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/General - silveiralab/Manuscripts/Walker & Varona et al 2025/Data/Varona_scripts/Data/all_votu_relative_fractional_transcription.tsv",
                        header=TRUE)

# Reshape to long format
vOTU_trxn_long <- vOTU_trxn %>%
  pivot_longer(
    cols = -Virus,                  # All columns except "Bin"
    names_to = "SampleLocation",  # New column for sample location
    values_to = "viralTranscription"       # New column for abundance
  )

#vOTU relative abundace

vOTU_abund <- read.table("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/General - silveiralab/Manuscripts/Walker & Varona et al 2025/Data/Varona_scripts/Data/all_votu_relative_fractional_abundance.tsv",
                         header=TRUE)

# Reshape to long format
vOTU_abund_long <- vOTU_abund %>%
  pivot_longer(
    cols = -Virus,                  # All columns except "Bin"
    names_to = "SampleLocation",  # New column for sample location
    values_to = "viralAbundance"       # New column for abundance
  )

#bMAG relativea bundance
host_abund <- read.table("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/General - silveiralab/Manuscripts/Walker & Varona et al 2025/Data/Varona_scripts/Data/all_mag_relative_fractional_abundance.tsv",
                         header=TRUE)
# Reshape to long format
host_abund_long <- host_abund %>%
  pivot_longer(
    cols = -bin,                  # All columns except "Bin"
    names_to = "SampleLocation",  # New column for sample location
    values_to = "hostAbundance"       # New column for abundance
  )


#merge the viral transcription and abundance together:
vOTU_abund_trxn_long <- vOTU_abund_long %>%
  left_join(vOTU_trxn_long,
            by = c("Virus", "SampleLocation"),
            suffix = c("_abund", "_trxn"))

vOTU_abund_trxn_long

#now merge the host info and only keep the rows that have the host info, if no host remove from df
virus_host_abund <- vOTU_abund_trxn_long %>% 
  full_join(hoges3 %>% select(virus,host,domain,phylum,class,order,family,genus,species,prophage_blast,crispr,hic,iphop),
            by = c("Virus"="virus"),relationship = "many-to-many") %>%
  mutate(across(c(crispr, hic, prophage_blast, iphop),
                ~ ifelse(. == 0, FALSE, TRUE))) %>%
  filter(!is.na(host))

#now add host abundance
virus_host_abund2 <- virus_host_abund %>%
  left_join(host_abund_long,
            by = c("host"="bin", "SampleLocation"="SampleLocation"))

#add the host transcription as well
host_trxn <- read.table("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/General - silveiralab/Manuscripts/Walker & Varona et al 2025/Data/Varona_scripts/Data/all_mag_relative_rpkm.tsv", header=TRUE)

host_trxn_long <- host_trxn %>%
  pivot_longer(
    cols = -bin,                  # All columns except "Bin"
    names_to = "SampleLocation",  # New column for sample location
    values_to = "hosttrxn"       # New column for abundance
  )

virus_host_abund3 <- virus_host_abund2 %>%
  left_join(host_trxn_long,
            by = c("host"="bin", "SampleLocation"="SampleLocation"))

#also want to add viral lifestyle info
viral_lifestyle <- read.table("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/General - silveiralab/Manuscripts/Walker & Varona et al 2025/Data/Varona_scripts/Data/viral_lifestyle_summary.tsv", header=TRUE, sep="\t")

virus_host_abund4 <- virus_host_abund3 %>%
  left_join(viral_lifestyle%>% select("Virus","n_lysogen_ids","n_provirus_ids"),
            by = c("Virus")) %>% 
  mutate(
    lifestyle = if_else(
      n_lysogen_ids != 0 | n_provirus_ids != 0,
      "temperate",
      "lytic"
    ))


#only rows where host and virus has an abundance
virus_host_abund5 <- virus_host_abund4 %>%
  filter(!is.na(hostAbundance) & hostAbundance != 0) %>%
  filter(!is.na(viralTranscription) & viralTranscription != 0) %>% 
  filter(!is.na(viralAbundance) & viralAbundance != 0)

 # for magma color scale

ggplot(
  virus_host_abund5 %>% filter(!is.na(hostAbundance)),
  aes(
    x = viralTranscription/hosttrxn,
    y = viralAbundance/hostAbundance,
    color = log10(hostAbundance),
    shape = lifestyle# continuous color mapping 
  )
) +
  #geom_point( aes(alpha = lifestyle), size = 3) +
#scale_alpha_manual(   values = c("temperate" = 1, "lytic" = 0.4),guide = "none"  # removes alpha legend) +
  geom_point(size = 3, alpha=0.7) +
  scale_shape_manual(
    values = c("temperate" = 16, "lytic" = 17)  # 16 = circle, 17 = triangle
  ) +
  scale_color_viridis_c(option = "magma") +  # continuous magma color scale
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() +
  labs(
    color = "Host abundance",
    x = "VHR transcription",
    y = "VHR abundance"
  )


#with quartiles:
library(dplyr)
library(ggplot2)
library(scales)
library(ggforce)
virus_host_abund5$VHR_abundance <- virus_host_abund5$viralAbundance/virus_host_abund5$hostAbundance
virus_host_abund5$VHR_transcription <- virus_host_abund5$viralTranscription/virus_host_abund5$hosttrxn


virus_host_abund5$log10_VHR_abundance <- log10(virus_host_abund5$VHR_abundance)
virus_host_abund5$log10_VHR_transcription <- log10(virus_host_abund5$VHR_transcription)

virus_host_abund6 <- virus_host_abund5 %>% filter(!is.na(VHR_abundance)) %>% filter(!is.na(VHR_transcription)) %>%
  filter((VHR_abundance!=0)) %>%   filter((VHR_transcription!=0))

# Filter out NAs and zeros
virus_host_abund6 <- virus_host_abund5 %>%
  filter(!is.na(VHR_abundance), !is.na(VHR_transcription)) %>%
  filter(VHR_abundance != 0, VHR_transcription != 0)

# Compute quartile thresholds once
abund_q25 <- quantile(virus_host_abund6$log10_VHR_abundance, 0.25, na.rm = TRUE) 
abund_q75 <- quantile(virus_host_abund6$log10_VHR_abundance, 0.75, na.rm = TRUE)
trxn_q25 <- quantile(virus_host_abund6$log10_VHR_transcription, 0.25, na.rm = TRUE)
trxn_q75 <- quantile(virus_host_abund6$log10_VHR_transcription, 0.75, na.rm = TRUE)

# Classify rows based on quartile cutoffs
virus_host_abund7 <- virus_host_abund6 %>%
  mutate(
    category = case_when(
      log10_VHR_abundance >= abund_q75 & log10_VHR_transcription >= trxn_q75 ~ "cat_2",  # high-high
      log10_VHR_abundance <= abund_q25 & log10_VHR_transcription <= trxn_q25 ~ "cat_3",  # low-low
      log10_VHR_abundance >= abund_q75 & log10_VHR_transcription <= trxn_q25 ~ "cat_1",  # high-low
      log10_VHR_abundance <= abund_q25 & log10_VHR_transcription >= trxn_q75 ~ "cat_4",  # low-high
      TRUE ~ "med"
    )
  )

virus_host_abund7 <- virus_host_abund7 %>% mutate(SB = case_when(order == "o__Chromatiales" ~ "PSB",
                                                                 order == "o__Chlorobiales" ~ "GSB",
                                                                 order %in% c("o__Desulfobacterales", "o__Desulfobulbales",
                                                                              "o__Desulfatiglandales") ~ "SRB",
                                                                 .default = "Other"),
                                                  chemo = case_when(SB == "Other" ~ "non", 
                                                                    grepl("_C", SampleLocation) ~ "chemo",
                                                                    SB != "Other" ~ "med",
                                                                    .default = "non"),
                                                  SB2 = paste0(SB, "_", chemo))


ggplot(
  virus_host_abund7 %>% filter(!is.na(hostAbundance)),
  aes(
    x = VHR_transcription,
    y = VHR_abundance,
    color = SB2,
    shape = lifestyle,
    alpha = chemo
  )
) +
  geom_point(aes(size = -log10(hostAbundance))) +
  geom_vline(xintercept = c(10^trxn_q25, 10^trxn_q75), linetype = "dotted") +
  geom_hline(yintercept = c(10^abund_q25, 10^abund_q75), linetype = "dotted") +
  scale_shape_manual(values = c("temperate" = 16, "lytic" = 17)) +
  scale_color_manual(values = c("green4", "green3", "gray75", "orchid3","orchid1", "khaki2", "khaki1")) +
  scale_alpha_manual(values = c(1, 0.4, 0.1))+
  scale_size(trans = 'reverse')+
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  theme_classic() +
  labs(
    color = "Host abundance (log10)",
    x = "VHR transcription",
    y = "VHR abundance"
  )

# Summarize relative abundance per category
pie_data <- virus_host_abund7 %>%
  group_by(category, class) %>%
  summarise(total_rel = sum(hostAbundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(category) %>%
  mutate(prop = total_rel / sum(total_rel))

# Plot pies


pie_data <- virus_host_abund7 %>%
  group_by(category, class) %>%
  summarise(total_host = sum(hostAbundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(category) %>%
  mutate(
    rel_abund = total_host / sum(total_host)  # sums to 1 (or 100%) within each category
  ) %>%
  ungroup()

# Pie chart (each facet = one category)
ggplot(pie_data, aes(x = "", y = rel_abund, fill = class)) +
  geom_col(color = "white", width = 1) +
  coord_polar(theta = "y") +
  facet_wrap(~ category) +
  scale_fill_viridis_d(option = "plasma") +
  theme_void() +
  labs(
    fill = "Bacterial Class",
    title = "Distribution of Virus-Associated Classes by VHR Category"
  )


#only top 14

library(dplyr)
library(ggplot2)
library(viridis)

# Step 1: Summarize host abundance and select top 14 overall
top14_hosts <- virus_host_abund7 %>%
  group_by(class) %>%
  summarise(global_abundance = sum(hostAbundance, na.rm = TRUE)) %>%
  arrange(desc(global_abundance)) %>%
  slice_head(n = 14) %>%
  pull(class)

# Step 2: Filter and recalculate relative abundance per category
pie_data <- virus_host_abund7 %>%
  mutate(class = if_else(class %in% top14_hosts, class, "Other")) %>%
  group_by(category, class) %>%
  summarise(total_host = sum(hostAbundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(category) %>%
  mutate(rel_abund = total_host / sum(total_host)) %>%
  ungroup()
# Step 3: Plot pie charts
ggplot(pie_data, aes(x = "", y = rel_abund, fill = class)) +
  geom_col(color = "white", width = 1) +
  coord_polar(theta = "y") +
  facet_wrap(~ category) +
  theme_void() +
  labs(
    fill = "Bacterial Class",
    title = "Top 14 Virus-Associated Host Classes by VHR Category"
  ) + 
  scale_fill_manual(values = my_colors)


cor.test(virus_host_abund7$VHR_abundance, virus_host_abund7$VHR_transcription, 
         method = 'spearman')
