library(dplyr)
library(tidyr)
library(stringi)
library(compositions)


setwd("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/input_data/")
omics <- read.delim("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/omics_read_counts.tsv")
c2b <- read.delim("hosts/contigs2bins.tsv", col.names = c("bin", "Contig"))

#vOTU fractional abundance calculations
dna_map_long <- read.delim("coverm_votu_qc_DNA.tsv")
dna_map_long <- dna_map_long %>% pivot_longer(cols = LB_B1.Length:POI_S5.Covered.Fraction,
                                              names_to = "temp", values_to = "num")
dna_map_long <- cbind(dna_map_long, stri_split_fixed(dna_map_long$temp, ".", n = 2, simplify = TRUE))
colnames(dna_map_long) <- c("Contig", "temp", "num", "sample", "metric")
dna_map_long <- dna_map_long %>% select(-temp) 
dna_map_long <- dna_map_long %>% pivot_wider(names_from = "metric", values_from = "num")
dna_map_long <- dna_map_long %>% 
#  mutate(Read.Count = case_when(Covered.Fraction < 0.05 ~ 0, .default = Read.Count)) %>% 
  select(-Covered.Fraction)
dna_map <- dna_map_long %>% pivot_wider(names_from = "sample", values_from = "Read.Count")
dna_mean_len <- dna_map_long %>% filter(Read.Count != 0) %>% group_by(sample) %>% summarise(mean_len = mean(Length))

dna_frac <- sweep(dna_map[ , 3:45], 2, dna_mean_len$mean_len, "*")
dna_frac <- sweep(dna_frac, 2, omics$reads, "/")
dna_frac <- sweep(dna_frac, 1, dna_map$Length, "/")
dna_per <- prop.table(as.matrix(dna_frac), 2) %>% as.data.frame()
dna_per <- data.frame(Virus = dna_map$Contig, dna_per)
dna_frac <- data.frame(Virus = dna_map$Contig, dna_frac)

dna_clr <- clr(dna_map[ , 3:45])
dna_clr <- data.frame(Virus = dna_map$Contig, dna_clr)

write.table(dna_map[, -2], 
            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/all_votu_read_counts.tsv",
            sep = "\t", row.names = F, quote = F)
write.table(dna_per, 
            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/all_votu_relative_fractional_abundance.tsv",
            sep = "\t", row.names = F, quote = F)
write.table(dna_frac, 
            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/all_votu_fractional_abundance.tsv",
            sep = "\t", row.names = F, quote = F)
write.table(dna_clr, 
            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/all_votu_clr_abundance.tsv",
            sep = "\t", row.names = F, quote = F)

#vOTU RPKM 
rna_raw <- read.delim("coverm_votu_qc_RNA.tsv")
rna_map_long <- rna_raw %>% pivot_longer(cols = LB_B1.Length:POI_S4.Covered.Fraction,
                                              names_to = "temp", values_to = "num")
rna_map_long <- cbind(rna_map_long, stri_split_fixed(rna_map_long$temp, ".", n = 2, simplify = TRUE))
colnames(rna_map_long) <- c("Contig", "temp", "num", "sample", "metric")
rna_map_long <- rna_map_long %>% select(-temp)
rna_map_long <- rna_map_long %>% pivot_wider(names_from = "metric", values_from = "num")
rna_map <- rna_map_long %>% select(-c(Covered.Fraction)) %>% 
  pivot_wider(names_from = "sample", values_from = "Read.Count")
rna_mean_len <- rna_map_long %>% 
  filter(Read.Count != 0) %>% 
  group_by(sample) %>% 
  summarise(mean_len = mean(Length))

rna_frac <- rna_map[ , 3:35] * rna_mean_len$mean_len
rna_frac <- sweep(rna_frac, 2, omics[!is.na(omics$transcripts), "reads"], "/")
rna_frac <- sweep(rna_frac, 1, rna_map$Length, "/")
rna_per <- prop.table(as.matrix(rna_frac), 2) %>% as.data.frame()
rna_per <- data.frame(Virus = rna_map$Contig, rna_per)
rna_frac <- data.frame(Virus = rna_map$Contig, rna_frac)

rna_clr <- clr(rna_map[ , 3:35])
rna_clr <- data.frame(Virus = rna_map$Contig, rna_clr)

write.table(rna_map[, -2], 
            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/all_votu_transcript_counts.tsv",
            sep = "\t", row.names = F, quote = F)
write.table(rna_per, 
            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/all_votu_relative_fractional_transcription.tsv",
            sep = "\t", row.names = F, quote = F)
write.table(rna_frac, 
            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/all_votu_fractional_transcription.tsv",
            sep = "\t", row.names = F, quote = F)
write.table(rna_clr, 
            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/all_votu_clr_transcription.tsv",
            sep = "\t", row.names = F, quote = F)
#vOTU cds DNA 
dna_cds_raw <- read.delim("coverm_votu_cds_DNA.tsv")
colnames(dna_cds_raw) <- gsub("_cds", "", colnames(dna_cds_raw))
dna_cds_map_long <- dna_cds_raw %>% pivot_longer(cols = LB_B1.Length:POI_S5.Covered.Fraction,
                                                 names_to = "temp", values_to = "num")
dna_cds_map_long <- cbind(dna_cds_map_long, stri_split_fixed(dna_cds_map_long$temp, ".", n = 2, simplify = TRUE))
colnames(dna_cds_map_long) <- c("Contig", "temp", "num", "sample", "metric")
dna_cds_map_long <- dna_cds_map_long %>% select(-temp)
dna_cds_map_long <- dna_cds_map_long %>% pivot_wider(names_from = "metric", values_from = "num")
dna_cds_map <- dna_cds_map_long %>% select(-c(Covered.Fraction)) %>% 
  pivot_wider(names_from = "sample", values_from = "Read.Count")

dna_cds_frac <- dna_cds_map[ , 3:45]
dna_cds_frac <- sweep(dna_cds_frac, 2, omics[!is.na(omics$reads), "reads"], "/")
dna_cds_frac <- sweep(dna_cds_frac, 1, dna_cds_map$Length, "/")
dna_cds_per <- dna_cds_frac * 10^9 
dna_cds_per <- data.frame(Virus = gsub(" .*", "", dna_cds_map$Contig), dna_cds_per)

write.table(dna_cds_per, 
            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/all_votu_rpkm_dna_cds.tsv",
            sep = "\t", row.names = F, quote = F)
#vOTU cds RPKM 
rna_cds_raw <- read.delim("coverm_votu_cds_RNA.tsv")
colnames(rna_cds_raw) <- gsub("_cds", "", colnames(rna_cds_raw))
colnames(rna_cds_raw) <- gsub("LIMB", "LB", colnames(rna_cds_raw))
rna_cds_map_long <- rna_cds_raw %>% pivot_longer(cols = LB_B1.Length:POI_S4.Covered.Fraction,
                                         names_to = "temp", values_to = "num")
rna_cds_map_long <- cbind(rna_cds_map_long, stri_split_fixed(rna_cds_map_long$temp, ".", n = 2, simplify = TRUE))
colnames(rna_cds_map_long) <- c("Contig", "temp", "num", "sample", "metric")
rna_cds_map_long <- rna_cds_map_long %>% select(-temp)
rna_cds_map_long <- rna_cds_map_long %>% pivot_wider(names_from = "metric", values_from = "num")
rna_cds_map <- rna_cds_map_long %>% select(-c(Covered.Fraction)) %>% 
  pivot_wider(names_from = "sample", values_from = "Read.Count")

rna_cds_frac <- rna_cds_map[ , 3:35] 
rna_cds_frac <- sweep(rna_cds_frac, 2, omics[!is.na(omics$transcripts), "reads"], "/")
rna_cds_frac <- sweep(rna_cds_frac, 1, rna_cds_map$Length, "/")
rna_cds_per <- rna_cds_frac * 10^9 
rna_cds_per <- data.frame(Virus = gsub(" .*", "", rna_cds_map$Contig), rna_cds_per)

write.table(rna_cds_per, 
            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/all_votu_rpkm_rna_cds.tsv",
            sep = "\t", row.names = F, quote = F)

#MAG fractional abundance
mag_dna_map <- read.delim("hosts/coverm_host_DNA.tsv")
host_dna <- mag_dna_map %>% 
  select(Contig, contains("Read.Count"))
colnames(host_dna) <- gsub(".Read.Count", "", colnames(host_dna))
colnames(host_dna) <- gsub("dna_map", "", colnames(host_dna))

mag_dna <- merge(host_dna, c2b, by = "Contig")
mag_dna_sum <- mag_dna %>% 
  group_by(bin) %>% 
  summarise_if(is.numeric, sum) %>%
  arrange(bin)
mag_info <- read.delim("hosts/genomeInformation.csv", sep = ",")
mag_info$genome <- gsub(".fna", "", mag_info$genome)
mag_info <- mag_info[mag_info$genome %in% mag_dna_sum$bin, ] %>%
  arrange(genome)
mag_mean_len <- mean(mag_info$length)

mag_dna_frac <- mag_dna_sum[ , 2:44] * mag_mean_len
mag_dna_frac <- sweep(mag_dna_frac, 2, omics$reads, "/")
mag_dna_frac <- sweep(mag_dna_frac, 1, mag_info$length, "/")
mag_dna_per <- prop.table(as.matrix(mag_dna_frac), 2) %>% as.data.frame()
mag_dna_per <- data.frame(bin = mag_dna_sum$bin, mag_dna_per)
mag_dna_frac <- data.frame(bin = mag_dna_sum$bin, mag_dna_frac)

mag_dna_clr <- clr(mag_dna_sum[ , 2:44])
mag_dna_clr <- data.frame(bin = mag_dna_sum$bin, mag_dna_clr)

write.table(mag_dna_per, 
          "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/all_mag_relative_fractional_abundance.tsv",
           sep = "\t", row.names = F, quote = F)
write.table(mag_dna_frac, 
            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/all_mag_fractional_abundance.tsv",
            sep = "\t", row.names = F, quote = F)
write.table(mag_dna_clr, 
            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/all_mag_clr_abundance.tsv",
            sep = "\t", row.names = F, quote = F)

#MAG fractional rpkm
mag_rna_map <- read.delim("hosts/coverm_host_RNA.tsv")
host_rna <- mag_rna_map %>% 
  select(Contig, contains("Read.Count"))
colnames(host_rna) <- gsub(".Read.Count", "", colnames(host_rna))
colnames(host_rna) <- gsub("rna_map", "", colnames(host_rna))
colnames(host_rna) <- gsub("\\.", "_", colnames(host_rna))
colnames(host_rna) <- gsub("LIMB", "LB", colnames(host_rna))
colnames(host_rna) <- gsub("MAH_C4_v2", "MAH_C2", colnames(host_rna))
host_rna <- host_rna[ , c("Contig", omics[omics$sample %in% colnames(host_rna), "sample"])]

mag_rna <- merge(host_rna, c2b, by = "Contig")
mag_rna_sum <- mag_rna %>% 
  group_by(bin) %>% 
  summarise_if(is.numeric, sum) %>%
  arrange(bin)
mag_info <- mag_info[mag_info$genome %in% mag_rna_sum$bin, ] %>%
  arrange(genome)

mag_rna_frac <- mag_rna_sum[ , 2:34] * mag_mean_len
mag_rna_frac <- sweep(mag_rna_frac, 2, omics[!is.na(omics$transcripts), "reads"], "/")
mag_rna_frac <- sweep(mag_rna_frac, 1, mag_info$length, "/")
mag_rna_per <- prop.table(as.matrix(mag_rna_frac), 2) %>% as.data.frame()
mag_rna_per <- data.frame(bin = mag_rna_sum$bin, mag_rna_per)
mag_rna_frac <- data.frame(bin = mag_rna_sum$bin, mag_rna_frac)

mag_rna_clr <- clr(mag_rna_sum[ , 2:34])
mag_rna_clr <- data.frame(bin = mag_rna_sum$bin, mag_rna_clr)

write.table(mag_rna_per, 
            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/all_mag_relative_transcription.tsv",
           sep = "\t", row.names = F, quote = F)
write.table(mag_rna_frac, 
            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/all_mag_fractional_transcription.tsv",
            sep = "\t", row.names = F, quote = F)
write.table(mag_rna_clr, 
            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/figures/data/all_mag_clr_transcription.tsv",
            sep = "\t", row.names = F, quote = F)

all_abund <- merge(mag_dna_per, mag_rna_per, by = "bin", suffixes = c("_dna", "_rna"))

write.table(all_abund, 
            "/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/input_data/hosts/all_mag_all_abundance_percent.tsv",
            sep = "\t", row.names = F, quote = F)
