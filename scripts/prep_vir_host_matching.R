library(dplyr)
library(stringr)
library(tidyr)
library(ggpubr)
library(ggsci)
library(ggplot2)


setwd("/Users/jordanwalker/Library/CloudStorage/OneDrive-UniversityofMiami/Research/SPACE_2023/viruses/")

# import vOTUs ------------------------------------------------------------
votu <- read.delim("input_data/votu_clusters.tsv", header = F, col.names = c("rep", "members"))
qc_vir <- data.frame(Virus = readLines("input_data/all_qc_viruses.txt")) #Read in the QC virus names
qc_vir$orig_contig <- gsub("_._checkv.*", "", qc_vir$Virus)


qc_votu <- merge(votu, qc_vir, by.x = "rep", by.y = "orig_contig") %>% #Filter the votus based on QC viruses
  select(-rep) %>%
  rename(rep = Virus)
qc_votu <- qc_votu %>% 
  mutate(members = strsplit(members, ",")) %>% #Make a long dataframe for easier merging
  unnest(members)

vitap_all <- read.delim("input_data/vitap_all_lineages.tsv")
vitap_best <- read.delim("input_data/vitap_best_determined_lineages.tsv")
vitap_best <- vitap_best[vitap_best$Genome_ID %in% qc_votu$rep, ]

vitap_left <- vitap_all[!vitap_all$Genome_ID %in% vitap_best$Genome_ID, ]
vitap_left_gp <- vitap_left %>%
  group_by(Genome_ID) %>%
  filter(lineage_score == max(lineage_score, na.rm = TRUE)) %>%
  ungroup()
vitap_left_gp <- vitap_left_gp %>% filter(grepl("MAH|POI|LB|LIMB|checkv", Genome_ID))
vitap_left_gp <- vitap_left_gp %>% filter(lineage_score >= 0.1)
vitap_final <- rbind(vitap_best[ , 1:2], vitap_left_gp[ , 1:2])
vitap_final$lineage <- gsub("-", "NA", vitap_final$lineage)
vitap_final <- vitap_final %>% separate(lineage, into = c("species", "genus", "family", "order", "class", "phylum", "kingdom", "realm"), 
                         sep = ";", convert = T)
votu_taxa <- merge(qc_vir, vitap_final, by.x = "Virus", by.y = "Genome_ID", all.x = T) %>%
  select(-orig_contig) %>% unique()

# import host information -------------------------------------------------

mag_clus <- read.delim("input_data/hosts/mag_clusters.tsv")
mag_taxa <- read.delim("input_data/hosts/mag_taxonomy.tsv")
mags <- merge(mag_clus, mag_taxa, by.x = "reps", by.y = "user_genome")
c2b <- read.delim("input_data/hosts/contigs2bins.tsv", col.names = c("mag", "Contig"))
rep_contigs <- c2b[c2b$mag %in% mag_taxa$user_genome, ]

# import virus host matching -------------------------------------------------------
blast_nm <- c("qseqid", "qlen", "sseqid", "slen", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
crispr_raw <- read.table("input_data/crispr_votu_noderep_100cov_align20bp_2mis.tsv", sep = "\t", col.names = blast_nm)
crispr_raw <- crispr_raw %>% mutate(contig = gsub("_CRISPR.*", "", qseqid))
crispr_raw <- crispr_raw[crispr_raw$contig %in% rep_contigs$Contig, ]

hic_raw <- data.frame(seq_name = readLines("input_data/hic_virus_contigs.txt"))
hic_filt <- data.frame(seq_name = hic_raw[hic_raw$seq_name %in% gsub("_1_checkv_.*", "", qc_votu$rep), ])

prophage_raw <- read.table("input_data/prophage_quality_bins_no_derep_filtered_blastn_95ID_100cov_500bpflank.tsv", 
                           sep = "\t", col.names = blast_nm)
prophage_raw <- merge(c2b, prophage_raw, by.x = "Contig", by.y = "sseqid")

iphop_raw <- read.table("input_data/Host_prediction_to_genome_m90.csv", sep = ",", header = T)
iphop_filt <- iphop_raw%>% filter(Confidence.score >= 95)

links <- list()
links[["crispr"]] <- data.frame(host = gsub(".*bin_", "", crispr_raw$qseqid),
                           virus = crispr_raw$sseqid)
links[["crispr"]]$host <- gsub("_CRISPR.*", "", links[["crispr"]]$host)

links[["hic"]] <- data.frame(host = paste("metator", 
                                     tolower(str_match(hic_filt$seq_name, "[MLP].*_[CBS].")), 
                                     gsub("_[MPL].*", "", hic_filt$seq_name),
                                     sep = "_"),
                           virus = hic_filt$seq_name)

links[["prophage_blast"]] <- data.frame(host = prophage_raw$mag,
                             virus = prophage_raw$qseqid)
links[["iphop"]] <- data.frame(host = iphop_filt$Host.genome,
                          virus = iphop_filt$Virus)
links[["iphop"]] <- links[["iphop"]][grepl("LIMB|POI|MAH|LB", links[["iphop"]]$host), ]

vir_hosts_long <- do.call("rbind", links)
vir_hosts_long$method <- gsub("\\.[0-9]*", "", rownames(vir_hosts_long))
vir_hosts_long <- unique(vir_hosts_long)

vir_hosts <- vir_hosts_long %>% mutate(val = 1) %>% 
  pivot_wider(names_from = method, values_from = val)
vir_hosts[is.na(vir_hosts)] <- 0
table(apply(vir_hosts[,3:6], 1, sum))
qc_vir_hosts <- merge(qc_votu, vir_hosts, by.x = "members", by.y = "virus") %>%
  select(-members) 
colnames(qc_vir_hosts)[1] <- "virus"

qc_vir_hosts <- merge(qc_vir_hosts, mags, by.x = "host", by.y = "members")
ag_qc_vir_hosts <- qc_vir_hosts %>%
  group_by_if(is.character) %>% 
  summarise_if(is.numeric, sum)
colnames(ag_qc_vir_hosts)[1:3] <- c("original_host", "virus", "rep_host")

qc_vh_long <- qc_vir_hosts %>% select(-c(reps)) %>% 
  pivot_longer(crispr:iphop, values_to = "match", names_to = "method") %>%
  select(-c(domain:species)) %>%
  filter(match > 0)

cyto <- merge(qc_votu, qc_vh_long, by.x = "members", by.y = "virus") %>%
  select(-c(members, match)) %>% unique ()
colnames(cyto)[1:2] <- c("virus", "host")

cyto_host <- data.frame(node = unique(cyto$host))
cyto_host <- merge(cyto_host, mag_taxa, by.x = "node", by.y = "user_genome") %>% 
  mutate(taxa = paste(class, order, sep = ","), vir_host = "host") %>%
  select(node, taxa, vir_host)
cyto_vir <- data.frame(node = unique(cyto$vir), seq_name = gsub("_ext.*|_rna", "", unique(cyto$vir)))
cyto_vir <- merge(cyto_vir, votu_taxa, by.x = "seq_name", by.y = "Virus") %>%
  mutate(taxa = paste(class, order, sep = ","), vir_host = "virus") %>%
  select(node, taxa, vir_host)

cyto_taxa <- rbind(cyto_host, cyto_vir)

write.table(ag_qc_vir_hosts, "figures/data/qc_votus_host_matches.tsv", sep = "\t", row.names = F, quote = F)
write.table(qc_votu, "figures/data/votu_clusters_long.tsv", sep = "\t", row.names = F, quote = F)
write.table(votu_taxa, "figures/data/qc_votu_taxonomy.tsv", sep = "\t", row.names = F, quote = F)

write.table(cyto, "figures/data/cyto_links.tsv", sep = "\t", row.names = F, quote = F)
write.table(cyto_taxa, "figures/data/cyto_taxa.tsv", sep = "\t", row.names = F, quote = F)

cyto_order <- merge(cyto, mag_taxa, by.x = "host", by.y = "user_genome") %>%
  mutate(taxa = paste(class, order, sep = ",")) %>%
  select(virus, taxa) %>% unique()
colnames(cyto_order) <- c("virus", "host")
write.table(cyto_order, "figures/data/cyto_order_links.tsv", sep = "\t", row.names = F, quote = F)

cyto_host <- data.frame(node = unique(cyto_order$host))
cyto_host$taxa <- gsub(",.*", "", cyto_host$node)
cyto_host <- merge(cyto_host, mag_taxa, by.x = "taxa", by.y = "class") %>% 
  mutate(taxa = paste(phylum, taxa, sep = ",")) %>% 
  select(node, taxa)
cyto_host$vir_host <- "host"
cyto_vir <- data.frame(node = unique(cyto$vir), seq_name = gsub("_ext.*|_rna", "", unique(cyto$vir)))
cyto_vir <- merge(cyto_vir, votu_taxa, by.x = "seq_name", by.y = "Virus") %>%
  mutate(taxa = paste(class, order, sep = ","), vir_host = "virus") %>%
  select(node, taxa, vir_host)
cyto_taxa <- rbind(cyto_host, cyto_vir)
write.table(cyto_taxa, "figures/data/cyto_order_taxa.tsv", sep = "\t", row.names = F, quote = F)

multis <- ag_qc_vir_hosts[duplicated(ag_qc_vir_hosts$virus), 2]
multi_infect <- cyto[cyto$virus %in% multis$virus, ]
multi_infect <- merge(multi_infect, mag_taxa, by.x = "host", by.y = "user_genome")

host_spec <- multi_infect %>% group_by(virus) %>% 
  summarise(how = paste(unique(method), collapse = ","), across(domain:genus, n_distinct))
host_spec <- host_spec %>% rowwise() %>% 
  filter(sum(domain, phylum, class, order, family, genus) > 6)
host_spec_long <- host_spec %>% 
  pivot_longer(domain:genus, names_to = "rank", values_to = "num_host") %>% 
  filter(num_host > 1) %>% 
  group_by(how, rank, num_host) %>%
  summarise(count = n())
host_spec_long$num_host <- as.factor(host_spec_long$num_host)
host_spec_long$rank <- factor(host_spec_long$rank, levels = c("domain", "phylum", "class", "order", "family", "genus"))
host_spec_long[host_spec_long$how == "iphop,hic", c("how")] <- "hic"
ggplot(host_spec_long, aes(fill=how, y=count, x=rank)) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(~ num_host)+
  scale_fill_futurama()

