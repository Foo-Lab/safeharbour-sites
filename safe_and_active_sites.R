pacman::p_load(tidyverse, janitor, skimr, biomaRt, GenomicFeatures,
               GenomicRanges, openxlsx, Repitools)
rm(list=ls())


hg38_chrom_lenght <- getChromInfoFromBiomart(dataset = "hsapiens_gene_ensembl", 
                                             host = "https://asia.ensembl.org")

chr_ID <- c((1:22), "X", "Y")

hg38_chrom_lenght <- hg38_chrom_lenght %>% subset(chrom %in% chr_ID)

seq_info <- Seqinfo(seqnames = chr_ID, seqlengths = hg38_chrom_lenght$length,
                    isCircular = c(rep(FALSE, 24)), genome = "hg38")

#########active regions

active_regions <- read.xlsx("HiC active compartments.xlsx")

active_regions$coordinates <- paste("chr", active_regions$Chromosome,
                                    ":", active_regions$Start_HiC_CompartmentA,
                                    ":", active_regions$End_HiC_CompartmentA,
                                    sep = "")

active_regions <- active_regions %>% mutate(coordinates = as.factor(coordinates))

active_regions <- distinct(active_regions, coordinates, .keep_all = TRUE)

active_regions <- active_regions %>% 
  dplyr::filter(`Fraction_of_HiC_samples_in_A.(21.in.total)` > 0.95) %>% 
  dplyr::filter(`Overlap_Compartment-Gene` == "full_overlap") %>% 
  dplyr::rename(seqnames = Chromosome, start = Start_HiC_CompartmentA,
                end = End_HiC_CompartmentA)

active_ranges <- makeGRangesFromDataFrame(active_regions,ignore.strand=T,
                                          keep.extra.columns=T, seqinfo = seq_info)

active_ranges <- trim(active_ranges)

length(active_ranges)


#####safe & active

hg38_full_filter <- read.xlsx("hg38_full_filter.xlsx")

hg38_full_filter_lnc100k_ranges <- makeGRangesFromDataFrame(hg38_full_filter,ignore.strand=T,
                                                            seqinfo = seq_info)

hits_safe_and_active_100_lnc <- findOverlaps(hg38_full_filter_lnc100k_ranges, active_ranges, 
                                             type = "any", maxgap = -1)

safe_100lnc_hits <- hg38_full_filter_lnc100k_ranges[queryHits(hits_safe_and_active_100_lnc)]
active_100lnc_hits <- active_ranges[subjectHits(hits_safe_and_active_100_lnc)]

safe_100lnc_hits_df <- annoGR2DF(safe_100lnc_hits)
active_100lnc_hits_df <- annoGR2DF(active_100lnc_hits)

safe_and_active_100lnc_hits_df <- active_100lnc_hits_df %>% dplyr::rename(
  active_start = start, active_end = end, active_width = width) %>% 
  mutate(start = safe_100lnc_hits_df$start, .before = active_start) %>% 
  mutate(end = safe_100lnc_hits_df$end, .before = active_start) %>% 
  mutate(width = safe_100lnc_hits_df$width, .before = active_start)

safe_and_active_100lnc_hits_df$locus <- paste(">", "chr", safe_and_active_100lnc_hits_df$chr,
                                              ":", safe_and_active_100lnc_hits_df$start,
                                              ":", safe_and_active_100lnc_hits_df$end, sep = "")

safe_and_active_100lnc_hits_df <- distinct(safe_and_active_100lnc_hits_df,
                                           locus, .keep_all = TRUE)

write.xlsx(safe_and_active_100lnc_hits_df,
           "safe_and_active_ranges_0.95active.xlsx")
