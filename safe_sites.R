pacman::p_load(tidyverse, janitor, skimr, biomaRt, GenomicFeatures,
               GenomicRanges, openxlsx, Repitools)
rm(list=ls())

####inputs

hg38_chrom_lenght <- getChromInfoFromBiomart(dataset = "hsapiens_gene_ensembl", 
                                             host = "https://asia.ensembl.org")

ensembl <- useEnsembl('genes',
                      dataset = "hsapiens_gene_ensembl",
                      version = 103)

#lista <- listAttributes(ensembl)


gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position",
                                  "end_position", "gene_biotype", "strand", "band", "ensembl_transcript_id",
                                  "transcript_start", "transcript_end", "transcription_start_site",
                                  "transcript_biotype", "description", "ensembl_transcript_id_version",
                                  "ensembl_gene_id_version"), mart = ensembl)

chr_ID <- c((1:22), "X", "Y")

gene_info <- gene_info %>% subset(chromosome_name %in% chr_ID)
hg38_chrom_lenght <- hg38_chrom_lenght %>% subset(chrom %in% chr_ID)

seq_info <- Seqinfo(seqnames = chr_ID, seqlengths = hg38_chrom_lenght$length,
                    isCircular = c(rep(FALSE, 24)), genome = "hg38")

rm(ensembl)

###chromosome ranges
hg38_ranges <- GRanges(seqnames = hg38_chrom_lenght$chrom, 
                       ranges = IRanges(start = 1,width = hg38_chrom_lenght$length), 
                       seqinfo = seq_info)

head(hg38_ranges)
rm(hg38_chrom_lenght)

###DNAse ENCODE3 meuleman2020 data https://www.encodeproject.org/files/ENCFF503GCK/
DHS_meuleman2020 <- read.delim("ENCFF503GCK.tsv")
temp <- str_split_fixed(DHS_meuleman2020$seqname, "r", 2)
temp <- as.data.frame(temp)
DHS_meuleman2020 <- DHS_meuleman2020 %>% mutate(seqname = temp$V2)
rm(temp)

dnase_ranges <- GRanges(seqnames= DHS_meuleman2020$seqname, 
                        IRanges(start= DHS_meuleman2020$start,end=DHS_meuleman2020$end),
                        seqinfo = seq_info)
head(dnase_ranges)
rm(DHS_meuleman2020)

dnase_ranges <- dnase_ranges + 2000


####UCR
UCR <- read.delim("ucr_hg38.txt")
temp <- str_split_fixed(UCR$seqname, "r", 2)
temp <- as.data.frame(temp)
UCR <- UCR %>% mutate(seqname = temp$V2)
ucr_ranges <- makeGRangesFromDataFrame(UCR,ignore.strand=T, seqinfo = seq_info)

head(ucr_ranges)

rm(temp, UCR)

####COSMIC

cancer_consensus <- read.xlsx("COSMIC.xlsx")
cancer_genes <- unique(cancer_consensus$Gene.Symbol)

cancer_coordinates <- gene_info %>% subset(external_gene_name %in% cancer_genes) %>% 
  dplyr::select(chromosome_name, ensembl_gene_id, external_gene_name, transcription_start_site) %>% 
  dplyr::rename(seqname = chromosome_name)

TSS_1 <- cancer_coordinates$transcription_start_site +1
cancer_coordinates <- cancer_coordinates %>% mutate(TSS_1 = TSS_1) 
cancer_coordinates <- cancer_coordinates%>% dplyr::rename(start = transcription_start_site) %>% 
  dplyr::rename(end = TSS_1)

cancer_ranges <- makeGRangesFromDataFrame(cancer_coordinates,ignore.strand=T,
                                          keep.extra.columns=T, seqinfo = seq_info)
cancer_ranges <- cancer_ranges + 300000

cancer_ranges <- trim(cancer_ranges)


head(cancer_ranges)

rm(cancer_consensus, cancer_coordinates, cancer_genes, TSS_1)

####miRNA

miRNA_coordinates <- gene_info %>% subset(transcript_biotype == "miRNA") %>% 
  dplyr::select(chromosome_name, ensembl_gene_id, external_gene_name, transcription_start_site) %>% 
  dplyr::rename(seqname = chromosome_name)

TSS_1 <- miRNA_coordinates$transcription_start_site +1
miRNA_coordinates <- miRNA_coordinates %>% mutate(TSS_1 = TSS_1) 
miRNA_coordinates <- miRNA_coordinates%>% dplyr::rename(start = transcription_start_site) %>% 
  dplyr::rename(end = TSS_1)

miRNA_ranges <- makeGRangesFromDataFrame(miRNA_coordinates,ignore.strand=T,
                                         keep.extra.columns=T, seqinfo = seq_info)
miRNA_ranges <- miRNA_ranges + 300000

miRNA_ranges <- trim(miRNA_ranges)

head(miRNA_ranges)

rm(miRNA_coordinates, TSS_1)

###lncRNA

lncRNA_coordinates <- gene_info %>% subset(transcript_biotype == "lncRNA") %>% 
  dplyr::select(chromosome_name, ensembl_gene_id, external_gene_name, transcription_start_site) %>% 
  dplyr::rename(seqname = chromosome_name)

TSS_1 <- lncRNA_coordinates$transcription_start_site +1
lncRNA_coordinates <- lncRNA_coordinates %>% mutate(TSS_1 = TSS_1) 
lncRNA_coordinates <- lncRNA_coordinates%>% dplyr::rename(start = transcription_start_site) %>% 
  dplyr::rename(end = TSS_1)

lncRNA_ranges <- makeGRangesFromDataFrame(lncRNA_coordinates,ignore.strand=T, 
                                          keep.extra.columns=T, seqinfo = seq_info)
lncRNA_ranges <- lncRNA_ranges + 100000

lncRNA_ranges <- trim(lncRNA_ranges)

head(lncRNA_ranges)

rm(lncRNA_coordinates, TSS_1)

###all TSS

TSS_coordinates <- gene_info %>% 
  dplyr::select(chromosome_name, ensembl_gene_id, external_gene_name, 
                ensembl_transcript_id, transcription_start_site) %>% 
  dplyr::rename(seqname = chromosome_name)

TSS_1 <- TSS_coordinates$transcription_start_site +1
TSS_coordinates <- TSS_coordinates %>% mutate(TSS_1 = TSS_1) 
TSS_coordinates <- TSS_coordinates%>% dplyr::rename(start = transcription_start_site) %>% 
  dplyr::rename(end = TSS_1)

TSS_ranges <- makeGRangesFromDataFrame(TSS_coordinates,ignore.strand=T, 
                                       keep.extra.columns=T, seqinfo = seq_info)
TSS_ranges <- TSS_ranges + 50000

TSS_ranges <- trim(TSS_ranges)

head(TSS_ranges)

rm(TSS_coordinates, TSS_1)


###gene coordinates
gene_coordinates <- gene_info %>% 
  dplyr::select(chromosome_name, ensembl_gene_id, external_gene_name, 
                start_position, end_position) %>% 
  dplyr::rename(seqname = chromosome_name, start = start_position, end = end_position)

gene_coordinates <- distinct(gene_coordinates)

gene_ranges <- makeGRangesFromDataFrame(gene_coordinates,ignore.strand=T,
                                        keep.extra.columns=T, seqinfo = seq_info)

head(gene_ranges)

rm(gene_coordinates, gene_info)


#### filtering safe coordinates

ranges_to_remove <- c(cancer_ranges, dnase_ranges, gene_ranges, lncRNA_ranges, miRNA_ranges,
                      TSS_ranges, ucr_ranges)
ranges_to_remove <- reduce(ranges_to_remove)


hg38_full_filter <- setdiff(hg38_ranges, ranges_to_remove)

write.xlsx(hg38_full_filter, "hg38_full_filter.xlsx")
