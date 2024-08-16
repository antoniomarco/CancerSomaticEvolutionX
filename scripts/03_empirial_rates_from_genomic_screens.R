# R
library(ggplot2)
library(patchwork)
library(GenomicFeatures)
library(ggbeeswarm)
library(scales)
# # # Build Homo.sapiens.hg38
# # Source of code, modified from
# # https://bioinformatics.stackexchange.com/questions/3950/how-to-get-results-from-homo-sapiens-package-in-bioconductor-for-a-specific-refe
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
# BiocManager::install("OrganismDbi")
# BiocManager::install("GO.db")
# BiocManager::install("org.Hs.eg.db")
# 
library(OrganismDbi)
# # instructions for hg19, for hg38 change accordingly
# gd <- list(join1 = c(GO.db="GOID", org.Hs.eg.db="GO"),
#              join2 = c(org.Hs.eg.db = "ENTREZID",
#                        TxDb.Hsapiens.UCSC.hg38.knownGene = "GENEID"))
# destination <- tempfile()
# dir.create("datasets")
# makeOrganismPackage(pkgname = "Homo.sapiens.hg38", graphData = gd,
#                       organism = "Homo sapiens", version = "1.0.0",
#                       maintainer = NA",
#                       author = "Author Name", destDir = "datasets",
#                       license = "Artistic-2.0")
# install.packages("datasets/Homo.sapiens.hg38", repos = NULL, type="source")
library(Homo.sapiens.hg38)


cosmic <- read.table(gzfile("datasets/Cosmic_GenomeScreen_parsed.tab.gz"), sep = "\t")

# Select uniquely 'missense_variant' (sometimes a change is various things at the same time)
cosmic_missense <- cosmic[cosmic$V3 == 'missense_variant',]


# Overlapping PUDMED COSstudies
#      51 15908952COSU34
#     125 16140923COSU22
#     128 16140923COSU34
#       3 16140923COSU36
#       5 16618716COSU34
#       2 16618716COSU38
#       3 16618716COSU7
#       3 20054297COSU255
#      14 20054297COSU679
#    7149 21642962COSU340
#   25190 21720365COSU331
#      82 21995386COSU351
#    2585 22608084COSU385
#    7986 22722201COSU385
#   15611 22810696COSU375
#   33453 22810696COSU376
#    1508 22832583COSU379
#    1585 23770606COSU486


# Read TSG/OG Census
census <- read.table(gzfile("datasets/Cosmic_CancerGeneCensus_v99_GRCh38.tsv.gz"), sep = "\t", header = TRUE, fill = TRUE)
## POG
POG_all <- census[grepl("oncogene", census$ROLE_IN_CANCER),]
# POG not TSG
POG_A <- POG_all[(!grepl("TSG", POG_all$ROLE_IN_CANCER)) & (POG_all$CHROMOSOME %in% 1:22),]$GENE_SYMBOL
POG_X <- POG_all[(!grepl("TSG", POG_all$ROLE_IN_CANCER)) & (POG_all$CHROMOSOME == "X"),]$GENE_SYMBOL
## TSG
TSG_all <- census[grepl("TSG", census$ROLE_IN_CANCER),]
# POG not TSG
TSG_A <- TSG_all[(!grepl("oncogene", TSG_all$ROLE_IN_CANCER)) & (TSG_all$CHROMOSOME %in% 1:22),]$GENE_SYMBOL
TSG_X <- TSG_all[(!grepl("oncogene", TSG_all$ROLE_IN_CANCER)) & (TSG_all$CHROMOSOME == "X"),]$GENE_SYMBOL
# Transcript legnths
txlens <- transcriptLengths(TxDb.Hsapiens.UCSC.hg38.knownGene, with.cds_len=TRUE)


# # For debugging
# my_test_list <- sample(unique(cosmic$V1), 100)


## Report by gene putting all studies altogether
# Function
merge_studies_report_genes <- function(my_test_list = my_test_list, input_dataset = cosmic_missense){
  cosmic_missense_subset <- input_dataset[input_dataset$V1 %in% my_test_list,]
  # Combined studies, for individual genes
  cosmic_missense_subset_splitGene <- split(cosmic_missense_subset$V4, cosmic_missense_subset$V1)
  cosmic_missense_subset_countGene <- data.frame(t(data.frame(lapply(cosmic_missense_subset_splitGene, sum))))
  enst2name <- data.frame(select(Homo.sapiens.hg38, keys = my_test_list, columns = c('SYMBOL', 'TXNAME'), keytype = 'SYMBOL'))
  enst2name <- enst2name[! is.na(enst2name$TXNAME),]
  rownames(enst2name) <- enst2name$TXNAME
  txlens_subset <- txlens[txlens$tx_name %in% rownames(enst2name),]
  txlens_subset$symbol <- enst2name[txlens_subset$tx_name,]$SYMBOL
  txlens_subset_split <- split(txlens_subset$cds_len, txlens_subset$symbol)
  txlens_subset_CDSlength <- t(data.frame(lapply(txlens_subset_split, max)))
  # And now combine (this is fo the gene only, do for the study!)
  ddd <- merge(cosmic_missense_subset_countGene, txlens_subset_CDSlength, by = "row.names")
  rownames(ddd) <- ddd[,1]
  ddd <- ddd[,c(2,3)]
  colnames(ddd) <- c("missense_mutations", "CDS_length")
  ddd$MS_NT <- ddd$missense_mutations / ddd$CDS_length
  return(ddd)
}

# analysis
POG_A_genes_rates <- merge_studies_report_genes(my_test_list = POG_A)
POG_X_genes_rates <- merge_studies_report_genes(my_test_list = POG_X)
TSG_A_genes_rates <- merge_studies_report_genes(my_test_list = TSG_A)
TSG_X_genes_rates <- merge_studies_report_genes(my_test_list = TSG_X)
POG_A_genes_rates$cat <- "POG_A"
POG_X_genes_rates$cat <- "POG_X"
TSG_A_genes_rates$cat <- "TSG_A"
TSG_X_genes_rates$cat <- "TSG_X"
Genes_Rates_POG <- rbind(POG_A_genes_rates, POG_X_genes_rates)
Genes_Rates_TSG <- rbind(TSG_A_genes_rates, TSG_X_genes_rates)
# change names:
Genes_Rates_TSG$cat <- gsub("TSG_", "", Genes_Rates_TSG$cat)
Genes_Rates_POG$cat <- gsub("POG_", "", Genes_Rates_POG$cat)

# POG
plots <- list()
plots[[1]] <- ggplot(data = Genes_Rates_POG, aes(x=cat, y=MS_NT)) + 
  geom_boxplot(outlier.shape = NA, col = "grey") +
  geom_beeswarm(size = 1) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  coord_flip() + 
  labs(x = "", y = "Mutations per nucleotide") +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  theme(plot.title = element_text(size=12)) +
  labs(title = "POG") +
  theme(plot.title = element_text(hjust=0.9, vjust = -50, margin=margin(t=40,b=-25)),
        plot.tag = element_text()) 


plots[[2]] <- ggplot(data = Genes_Rates_TSG, aes(x=cat, y=MS_NT)) + 
  geom_boxplot(outlier.shape = NA, col = "grey") +
  geom_beeswarm() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  coord_flip() + 
  labs(x = "", y = "Mutations per nucleotide") +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  theme(plot.title = element_text(size=14)) +
  labs(title = "TSG") +
  theme(plot.title = element_text(hjust=0.9, vjust = -50, margin=margin(t=40,b=-25)),
        plot.tag = element_text()) 



## Report for each study
# Function
merge_genes_report_studies <- function(my_test_list = my_test_list, input_dataset = cosmic_missense){
  cosmic_missense_subset <- input_dataset[input_dataset$V1 %in% my_test_list,]
  # Combined genes, for individual studies
  cosmic_missense_subset_splitStudy <- split(cosmic_missense_subset$V4, cosmic_missense_subset$V2)
  cosmic_missense_subset_countStudy <- data.frame(t(data.frame(lapply(cosmic_missense_subset_splitStudy, sum))))
  enst2name <- data.frame(select(Homo.sapiens.hg38, keys = my_test_list, columns = c('SYMBOL', 'TXNAME'), keytype = 'SYMBOL'))
  enst2name <- enst2name[! is.na(enst2name$TXNAME),]
  rownames(enst2name) <- enst2name$TXNAME
  # total length
  txlens_subset <- txlens[txlens$tx_name %in% rownames(enst2name),]
  txlens_subset$symbol <- enst2name[txlens_subset$tx_name,]$SYMBOL
  txlens_subset_split <- split(txlens_subset$cds_len, txlens_subset$symbol)
  txlens_subset_CDSlength <- t(data.frame(lapply(txlens_subset_split, max)))
  cosmic_missense_subset_countStudy$V2 <- sum(txlens_subset_CDSlength)
  cosmic_missense_subset_countStudy$V3 <- cosmic_missense_subset_countStudy[,1]/cosmic_missense_subset_countStudy[,2]
  colnames(cosmic_missense_subset_countStudy) <- c("MUT","LGTH", "MS_NT")
  return(cosmic_missense_subset_countStudy)
}

# Analysis
POG_A_studies_rates <- merge_genes_report_studies(POG_A)
POG_X_studies_rates <- merge_genes_report_studies(POG_X)
TSG_A_studies_rates <- merge_genes_report_studies(TSG_A)
TSG_X_studies_rates <- merge_genes_report_studies(TSG_X)
POG_A_studies_rates$cat <- "POG_A"
POG_X_studies_rates$cat <- "POG_X"
TSG_A_studies_rates$cat <- "TSG_A"
TSG_X_studies_rates$cat <- "TSG_X"
Genes_RatesS_POG <- rbind(POG_A_studies_rates, POG_X_studies_rates)
Genes_RatesS_TSG <- rbind(TSG_A_studies_rates, TSG_X_studies_rates)


# paired
# Unequal size!!!
POG_paired_studies_rates <- merge(POG_A_studies_rates, POG_X_studies_rates, by = 0)
TSG_paired_studies_rates <- merge(TSG_A_studies_rates, TSG_X_studies_rates, by = 0)


# PLOT
# combine
POG_paired_studies_rates$LOR <- log2(POG_paired_studies_rates$MS_NT.x/POG_paired_studies_rates$MS_NT.y)
POG_paired_studies_rates$TYPE <- "POG"
TSG_paired_studies_rates$LOR <- log2(TSG_paired_studies_rates$MS_NT.x/TSG_paired_studies_rates$MS_NT.y)
TSG_paired_studies_rates$TYPE <- "TSG"
Paired_studies <- rbind(POG_paired_studies_rates, TSG_paired_studies_rates)

plots[[3]] <- ggplot(Paired_studies, aes(LOR, fill = TYPE)) + geom_density(alpha = 0.2) + 
  labs(x = "log2-fold ratio (A/X)", y = "density") +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  theme(plot.title = element_text(size=12)) +
  # labs(title = "By study") +
  theme(plot.title = element_text(hjust=0.9, margin=margin(t=40,b=-30)),
        plot.tag = element_text())           


# OUTPUTS
writeLines(capture.output(wilcox.test(POG_paired_studies_rates$MS_NT.x, POG_paired_studies_rates$MS_NT.y, paired = TRUE)), con = file("results/POG_AvsX_studies_wilcox.txt"))
writeLines(capture.output(wilcox.test(TSG_paired_studies_rates$MS_NT.x, TSG_paired_studies_rates$MS_NT.y, paired = TRUE)), con = file("results/TSG_AvsX_studies_wilcox.txt"))
writeLines(capture.output(wilcox.test(POG_A_genes_rates$MS_NT, POG_X_genes_rates$MS_NT)), con = file("results/POG_AvsX_wilcox.txt"))
writeLines(capture.output(wilcox.test(TSG_A_genes_rates$MS_NT, TSG_X_genes_rates$MS_NT)), con = file("results/TSG_AvsX_wilcox.txt"))
# plot
ggsave(plots[[1]] / plots[[2]] / plots[[3]] + geom_density(alpha = 0.2) + plot_annotation(tag_levels = "A"), width = 9, units = "cm", file = "plots/Figure_5.png")



## Split by sex
## by sex

# Load samples
fsamples <- read.table(gzfile("datasets/female_samples.txt.gz"))$V1
msamples <- read.table(gzfile("datasets/male_samples.txt.gz"))$V1

# Cosmic filter
cosmic_missense_female <- cosmic_missense[cosmic_missense$V5 %in% fsamples,]
cosmic_missense_male <- cosmic_missense[cosmic_missense$V5 %in% msamples,]

# Female by study
POGf_A_studies_rates <- merge_genes_report_studies(my_test_list = POG_A, input_dataset = cosmic_missense_female)
POGf_X_studies_rates <- merge_genes_report_studies(my_test_list = POG_X, input_dataset = cosmic_missense_female)
TSGf_A_studies_rates <- merge_genes_report_studies(my_test_list = TSG_A, input_dataset = cosmic_missense_female)
TSGf_X_studies_rates <- merge_genes_report_studies(my_test_list = TSG_X, input_dataset = cosmic_missense_female)
POGf_paired_studies_rates <- merge(POGf_A_studies_rates, POGf_X_studies_rates, by = 0)
TSGf_paired_studies_rates <- merge(TSGf_A_studies_rates, TSGf_X_studies_rates, by = 0)
POGf_paired_studies_rates$LOR <- log2(POGf_paired_studies_rates$MS_NT.x/POGf_paired_studies_rates$MS_NT.y)
POGf_paired_studies_rates$TYPE <- "POG"
TSGf_paired_studies_rates$LOR <- log2(TSGf_paired_studies_rates$MS_NT.x/TSGf_paired_studies_rates$MS_NT.y)
TSGf_paired_studies_rates$TYPE <- "TSG"
Paired_studies_female <- rbind(POGf_paired_studies_rates, TSGf_paired_studies_rates)


# Male by study
POGm_A_studies_rates <- merge_genes_report_studies(my_test_list = POG_A, input_dataset = cosmic_missense_male)
POGm_X_studies_rates <- merge_genes_report_studies(my_test_list = POG_X, input_dataset = cosmic_missense_male)
TSGm_A_studies_rates <- merge_genes_report_studies(my_test_list = TSG_A, input_dataset = cosmic_missense_male)
TSGm_X_studies_rates <- merge_genes_report_studies(my_test_list = TSG_X, input_dataset = cosmic_missense_male)
POGm_paired_studies_rates <- merge(POGm_A_studies_rates, POGm_X_studies_rates, by = 0)
TSGm_paired_studies_rates <- merge(TSGm_A_studies_rates, TSGm_X_studies_rates, by = 0)
POGm_paired_studies_rates$LOR <- log2(POGm_paired_studies_rates$MS_NT.x/POGm_paired_studies_rates$MS_NT.y)
POGm_paired_studies_rates$TYPE <- "POG"
TSGm_paired_studies_rates$LOR <- log2(TSGm_paired_studies_rates$MS_NT.x/TSGm_paired_studies_rates$MS_NT.y)
TSGm_paired_studies_rates$TYPE <- "TSG"
Paired_studies_male <- rbind(POGm_paired_studies_rates, TSGm_paired_studies_rates)



plot_hist_1 <- ggplot(Paired_studies_male, aes(LOR, fill = TYPE)) + geom_density(alpha = 0.2) + 
  labs(x = "log2-fold ratio (A/X)", y = "density") +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  theme(plot.title = element_text(size=20)) +
  labs(title = "\U2642") +
  theme(plot.title = element_text(hjust=0, margin=margin(t=40, b=-80)),
        plot.tag = element_text())           

plot_hist_2 <- ggplot(Paired_studies_female, aes(LOR, fill = TYPE)) + geom_density(alpha = 0.2) + 
  labs(x = "log2-fold ratio (A/X)", y = "density") +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  theme(plot.title = element_text(size=20)) +
  labs(title = "\U2640") +
  theme(plot.title = element_text(hjust=0, margin=margin(t=40, b=-80)),
        plot.tag = element_text())           

# OUTPUTS
writeLines(capture.output(wilcox.test(Paired_studies_male[Paired_studies_male$TYPE == "POG",]$MS_NT.x, Paired_studies_male[Paired_studies_male$TYPE == "POG",]$MS_NT.y, paired = TRUE)), con = file("results/POG_AvsX_studies_male_wilcox.txt"))
writeLines(capture.output(wilcox.test(Paired_studies_male[Paired_studies_male$TYPE == "TSG",]$MS_NT.x, Paired_studies_male[Paired_studies_male$TYPE == "TSG",]$MS_NT.y, paired = TRUE)), con = file("results/TSG_AvsX_studies_male_wilcox.txt"))
writeLines(capture.output(wilcox.test(Paired_studies_female[Paired_studies_female$TYPE == "POG",]$MS_NT.x, Paired_studies_female[Paired_studies_female$TYPE == "POG",]$MS_NT.y, paired = TRUE)), con = file("results/POG_AvsX_studies_female_wilcox.txt"))
writeLines(capture.output(wilcox.test(Paired_studies_female[Paired_studies_female$TYPE == "TSG",]$MS_NT.x, Paired_studies_female[Paired_studies_female$TYPE == "TSG",]$MS_NT.y, paired = TRUE)), con = file("results/TSG_AvsX_studies_female_wilcox.txt"))


## STOP GAINED

# Select uniquely 'stop_gained'
cosmic_stop <- cosmic[cosmic$V3 == 'stop_gained',]

# Cosmic filter
cosmic_stop_female <- cosmic_stop[cosmic_stop$V5 %in% fsamples,]
cosmic_stop_male <- cosmic_stop[cosmic_stop$V5 %in% msamples,]

# Female by study
POGf_A_studies_rates_stop <- merge_genes_report_studies(my_test_list = POG_A, input_dataset = cosmic_stop_female)
POGf_X_studies_rates_stop <- merge_genes_report_studies(my_test_list = POG_X, input_dataset = cosmic_stop_female)
TSGf_A_studies_rates_stop <- merge_genes_report_studies(my_test_list = TSG_A, input_dataset = cosmic_stop_female)
TSGf_X_studies_rates_stop <- merge_genes_report_studies(my_test_list = TSG_X, input_dataset = cosmic_stop_female)
POGf_paired_studies_rates_stop <- merge(POGf_A_studies_rates_stop, POGf_X_studies_rates_stop, by = 0)
TSGf_paired_studies_rates_stop <- merge(TSGf_A_studies_rates_stop, TSGf_X_studies_rates_stop, by = 0)
POGf_paired_studies_rates_stop$LOR <- log2(POGf_paired_studies_rates_stop$MS_NT.x/POGf_paired_studies_rates_stop$MS_NT.y)
POGf_paired_studies_rates_stop$TYPE <- "POG"
TSGf_paired_studies_rates_stop$LOR <- log2(TSGf_paired_studies_rates_stop$MS_NT.x/TSGf_paired_studies_rates_stop$MS_NT.y)
TSGf_paired_studies_rates_stop$TYPE <- "TSG"
Paired_studies_female_stop <- rbind(POGf_paired_studies_rates_stop, TSGf_paired_studies_rates_stop)


# Male by study
POGm_A_studies_rates_stop <- merge_genes_report_studies(my_test_list = POG_A, input_dataset = cosmic_stop_male)
POGm_X_studies_rates_stop <- merge_genes_report_studies(my_test_list = POG_X, input_dataset = cosmic_stop_male)
TSGm_A_studies_rates_stop <- merge_genes_report_studies(my_test_list = TSG_A, input_dataset = cosmic_stop_male)
TSGm_X_studies_rates_stop <- merge_genes_report_studies(my_test_list = TSG_X, input_dataset = cosmic_stop_male)
POGm_paired_studies_rates_stop <- merge(POGm_A_studies_rates_stop, POGm_X_studies_rates_stop, by = 0)
TSGm_paired_studies_rates_stop <- merge(TSGm_A_studies_rates_stop, TSGm_X_studies_rates_stop, by = 0)
POGm_paired_studies_rates_stop$LOR <- log2(POGm_paired_studies_rates_stop$MS_NT.x/POGm_paired_studies_rates_stop$MS_NT.y)
POGm_paired_studies_rates_stop$TYPE <- "POG"
TSGm_paired_studies_rates_stop$LOR <- log2(TSGm_paired_studies_rates_stop$MS_NT.x/TSGm_paired_studies_rates_stop$MS_NT.y)
TSGm_paired_studies_rates_stop$TYPE <- "TSG"
Paired_studies_male_stop <- rbind(POGm_paired_studies_rates_stop, TSGm_paired_studies_rates_stop)

plot_hist_1_stop <- ggplot(Paired_studies_male_stop, aes(LOR, fill = TYPE)) + geom_density(alpha = 0.2) + 
  labs(x = "log2-fold ratio (A/X)", y = "density") +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  theme(plot.title = element_text(size=20)) +
  labs(title = "\U2642") +
  theme(plot.title = element_text(hjust=0, margin=margin(t=40, b=-80)),
        plot.tag = element_text())           

plot_hist_2_stop <- ggplot(Paired_studies_female_stop, aes(LOR, fill = TYPE)) + geom_density(alpha = 0.2) + 
  labs(x = "log2-fold ratio (A/X)", y = "density") +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  theme(plot.title = element_text(size=20)) +
  labs(title = "\U2640") +
  theme(plot.title = element_text(hjust=0, margin=margin(t=40, b=-80)),
        plot.tag = element_text())           



# Select uniquely 'synonymous_variant'
cosmic_syn <- cosmic[cosmic$V3 == 'synonymous_variant',]

# Cosmic filter
cosmic_syn_female <- cosmic_syn[cosmic_syn$V5 %in% fsamples,]
cosmic_syn_male <- cosmic_syn[cosmic_syn$V5 %in% msamples,]

# Female by study
POGf_A_studies_rates_syn <- merge_genes_report_studies(my_test_list = POG_A, input_dataset = cosmic_syn_female)
POGf_X_studies_rates_syn <- merge_genes_report_studies(my_test_list = POG_X, input_dataset = cosmic_syn_female)
TSGf_A_studies_rates_syn <- merge_genes_report_studies(my_test_list = TSG_A, input_dataset = cosmic_syn_female)
TSGf_X_studies_rates_syn <- merge_genes_report_studies(my_test_list = TSG_X, input_dataset = cosmic_syn_female)
POGf_paired_studies_rates_syn <- merge(POGf_A_studies_rates_syn, POGf_X_studies_rates_syn, by = 0)
TSGf_paired_studies_rates_syn <- merge(TSGf_A_studies_rates_syn, TSGf_X_studies_rates_syn, by = 0)
POGf_paired_studies_rates_syn$LOR <- log2(POGf_paired_studies_rates_syn$MS_NT.x/POGf_paired_studies_rates_syn$MS_NT.y)
POGf_paired_studies_rates_syn$TYPE <- "POG"
TSGf_paired_studies_rates_syn$LOR <- log2(TSGf_paired_studies_rates_syn$MS_NT.x/TSGf_paired_studies_rates_syn$MS_NT.y)
TSGf_paired_studies_rates_syn$TYPE <- "TSG"
Paired_studies_female_syn <- rbind(POGf_paired_studies_rates_syn, TSGf_paired_studies_rates_syn)


# Male by study
POGm_A_studies_rates_syn <- merge_genes_report_studies(my_test_list = POG_A, input_dataset = cosmic_syn_male)
POGm_X_studies_rates_syn <- merge_genes_report_studies(my_test_list = POG_X, input_dataset = cosmic_syn_male)
TSGm_A_studies_rates_syn <- merge_genes_report_studies(my_test_list = TSG_A, input_dataset = cosmic_syn_male)
TSGm_X_studies_rates_syn <- merge_genes_report_studies(my_test_list = TSG_X, input_dataset = cosmic_syn_male)
POGm_paired_studies_rates_syn <- merge(POGm_A_studies_rates_syn, POGm_X_studies_rates_syn, by = 0)
TSGm_paired_studies_rates_syn <- merge(TSGm_A_studies_rates_syn, TSGm_X_studies_rates_syn, by = 0)
POGm_paired_studies_rates_syn$LOR <- log2(POGm_paired_studies_rates_syn$MS_NT.x/POGm_paired_studies_rates_syn$MS_NT.y)
POGm_paired_studies_rates_syn$TYPE <- "POG"
TSGm_paired_studies_rates_syn$LOR <- log2(TSGm_paired_studies_rates_syn$MS_NT.x/TSGm_paired_studies_rates_syn$MS_NT.y)
TSGm_paired_studies_rates_syn$TYPE <- "TSG"
Paired_studies_male_syn <- rbind(POGm_paired_studies_rates_syn, TSGm_paired_studies_rates_syn)

plot_hist_1_syn <- ggplot(Paired_studies_male_syn, aes(LOR, fill = TYPE)) + geom_density(alpha = 0.2) + 
  labs(x = "log2-fold ratio (A/X)", y = "density") +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  theme(plot.title = element_text(size=20)) +
  labs(title = "\U2642") +
  theme(plot.title = element_text(hjust=0, margin=margin(t=40, b=-80)),
        plot.tag = element_text())           

plot_hist_2_syn <- ggplot(Paired_studies_female_syn, aes(LOR, fill = TYPE)) + geom_density(alpha = 0.2) + 
  labs(x = "log2-fold ratio (A/X)", y = "density") +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  theme(plot.title = element_text(size=20)) +
  labs(title = "\U2640") +
  theme(plot.title = element_text(hjust=0, margin=margin(t=40, b=-80)),
        plot.tag = element_text())           


# plot
ggsave((plot_hist_1 + plot_hist_2) / (plot_hist_1_stop + plot_hist_2_stop) / (plot_hist_1_syn + plot_hist_2_syn) + plot_annotation(tag_levels = "A") + plot_layout(guides = "collect"), width = 18, units = "cm", file = "plots/Figure_6.png")



# BREAST
# 5 genome screens in COSMIC
genomics_screens_breast <- c("COSU652", "COSU414", "22608084COSU385", "22722201COSU385", "COSU385", "COSU669", "COSU668")
cosmic_missense_breast <- cosmic_missense[which(cosmic_missense$V2 %in% genomics_screens_breast),]
POG_A_breast <- merge_studies_report_genes(my_test_list = POG_A, input_dataset = cosmic_missense_breast)
POG_X_breast <- merge_studies_report_genes(my_test_list = POG_X, input_dataset = cosmic_missense_breast)
TSG_A_breast <- merge_studies_report_genes(my_test_list = TSG_A, input_dataset = cosmic_missense_breast)
TSG_X_breast <- merge_studies_report_genes(my_test_list = TSG_X, input_dataset = cosmic_missense_breast)
POG_A_breast$cat <- "A"
POG_X_breast$cat <- "X"
TSG_A_breast$cat <- "A"
TSG_X_breast$cat <- "X"
Genes_Rates_breast_POG <- rbind(POG_A_breast, POG_X_breast)
Genes_Rates_breast_TSG <- rbind(TSG_A_breast, TSG_X_breast)

# POG
plots <- list()
plots[[1]] <- ggplot(data = Genes_Rates_breast_POG, aes(x=cat, y=MS_NT)) + 
  geom_boxplot(outlier.shape = NA, col = "grey") +
  geom_beeswarm(size = 1) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  coord_flip() + 
  labs(x = "", y = "Mutations per nucleotide") +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  theme(plot.title = element_text(size=14)) +
  labs(title = "OG-BC") +
  theme(plot.title = element_text(hjust=0.9, margin=margin(t=40,b=-25)),
        plot.tag = element_text()) 
plots[[2]] <- ggplot(data = Genes_Rates_breast_TSG, aes(x=cat, y=MS_NT, label=rownames(Genes_Rates_breast_TSG))) + 
  geom_boxplot(outlier.shape = NA, col = "grey") +
  geom_beeswarm(size = 1) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_text(aes(label=ifelse((MS_NT>0.1) & (cat == "X"),rownames(Genes_Rates_breast_TSG),'')),hjust=0.5, vjust=2) +
  coord_flip() + 
  labs(x = "", y = "Mutations per nucleotide") +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  theme(plot.title = element_text(size=14)) +
  labs(title = "TSG-BC") +
  theme(plot.title = element_text(hjust=0.9, margin=margin(t=40,b=-25)),
        plot.tag = element_text())           

# PROSTATE
# 5 genome screens in COSMIC
genomics_screens_prostate <- c("COSU435", "COSU538", "COSU675", "COSU537", "COSU534")
cosmic_missense_prostate <- cosmic_missense[which(cosmic_missense$V2 %in% genomics_screens_prostate),]
POG_A_prostate <- merge_studies_report_genes(my_test_list = POG_A, input_dataset = cosmic_missense_prostate)
POG_X_prostate <- merge_studies_report_genes(my_test_list = POG_X, input_dataset = cosmic_missense_prostate)
TSG_A_prostate <- merge_studies_report_genes(my_test_list = TSG_A, input_dataset = cosmic_missense_prostate)
TSG_X_prostate <- merge_studies_report_genes(my_test_list = TSG_X, input_dataset = cosmic_missense_prostate)
POG_A_prostate$cat <- "A"
POG_X_prostate$cat <- "X"
TSG_A_prostate$cat <- "A"
TSG_X_prostate$cat <- "X"
Genes_Rates_prostate_POG <- rbind(POG_A_prostate, POG_X_prostate)
Genes_Rates_prostate_TSG <- rbind(TSG_A_prostate, TSG_X_prostate)

plots[[3]] <- ggplot(data = Genes_Rates_prostate_POG, aes(x=cat, y=MS_NT)) + 
  geom_boxplot(outlier.shape = NA, col = "grey") +
  geom_beeswarm(size = 1) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  coord_flip() + 
  labs(x = "", y = "Mutations per nucleotide") +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  theme(plot.title = element_text(size=14)) +
  labs(title = "OG-PCa") +
  theme(plot.title = element_text(hjust=0.9, margin=margin(t=40,b=-25)),
        plot.tag = element_text()) 
plots[[4]] <- ggplot(data = Genes_Rates_prostate_TSG, aes(x=cat, y=MS_NT, label=rownames(Genes_Rates_prostate_TSG))) + 
  geom_boxplot(outlier.shape = NA, col = "grey") +
  geom_beeswarm(size = 1) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_text(aes(label=ifelse((MS_NT>0.1) & (cat == "X"),rownames(Genes_Rates_prostate_TSG),'')),hjust=0.5, vjust=2) +
  coord_flip() + 
  labs(x = "", y = "Mutations per nucleotide") +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  theme(plot.title = element_text(size=14)) +
  labs(title = "TSG-PCa") +
  theme(plot.title = element_text(hjust=0.9, margin=margin(t=40,b=-25)),
        plot.tag = element_text()) 


# OUTPUTS
writeLines(capture.output(wilcox.test(POG_A_breast$MS_NT, POG_X_breast$MS_NT)), con = file("results/POG_AvsX_studies_breast_wilcox.txt"))
writeLines(capture.output(wilcox.test(TSG_A_breast$MS_NT, TSG_X_breast$MS_NT)), con = file("results/TSG_AvsX_studies_breast_wilcox.txt"))
writeLines(capture.output(wilcox.test(POG_A_prostate$MS_NT, POG_X_prostate$MS_NT)), con = file("results/POG_AvsX_studies_prostate_wilcox.txt"))
writeLines(capture.output(wilcox.test(TSG_A_prostate$MS_NT, TSG_X_prostate$MS_NT)), con = file("results/TSG_AvsX_studies_prostate_wilcox.txt"))
# plot
ggsave(((plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]]) ) + plot_annotation(tag_levels = "A"), width = 18, units = "cm", file = "plots/Figure_7.png")

#                missense_mutations CDS_length       MS_NT cat
# AMER1                 140       3408 0.041079812   X
# ATP2B3                 53       3663 0.014469014   X
# ATRX                  208       7479 0.027811205   X
# BCOR                  173       5268 0.032839787   X
# DDX3X                 333       2202 0.151226158   X
# KDM5C                  97       4683 0.020713218   X
# MED12                 106       6534 0.016222834   X
# ZMYM3                  89       4119 0.021607186   X
# ZRSR2                   6       1449 0.004140787   X

# Still DDX3X fast evolving gene        

# exit
q()
