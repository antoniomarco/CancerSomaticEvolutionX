# R
library(tidyverse)
library(patchwork)
library(scales)


cosmic <- read.table(gzfile("datasets/Cosmic_GenomeScreen_parsed.tab.gz"), sep = "\t")

# Nei-Gojobori syn/nonsyn sites
NG86_sites <- read.table(gzfile("datasets/GENE_SN_NG86.tab.gz"), sep = "\t", header = TRUE, fill = TRUE)

# Select uniquely 'missense_variant' (sometimes a change is various things at the same time)
cosmic_missense <- cosmic[cosmic$V3 == 'missense_variant',]
cosmic_synonymous <- cosmic[cosmic$V3 == 'synonymous_variant',]
cosmic_nonsense <- cosmic[cosmic$V3 == 'stop_gained',]
# It is reasonable to use the number of codons as a proxy of the number of STOPpable mutations per sequence:
# 18 codons can mutate into a STOP codon via 1 of 9 possible single nucleotide substitutions.
# 3 codons can mutate into a STOP codon via 2 such substitutions.
# 1 codon (possibly UAC) can mutate into a STOP codon via 3 different single substitutions.


# Number of missense, nonsense and synonymous variants per gene
cosmic_all_syn <- cosmic_synonymous[,c(1,4)] |>
  group_by(V1) |>
  summarise(syn = sum(V4))
cosmic_all_nonsyn <- cosmic_missense[,c(1,4)] |>
  group_by(V1) |>
  summarise(nonsyn = sum(V4))
cosmic_all_nonsense <- cosmic_nonsense[,c(1,4)] |>
  group_by(V1) |>
  summarise(nonsense = sum(V4))
cosmic_all_NS <- cosmic_all_syn |>
  left_join(cosmic_all_nonsyn, by = "V1") |>
  left_join(cosmic_all_nonsense, by = "V1") |>
  left_join(NG86_sites, join_by(V1 == name)) |>
  rename(gene = V1, syn = syn.x, nonsyn = nonsyn.x, synSites = syn.y, nonsynSites = nonsyn.y) |>
  mutate(LOR = log((nonsyn/nonsynSites)/(syn/synSites))) |>
  mutate(SD = sqrt((1/nonsyn)+(1/nonsynSites)+(1/syn)+(1/synSites))) |>
  mutate(zscore = LOR/SD) |>
  mutate(p = pnorm(abs(zscore), mean = 0, sd = 1, lower.tail = FALSE)*2)
cosmic_all_NS$q <- p.adjust(cosmic_all_NS$p, method = "fdr")


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

# ESCAPERS, from table Suppl 1, Slavney et al. MBE 33:384, 2015
XCI_escapers <- c("ABCB7", "ABCD1", "ACOT9", "ACRC", "AKAP17A", "ALG13",
                  "AMMECR1", "AP1S2", "ARHGAP4", "ARHGAP6", "ARHGEF9", 
                  "ARMCX1", "ARMCX2", "ARMCX4", "ARSD", "ARSE", "ARX", 
                  "ASB11", "ASMTL", "ATG4A", "ATP11C", "ATP2B3", "ATP6AP1", 
                  "ATP6AP2", "ATP7A", "AVPR2", "BCAP31", "BCOR", "BEX1",
                  "BRCC3", "BRS3", "BTK", "C1GALT1C1", "CA5B", "CA5BP1", 
                  "CCDC120", "CCNB3", "CD40LG", "CD99", "CD99L2", "CD99P1",
                  "CDK16", "CDKL5", "CDR1", "CHM", "CHRDL1", "CITED1", 
                  "CLCN4", "CLDN2", "CLIC2", "COL4A6", "COX7B", "CTPS2",
                  "CUL4B", "CXORF21", "CXORF23", "CXORF38", "CXORF56", 
                  "DDX3X", "DHRSX", "DIAPH2", "DLG3", "DMD", "DMRTC1", 
                  "DOCK11", "DRP2", "DUSP21", "EBP", "EDA2R", "EFNB1", 
                  "EIF1AX", "EIF2S3", "ENOX2", "ERCC6L", "FAAH2", "FAM104B", 
                  "FAM122B", "FAM50A", "FAM58A", "FAM9C", "FGD1", "FHL1",
                  "FLNA", "FOXO4", "FRAXA", "FRMPD3", "FRMPD4", "FTSJ1", 
                  "FUNDC1", "GAB3", "GATA1", "GDI1", "GEMIN8", "GLA", "GLRA2",
                  "GPC4", "GPKOW", "GPM6B", "GPR143", "GPR173", "GPR50",
                  "GRIA3", "GRIPAP1", "GRPR", "GTPBP6", "GYG2", "HAUS7", 
                  "HCFC1", "HDHD1", "HMGN5", "HPRT1", "HSD17B10", "HTR2C",
                  "HUWE1", "IDH3G", "IKBKG", "IL9R", "INE1", "IQSEC2", 
                  "IRAK1", "ITGB1BP2", "ITM2A", "KAL1", "KCND1", "KDM5C", 
                  "KDM6A", "KLHL13", "L1CAM", "LANCL3", "LAS1L", "LDOC1",
                  "LHFPL1", "LOC389895", "LOC92249", "LPAR4", "MAGEA8", 
                  "MAGED4", "MAGEH1", "MAGIX", "MAOA", "MAP7D2", "MCTS1", 
                  "MECP2", "MED14", "MID1", "MKRN4P", "MORF4L2", "MOSPD1", 
                  "MPC1L", "MSL3", "MSN", "MUM1L1", "MXRA5", "MYCLP1", 
                  "NAA10", "NAP1L3", "NDUFA1", "NHS", "NKAP", "NKAPP1", 
                  "NKRF", "NLGN3", "NLGN4X", "NXF2B", "NXF3", "NXT2",
                  "OFD1", "OTUD5", "PABPC5", "PAK3", "PCDH11X", "PDZD11",
                  "PDZD4", "PGRMC1", "PHEX", "PHF6", "PHKA2", "PIM2", "PIN4",
                  "PIR", "PJA1", "PLAC1", "PLCXD1", "PLP1", "PLS3", "PLXNB3",
                  "PNMA5", "PNPLA4", "PPP2R3B", "PRKX", "PRPS1", "PSMD10", 
                  "RAB40A", "RAB9A", "RAB9B", "RBBP7", "RBM10", "RBM3", "RBM41",
                  "RENBP", "REPS2", "RGN", "RHOXF1", "RHOXF2", "RIBC1", "RPA4",
                  "RPL36A", "RPS4X", "RS1", "S100G", "SASH3", "SAT1", "SH3BGRL",
                  "SH3KBP1", "SHOX", "SHROOM2", "SLC25A43", "SLC25A6", 
                  "SLC35A2", "SLC38A5", "SLC6A8", "SLC9A7", "SLITRK2", "SMC1A",
                  "SPANXN1", "SPIN3", "SRPX2", "SSR4", "STARD8", "STS", 
                  "SUV39H1", "SYAP1", "SYP", "SYTL4", "TAB3", "TAF9B", "TAZ",
                  "TBC1D25", "TBL1X", "TCEAL4", "TCEAL7", "TCEANC", "TENM1",
                  "TFE3", "TIMM8A", "TIMP1", "TKTL1", "TMEM164", "TMEM187",
                  "TMEM27", "TMSB15A", "TRAPPC2", "TREX2", "TRPC5", "TSPAN6",
                  "TSPAN7", "TSR2", "TXLNG", "UBA1", "UBE2A", "UBQLN2",
                  "UPF3B", "USP11", "USP9X", "UTP14A", "VAMP7", "VBP1", "WAS",
                  "WBP5", "WWC3", "XG", "XIST", "XKRX", "XPNPEP2", "ZBED1", 
                  "ZC4H2", "ZFX", "ZNF185", "ZRSR2")
XCI_inactivated <- c("SEPT6", "ACE2", "ACSL4", "ACTRT1", "AFF2", "AGTR2",
                     "AIFM1", "AKAP14", "AKAP4", "ALAS2", "AMELX", "AMER1", 
                     "AMOT", "APEX2", "APLN", "APOO", "APOOL", "AR", "ARAF",
                     "ARHGAP36", "ARHGEF6", "ARMCX3", "ARMCX5", "ARMCX6", 
                     "ARR3", "ARSF", "ARSH", "ASB12", "ASB9", "ASMT", "ATP1B4",
                     "ATRX", "AWAT1", "BCORL1", "BEND2", "BEX2", "BEX4", 
                     "BEX5", "BGN", "BHLHB9", "BMP15", "BMX", "BRWD3", 
                     "CACNA1F", "CAPN6", "CASK", "CCDC22", "CDX4", "CENPI", 
                     "CETN2", "CFP", "CHDC2", "CHIC1", "CHST7", "CLCN5", 
                     "CMC4", "CNGA2", "CNKSR2", "COL4A5", "CPXCR1", "CSAG1", 
                     "CSAG2", "CSF2RA", "CSTF2", "CT45A2", "CT45A3", "CT45A4",
                     "CT45A5", "CT45A6", "CT47A1", "CT47A10", "CT47A11", 
                     "CT47A2", "CT47A3", "CT47A4", "CT47A5", "CT47A6", 
                     "CT47A7", "CT47A8", "CT47A9", "CT55", "CT83", "CTAG1A",
                     "CTAG1B", "CTAG2", "CXCR3", "CXORF22", "CXORF27", 
                     "CXORF36", "CXORF40A", "CXORF40B", "CXORF57", "CXORF58", 
                     "CXORF65", "CXORF66", "CXORF67", "CYBB", "CYLC1", 
                     "CYSLTR1", "DACH2", "DCAF12L1", "DCAF12L2", "DCAF8L1", 
                     "DCX", "DDX26B", "DDX53", "DGAT2L6", "DKC1", "DMRTC1B",
                     "DNASE1L1", "DUSP9", "DYNLT3", "EDA", "EGFL6", "ELF4", 
                     "ELK1", "EMD", "ERAS", "ESX1", "F8", "F8A1", "F9", 
                     "FAM120C", "FAM122C", "FAM127A", "FAM127B", "FAM127C", 
                     "FAM133A", "FAM155B", "FAM156A", "FAM199X", "FAM3A",
                     "FAM46D", "FAM47A", "FAM47B", "FAM47C", "FAM9A", "FAM9B",
                     "FANCB", "FATE1", "FGF13", "FIGF", "FMR1NB", "FOXP3",
                     "FRMD7", "FTHL17", "FUNDC2", "GABRA3", "GABRE", "GABRQ", 
                     "GAGE1", "GAGE10", "GAGE12C", "GAGE12E", "GAGE12F", 
                     "GAGE12G", "GAGE12H", "GAGE12J", "GAGE2C", "GAGE2D", 
                     "GDPD2", "GJB1", "GK", "GLRA4", "GLUD2", "GNL3L", 
                     "GPR101", "GPR112", "GPR119", "GPR174", "GPR34", "GPR64", 
                     "GPR82", "GPRASP1", "GPRASP2", "GSPT2", "GUCY2F", 
                     "H2AFB1", "H2AFB2", "H2AFB3", "H2BFWT", "HCCS", "HDAC6",
                     "HDAC8", "HDX", "HEPH", "HMGB3", "HNRNPH2", "HSFX1", 
                     "HTATSF1", "IDS", "IGBP1", "IGSF1", "IL13RA1", "IL13RA2", 
                     "IL1RAPL1", "IL1RAPL2", "IL2RG", "IL3RA", "IRS4", "ITIH6",
                     "JADE3", "KCNE1L", "KIAA2022", "KIF4A", "KLF8", "KLHL15", 
                     "KLHL34", "KLHL4", "KRBOX4", "LAGE3", "LAMP2", 
                     "LINC00894", "LOC105373383", "LONRF3", "LRCH2", "LUZP4", 
                     "MAGEA1", "MAGEA10", "MAGEA12", "MAGEA2B", "MAGEA3",
                     "MAGEA5", "MAGEA6", "MAGEA9", "MAGEA9B", "MAGEB1", 
                     "MAGEB10", "MAGEB16", "MAGEB18", "MAGEB2", "MAGEB3", 
                     "MAGEB4", "MAGEB6", "MAGEC1", "MAGEC2", "MAGEC3", "MAGED1",
                     "MAGED2", "MAGED4B", "MAGEE2", "MAGT1", "MAMLD1", "MAOB",
                     "MAP3K15", "MAP7D3", "MBNL3", "MBTPS2", "MCF2", "MED12", 
                     "MID1IP1", "MID2", "MIR503HG", "MMGT1", "MORC4", "MOSPD2", 
                     "MPP1", "MST4", "MTM1", "MTMR1", "MTMR8", "NAP1L2", "NDP", 
                     "NDUFB11", "NGFRAP1", "NHSL2", "NONO", "NOX1", "NR0B1",
                     "NSDHL", "NUDT10", "NUDT11", "NUP62CL", "NXF2", "NXF5", 
                     "NYX", "OCRL", "OGT", "OPHN1", "OPN1LW", "OPN1MW",
                     "OPN1MW2", "OR13H1", "OTC", "OTUD6A", "P2RY10", "P2RY4",
                     "P2RY8", "PABPC1L2A", "PABPC1L2B", "PAGE1", "PAGE2", 
                     "PAGE2B", "PAGE4", "PAGE5", "PASD1", "PBDC1", "PCDH19",
                     "PCSK1N", "PCYT1B", "PDHA1", "PDK3", "PFKFB1", "PGAM4", 
                     "PGK1", "PHF8", "PHKA1", "PIGA", "PIH1D3", "PLP2",
                     "PLXNA3", "PNCK", "PNMA3", "PNMA6A", "POF1B", "POLA1",
                     "PORCN", "POU3F4", "PPEF1", "PPP1R3F", "PQBP1", "PRAF2",
                     "PRDX4", "PRICKLE3", "PRPS2", "PRRG1", "PRRG3", "PTCHD1", 
                     "RAB33A", "RAB39B", "RAB40AL", "RAB41", "RAI2", "RAP2C", 
                     "RBMX", "RBMX2", "RGAG4", "RHOXF2B", "RIPPLY1", "RLIM",
                     "RNF113A", "RNF128", "RP2", "RPGR", "RPL10", "RPL39", 
                     "RPS6KA3", "RPS6KA6", "RRAGB", "SAGE1", "SATL1", "SCML1",
                     "SCML2", "SH2D1A", "SHROOM4", "SLC10A3", "SLC16A2", 
                     "SLC25A14", "SLC25A5", "SLC25A53", "SLC6A14", "SLC7A3", 
                     "SLC9A6", "SLITRK4", "SMARCA1", "SMIM10", "SMPX", "SMS",
                     "SNX12", "SOWAHD", "SOX3", "SPACA5", "SPACA5B", "SPANXA1",
                     "SPANXA2", "SPANXB1", "SPANXC", "SPANXD", "SPANXN2",
                     "SPANXN3", "SPANXN5", "SPIN2A", "SPIN2B", "SPIN4", "SPRY3",
                     "SRPK3", "SRPX", "SSX1", "SSX2", "SSX3", "SSX4", "SSX4B", 
                     "SSX5", "SSX6", "STAG2", "SYN1", "SYTL5", "TAF1", "TAF7L",
                     "TBC1D8B", "TBX22", "TCEAL1", "TCEAL2", "TCEAL5", 
                     "TCEAL6", "TCEAL8", "TEX11", "TEX13B", "TEX28", "TGIF2LX",
                     "THOC2", "TIMM17B", "TLR7", "TLR8", "TMEM185A", 
                     "TMEM255A", "TMEM257", "TMEM31", "TMEM35", "TMEM47",
                     "TMLHE", "TMSB15B", "TMSB4X", "TNMD", "TRMT2B", "TRO", 
                     "TSC22D3", "TSPYL2", "UBE2NL", "UBL4A", "UPRT", "USP26", 
                     "USP27X", "UXT", "VCX", "VCX2", "VCX3A", "VGLL1", "VMA21",
                     "VSIG1", "VSIG4", "WDR13", "WDR44", "WDR45", "WNK3",
                     "XAGE1C", "XAGE1D", "XAGE2", "XAGE2B", "XAGE3", "XAGE5",
                     "XIAP", "XK", "YIPF6", "YY2", "ZBTB33", "ZCCHC12",
                     "ZCCHC13", "ZCCHC16", "ZCCHC5", "ZDHHC15", "ZDHHC9", 
                     "ZIC3", "ZMAT1", "ZMYM3", "ZNF157", "ZNF182", "ZNF275",
                     "ZNF280C", "ZNF41", "ZNF449", "ZNF630", "ZNF645", 
                     "ZNF711", "ZNF75D", "ZNF81", "ZXDA", "ZXDB")

# Parse to working table
cosmic_all_NS$type <- ifelse(cosmic_all_NS$gene %in% c(POG_A,POG_X), "POG", "unknown")
cosmic_all_NS$type <- ifelse(cosmic_all_NS$gene %in% c(TSG_A,TSG_X), "TSG", cosmic_all_NS$type)
cosmic_all_NS$chr <- ifelse(cosmic_all_NS$gene %in% c(TSG_X,POG_X), "X", "unknown")
cosmic_all_NS$chr <- ifelse(cosmic_all_NS$gene %in% c(TSG_A,POG_A), "A", cosmic_all_NS$chr)
cosmic_summary <- cosmic_all_NS[(cosmic_all_NS$type != "unknown") | (cosmic_all_NS$chr != "unknown"),]
cosmic_summary$Xinactive <- ifelse(cosmic_summary$gene %in% XCI_inactivated, "yes", "no")
cosmic_summary$escape <- ifelse(cosmic_summary$gene %in% XCI_escapers, "yes", "no")

# TODO FIGURE 1
# Test nonsense versus missense
plot_nonsense_vs_nonsyn <- cosmic_summary |>
  rename(Category = type) |> 
  ggplot(aes(x = log2(nonsyn/nonsynSites), 
             y = log2(nonsense/codons), 
             colour = Category)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("log2[NonSynonymous substitutions / NonSynonymous positions]") +
  ylab("log2[Nonsense substitutions / # codons]") +
  theme_bw()

summary(lm(log2(nonsense/codons) ~ log2(nonsyn/nonsynSites), data = cosmic_summary[cosmic_summary$type == "POG",]))
summary(lm(log2(nonsense/codons) ~ log2(nonsyn/nonsynSites), data = cosmic_summary[cosmic_summary$type == "TSG",]))

summary(lm(log2(nonsense/codons) ~ log2(nonsyn/nonsynSites) * type, data = cosmic_summary))
anova(lm(log2(nonsense/codons) ~ log2(nonsyn/nonsynSites) * type, data = cosmic_summary))
# Justifying the use of missense for OG and nonsense for TSG

ggsave("plots/Mutations.png", plot_nonsense_vs_nonsyn)


# Nonsense for TSG
cosmic_summary |>
  filter(type == "TSG") |>
  ggplot(aes(x = chr,
             y = log2(nonsense/codons))) +
  geom_boxplot() +
  xlab("Chromosome") +
  ylab("log2[Nonsense substitutions / # codons]") +
  theme_bw()
  
  
miss_tsg_x <- cosmic_summary |>  filter((type == "TSG") & (chr == "X")) |>
  mutate(x = log2(nonsense/codons)) |> select(x) |> drop_na() |> unlist() |> as.vector()
miss_tsg_a <- cosmic_summary |>  filter((type == "TSG") & (chr == "A")) |>
  mutate(x = log2(nonsense/codons)) |> select(x) |> drop_na() |> unlist() |> as.vector()
wilcox.test(miss_tsg_x, miss_tsg_a, alternative = "greater")

# Nonsense for TSG, scaled for mutation rate (syn)
cosmic_summary |>
  filter(type == "TSG") |>
  ggplot(aes(x = chr,
             y = log((nonsense/codons)/(syn/synSites)))) +
  geom_boxplot() +
  xlab("Chromosome") +
  ylab("Log-odds ratio (nonsense)") +
  theme_bw()



miss_tsg_x_s <- cosmic_summary |>  filter((type == "TSG") & (chr == "X")) |>
  mutate(x = log((nonsense/codons)/(syn/synSites))) |> select(x) |> drop_na() |> unlist() |> as.vector()
miss_tsg_a_s <- cosmic_summary |>  filter((type == "TSG") & (chr == "A")) |>
  mutate(x = log((nonsense/codons)/(syn/synSites))) |> select(x) |> drop_na() |> unlist() |> as.vector()
wilcox.test(miss_tsg_x_s, miss_tsg_a_s, alternative = "greater")

# Missense for POG
cosmic_summary |>
  filter(type == "POG") |>
  ggplot(aes(x = chr,
             y = log2(nonsyn/nonsynSites))) +
  geom_boxplot() +
  xlab("Chromosome") +
  ylab("Log-odds ratio (nonsynonynous)") +
  theme_bw()


miss_pog_x <- cosmic_summary |>  filter((type == "POG") & (chr == "X")) |>
  mutate(x = log2(nonsense/codons)) |> select(x) |> drop_na() |> unlist() |> as.vector()
miss_pog_a <- cosmic_summary |>  filter((type == "POG") & (chr == "A")) |>
  mutate(x = log2(nonsense/codons)) |> select(x) |> drop_na() |> unlist() |> as.vector()
wilcox.test(miss_pog_x, miss_pog_a, alternative = "less")

# Missense for POG, scaled for mutation rate (syn)
cosmic_summary |>
  filter(type == "POG") |>
  ggplot(aes(x = chr,
             y = LOR)) +
  geom_boxplot() +
  xlab("Chromosome") +
  ylab("Log-odds ratio (nonsynonymous)") +
  theme_bw()


miss_pog_x_s <- cosmic_summary |>  filter((type == "POG") & (chr == "X")) |>
  mutate(x = LOR) |> select(x) |> drop_na() |> unlist() |> as.vector()
miss_pog_a_s <- cosmic_summary |>  filter((type == "POG") & (chr == "A")) |>
  mutate(x = LOR) |> select(x) |> drop_na() |> unlist() |> as.vector()
wilcox.test(miss_pog_x_s, miss_pog_a_s, alternative = "less")




## Split by sex [PAPER]

# Load samples
fsamples <- read.table(gzfile("datasets/female_samples.txt.gz"))$V1
msamples <- read.table(gzfile("datasets/male_samples.txt.gz"))$V1

# Cosmic filter by sex
cosmic_missense_female <- cosmic_missense[cosmic_missense$V5 %in% fsamples,]
cosmic_missense_male <- cosmic_missense[cosmic_missense$V5 %in% msamples,]
cosmic_synonymous_female <- cosmic_synonymous[cosmic_synonymous$V5 %in% fsamples,]
cosmic_synonymous_male <- cosmic_synonymous[cosmic_synonymous$V5 %in% msamples,]
cosmic_nonsense_female <- cosmic_nonsense[cosmic_nonsense$V5 %in% fsamples,]
cosmic_nonsense_male <- cosmic_nonsense[cosmic_nonsense$V5 %in% msamples,]

# Number of missense, nonsense and synonymous variants per gene
# Number of studies with missense, nonsense and synonymous
# Female
cosmic_female_syn <- cosmic_synonymous_female[,c(1,4)] |>
  group_by(V1) |>
  summarise(syn = sum(V4))
cosmic_female_synStudies <- cosmic_synonymous_female[,c(1,4)] |>
  group_by(V1) |> 
  mutate(V4 = 1) |>
  summarise(synStudies = sum(V4))
cosmic_female_nonsyn <- cosmic_missense_female[,c(1,4)] |>
  group_by(V1) |>
  summarise(nonsyn = sum(V4))
cosmic_female_nonsynStudies <- cosmic_missense_female[,c(1,4)] |>
  group_by(V1) |>
  mutate(V4 = 1) |>
  summarise(nonsynStudies = sum(V4))
cosmic_female_nonsense <- cosmic_nonsense_female[,c(1,4)] |>
  group_by(V1) |>
  summarise(nonsense = sum(V4))
cosmic_female_nonsenseStudies <- cosmic_nonsense_female[,c(1,4)] |>
  group_by(V1) |>
  mutate(V4 = 1) |>
  summarise(nonsenseStudies = sum(V4))
cosmic_female_NS <- cosmic_female_syn |>
  left_join(cosmic_female_synStudies, by = "V1") |>
  left_join(cosmic_female_nonsyn, by = "V1") |>
  left_join(cosmic_female_nonsynStudies, by = "V1") |>
  left_join(cosmic_female_nonsense, by = "V1") |>
  left_join(cosmic_female_nonsenseStudies, by = "V1") |>
  left_join(NG86_sites, join_by(V1 == name)) |>
  rename(gene = V1, syn = syn.x, nonsyn = nonsyn.x, synSites = syn.y, nonsynSites = nonsyn.y) |>
  mutate(LOR = log((nonsyn/nonsynSites)/(syn/synSites))) |>
  mutate(SD = sqrt((1/nonsyn)+(1/nonsynSites)+(1/syn)+(1/synSites))) |>
  mutate(zscore = LOR/SD) |>
  mutate(p = pnorm(abs(zscore), mean = 0, sd = 1, lower.tail = FALSE)*2)
cosmic_female_NS$q <- p.adjust(cosmic_female_NS$p, method = "fdr")
# Male
cosmic_male_syn <- cosmic_synonymous_male[,c(1,4)] |>
  group_by(V1) |>
  summarise(syn = sum(V4))
cosmic_male_synStudies <- cosmic_synonymous_male[,c(1,4)] |>
  group_by(V1) |>
  mutate(V4 = 1) |>
  summarise(synStudies = sum(V4))
cosmic_male_nonsyn <- cosmic_missense_male[,c(1,4)] |>
  group_by(V1) |>
  summarise(nonsyn = sum(V4))
cosmic_male_nonsynStudies <- cosmic_missense_male[,c(1,4)] |>
  group_by(V1) |>
  mutate(V4 = 1) |>
  summarise(nonsynStudies = sum(V4))
cosmic_male_nonsense <- cosmic_nonsense_male[,c(1,4)] |>
  group_by(V1) |>
  summarise(nonsense = sum(V4))
cosmic_male_nonsenseStudies <- cosmic_nonsense_male[,c(1,4)] |>
  group_by(V1) |>
  mutate(V4 = 1) |>
  summarise(nonsenseStudies = sum(V4))
cosmic_male_NS <- cosmic_male_syn |>
  left_join(cosmic_male_synStudies, by = "V1") |>
  left_join(cosmic_male_nonsyn, by = "V1") |>
  left_join(cosmic_male_nonsynStudies, by = "V1") |>
  left_join(cosmic_male_nonsense, by = "V1") |>
  left_join(cosmic_male_nonsenseStudies, by = "V1") |>
  left_join(NG86_sites, join_by(V1 == name)) |>
  rename(gene = V1, syn = syn.x, nonsyn = nonsyn.x, synSites = syn.y, nonsynSites = nonsyn.y) |>
  mutate(LOR = log((nonsyn/nonsynSites)/(syn/synSites))) |>
  mutate(SD = sqrt((1/nonsyn)+(1/nonsynSites)+(1/syn)+(1/synSites))) |>
  mutate(zscore = LOR/SD) |>
  mutate(p = pnorm(abs(zscore), mean = 0, sd = 1, lower.tail = FALSE)*2)
cosmic_male_NS$q <- p.adjust(cosmic_male_NS$p, method = "fdr")

# Parse to working table
# female
cosmic_female_NS$type <- ifelse(cosmic_female_NS$gene %in% c(POG_A,POG_X), "POG", "unknown")
cosmic_female_NS$type <- ifelse(cosmic_female_NS$gene %in% c(TSG_A,TSG_X), "TSG", cosmic_female_NS$type)
cosmic_female_NS$chr <- ifelse(cosmic_female_NS$gene %in% c(TSG_X,POG_X), "X", "unknown")
cosmic_female_NS$chr <- ifelse(cosmic_female_NS$gene %in% c(TSG_A,POG_A), "A", cosmic_female_NS$chr)
cosmic_female_summary <- cosmic_female_NS[(cosmic_female_NS$type != "unknown") | (cosmic_female_NS$chr != "unknown"),]
cosmic_female_summary$Xinactive <- ifelse(cosmic_female_summary$gene %in% XCI_inactivated, "yes", "no")
cosmic_female_summary$escape <- ifelse(cosmic_female_summary$gene %in% XCI_escapers, "yes", "no")
# male
cosmic_male_NS$type <- ifelse(cosmic_male_NS$gene %in% c(POG_A,POG_X), "POG", "unknown")
cosmic_male_NS$type <- ifelse(cosmic_male_NS$gene %in% c(TSG_A,TSG_X), "TSG", cosmic_male_NS$type)
cosmic_male_NS$chr <- ifelse(cosmic_male_NS$gene %in% c(TSG_X,POG_X), "X", "unknown")
cosmic_male_NS$chr <- ifelse(cosmic_male_NS$gene %in% c(TSG_A,POG_A), "A", cosmic_male_NS$chr)
cosmic_male_summary <- cosmic_male_NS[(cosmic_male_NS$type != "unknown") | (cosmic_male_NS$chr != "unknown"),]
cosmic_male_summary$Xinactive <- ifelse(cosmic_male_summary$gene %in% XCI_inactivated, "yes", "no")
cosmic_male_summary$escape <- ifelse(cosmic_male_summary$gene %in% XCI_escapers, "yes", "no")


# TODO FIG 2AB
# Nonsense for TSG, scaled for mutation rate (syn)
# Female
female_TSG <- cosmic_female_summary |>
  filter(type == "TSG") |>
  ggplot(aes(x = chr,
             y = log((nonsense/codons)/(syn/synSites)))) +
  geom_boxplot() +
  xlab("Chromosome") +
  ylab("Log-odds ratio (nonsense)") +
  annotate("text", x = Inf, y = Inf, label = "females",
           hjust = 1.1, vjust = 1.5, size = 5)+
  theme_bw()


miss_tsg_x_s_female <- cosmic_female_summary |>  filter((type == "TSG") & (chr == "X")) |>
  mutate(x = log((nonsense/codons)/(syn/synSites))) |> select(x) |> drop_na() |> unlist() |> as.vector()
miss_tsg_a_s_female <- cosmic_female_summary |>  filter((type == "TSG") & (chr == "A")) |>
  mutate(x = log((nonsense/codons)/(syn/synSites))) |> select(x) |> drop_na() |> unlist() |> as.vector()
wilcox.test(miss_tsg_x_s_female, miss_tsg_a_s_female, alternative = "greater")

# Male
male_TSG <- cosmic_male_summary |>
  filter(type == "TSG") |>
  ggplot(aes(x = chr,
             y = log((nonsense/codons)/(syn/synSites)))) +
  geom_boxplot() +
  xlab("Chromosome") +
  ylab("Log-odds ratio (nonsense)") +
  annotate("text", x = Inf, y = Inf, label = "males",
           hjust = 1.1, vjust = 1.5, size = 5)+
  theme_bw()


miss_tsg_x_s_male <- cosmic_male_summary |>  filter((type == "TSG") & (chr == "X")) |>
  mutate(x = log((nonsense/codons)/(syn/synSites))) |> select(x) |> drop_na() |> unlist() |> as.vector()
miss_tsg_a_s_male <- cosmic_male_summary |>  filter((type == "TSG") & (chr == "A")) |>
  mutate(x = log((nonsense/codons)/(syn/synSites))) |> select(x) |> drop_na() |> unlist() |> as.vector()
wilcox.test(miss_tsg_x_s_male, miss_tsg_a_s_male, alternative = "greater")


ggsave("plots/TSG_female_male.png", female_TSG + male_TSG + plot_annotation(tag_levels = "A"))


# TODO FIG 3AB
# Minsense for OG, scaled for mutation rate (syn)
# Female
female_POG <- cosmic_female_summary |>
  filter(type == "POG") |>
  ggplot(aes(x = chr,
             y = LOR)) +
  geom_boxplot() +
  xlab("Chromosome") +
  ylab("Log-odds ratio (nonsynonymous)") +
  annotate("text", x = Inf, y = Inf, label = "females",
           hjust = 1.1, vjust = 1.5, size = 5)+
  theme_bw()

miss_pog_x_s_female <- cosmic_female_summary |>  filter((type == "POG") & (chr == "X")) |>
  mutate(x = LOR) |> select(x) |> drop_na() |> unlist() |> as.vector()
miss_pog_a_s_female <- cosmic_female_summary |>  filter((type == "POG") & (chr == "A")) |>
  mutate(x = LOR) |> select(x) |> drop_na() |> unlist() |> as.vector()
wilcox.test(miss_pog_x_s_female, miss_pog_a_s_female, alternative = "greater")
# Male
male_POG <- cosmic_male_summary |>
  filter(type == "POG") |>
  ggplot(aes(x = chr,
             y = LOR)) +
  geom_boxplot() +
  xlab("Chromosome") +
  ylab("Log-odds ratio (nonsynonymous)") +
  annotate("text", x = Inf, y = Inf, label = "males",
           hjust = 1.1, vjust = 1.5, size = 5)+
  theme_bw()

miss_pog_x_s_male <- cosmic_male_summary |>  filter((type == "POG") & (chr == "X")) |>
  mutate(x = LOR) |> select(x) |> drop_na() |> unlist() |> as.vector()
miss_pog_a_s_male <- cosmic_male_summary |>  filter((type == "POG") & (chr == "A")) |>
  mutate(x = LOR) |> select(x) |> drop_na() |> unlist() |> as.vector()
wilcox.test(miss_pog_x_s_male, miss_pog_a_s_male, alternative = "greater")


ggsave("plots/POG_female_male.png", female_POG + male_POG + plot_annotation(tag_levels = "A"))



# TODO high prolif tissue

# Median proliferation score under 0 (Thorsson)
# KICH  COSU634 Kidney Chromophobe
# PCPG  COSU664 Pheochromocytoma and Paraganglioma
# KIRP  COSU543 Kidney renal papillary cell carcinoma
# THCA  COSU542 Thyroid carcinoma
# PRAD  COSU435 [only male] EXCLUDED (Prostate adenocarcinoma)
# LGG   COSU545 Brain Lower Grade Glioma
# UVM   COSU667 Uveal Melanoma
# KIRC  COSU416 Kidney renal clear cell carcinoma
# ACC   COSU631 Adrenocortical carcinoma
# LIHC  COSU628 Liver hepatocellular carcinoma
# PAAD  COSU629 [only female] EXCLUDED (Ovarian serous cystadenocarcinoma)
# CHOL  COSU662 Cholangiocarcinoma
# MESO  COSU663 Mesothelioma
# LUAD  COSU417 Lung adenocarcinoma

low_proliferation <- c("COSU634", "COSU664", "COSU543",
                       "COSU542", "COSU545", 
                       "COSU667", "COSU416", "COSU631",
                       "COSU628", "COSU662",
                       "COSU663", "COSU417")

# Cosmic filter by sex
cosmic_missense_femaleLP <- cosmic_missense[(cosmic_missense$V5 %in% fsamples) & (cosmic_missense$V2 %in% low_proliferation),]
cosmic_missense_maleLP <- cosmic_missense[(cosmic_missense$V5 %in% msamples) & (cosmic_missense$V2 %in% low_proliferation),]
cosmic_synonymous_femaleLP <- cosmic_synonymous[(cosmic_synonymous$V5 %in% fsamples) & (cosmic_synonymous$V2 %in% low_proliferation),]
cosmic_synonymous_maleLP <- cosmic_synonymous[(cosmic_synonymous$V5 %in% msamples) & (cosmic_synonymous$V2 %in% low_proliferation),]
cosmic_nonsense_femaleLP <- cosmic_nonsense[(cosmic_nonsense$V5 %in% fsamples) & (cosmic_nonsense$V2 %in% low_proliferation),]
cosmic_nonsense_maleLP <- cosmic_nonsense[(cosmic_nonsense$V5 %in% msamples) & (cosmic_nonsense$V2 %in% low_proliferation),]

# Number of missense, nonsense and synonymous variants per gene
# Number of studies with missense, nonsense and synonymous
# FemaleLP
cosmic_femaleLP_syn <- cosmic_synonymous_femaleLP[,c(1,4)] |>
  group_by(V1) |>
  summarise(syn = sum(V4))
cosmic_femaleLP_synStudies <- cosmic_synonymous_femaleLP[,c(1,4)] |>
  group_by(V1) |> 
  mutate(V4 = 1) |>
  summarise(synStudies = sum(V4))
cosmic_femaleLP_nonsyn <- cosmic_missense_femaleLP[,c(1,4)] |>
  group_by(V1) |>
  summarise(nonsyn = sum(V4))
cosmic_femaleLP_nonsynStudies <- cosmic_missense_femaleLP[,c(1,4)] |>
  group_by(V1) |>
  mutate(V4 = 1) |>
  summarise(nonsynStudies = sum(V4))
cosmic_femaleLP_nonsense <- cosmic_nonsense_femaleLP[,c(1,4)] |>
  group_by(V1) |>
  summarise(nonsense = sum(V4))
cosmic_femaleLP_nonsenseStudies <- cosmic_nonsense_femaleLP[,c(1,4)] |>
  group_by(V1) |>
  mutate(V4 = 1) |>
  summarise(nonsenseStudies = sum(V4))
cosmic_femaleLP_NS <- cosmic_femaleLP_syn |>
  left_join(cosmic_femaleLP_synStudies, by = "V1") |>
  left_join(cosmic_femaleLP_nonsyn, by = "V1") |>
  left_join(cosmic_femaleLP_nonsynStudies, by = "V1") |>
  left_join(cosmic_femaleLP_nonsense, by = "V1") |>
  left_join(cosmic_femaleLP_nonsenseStudies, by = "V1") |>
  left_join(NG86_sites, join_by(V1 == name)) |>
  rename(gene = V1, syn = syn.x, nonsyn = nonsyn.x, synSites = syn.y, nonsynSites = nonsyn.y) |>
  mutate(LOR = log((nonsyn/nonsynSites)/(syn/synSites))) |>
  mutate(SD = sqrt((1/nonsyn)+(1/nonsynSites)+(1/syn)+(1/synSites))) |>
  mutate(zscore = LOR/SD) |>
  mutate(p = pnorm(abs(zscore), mean = 0, sd = 1, lower.tail = FALSE)*2)
cosmic_femaleLP_NS$q <- p.adjust(cosmic_femaleLP_NS$p, method = "fdr")
# Male
cosmic_maleLP_syn <- cosmic_synonymous_maleLP[,c(1,4)] |>
  group_by(V1) |>
  summarise(syn = sum(V4))
cosmic_maleLP_synStudies <- cosmic_synonymous_maleLP[,c(1,4)] |>
  group_by(V1) |>
  mutate(V4 = 1) |>
  summarise(synStudies = sum(V4))
cosmic_maleLP_nonsyn <- cosmic_missense_maleLP[,c(1,4)] |>
  group_by(V1) |>
  summarise(nonsyn = sum(V4))
cosmic_maleLP_nonsynStudies <- cosmic_missense_maleLP[,c(1,4)] |>
  group_by(V1) |>
  mutate(V4 = 1) |>
  summarise(nonsynStudies = sum(V4))
cosmic_maleLP_nonsense <- cosmic_nonsense_maleLP[,c(1,4)] |>
  group_by(V1) |>
  summarise(nonsense = sum(V4))
cosmic_maleLP_nonsenseStudies <- cosmic_nonsense_maleLP[,c(1,4)] |>
  group_by(V1) |>
  mutate(V4 = 1) |>
  summarise(nonsenseStudies = sum(V4))
cosmic_maleLP_NS <- cosmic_maleLP_syn |>
  left_join(cosmic_maleLP_synStudies, by = "V1") |>
  left_join(cosmic_maleLP_nonsyn, by = "V1") |>
  left_join(cosmic_maleLP_nonsynStudies, by = "V1") |>
  left_join(cosmic_maleLP_nonsense, by = "V1") |>
  left_join(cosmic_maleLP_nonsenseStudies, by = "V1") |>
  left_join(NG86_sites, join_by(V1 == name)) |>
  rename(gene = V1, syn = syn.x, nonsyn = nonsyn.x, synSites = syn.y, nonsynSites = nonsyn.y) |>
  mutate(LOR = log((nonsyn/nonsynSites)/(syn/synSites))) |>
  mutate(SD = sqrt((1/nonsyn)+(1/nonsynSites)+(1/syn)+(1/synSites))) |>
  mutate(zscore = LOR/SD) |>
  mutate(p = pnorm(abs(zscore), mean = 0, sd = 1, lower.tail = FALSE)*2)
cosmic_maleLP_NS$q <- p.adjust(cosmic_maleLP_NS$p, method = "fdr")

# Parse to working table
# femaleLP
cosmic_femaleLP_NS$type <- ifelse(cosmic_femaleLP_NS$gene %in% c(POG_A,POG_X), "POG", "unknown")
cosmic_femaleLP_NS$type <- ifelse(cosmic_femaleLP_NS$gene %in% c(TSG_A,TSG_X), "TSG", cosmic_femaleLP_NS$type)
cosmic_femaleLP_NS$chr <- ifelse(cosmic_femaleLP_NS$gene %in% c(TSG_X,POG_X), "X", "unknown")
cosmic_femaleLP_NS$chr <- ifelse(cosmic_femaleLP_NS$gene %in% c(TSG_A,POG_A), "A", cosmic_femaleLP_NS$chr)
cosmic_femaleLP_summary <- cosmic_femaleLP_NS[(cosmic_femaleLP_NS$type != "unknown") | (cosmic_femaleLP_NS$chr != "unknown"),]
cosmic_femaleLP_summary$Xinactive <- ifelse(cosmic_femaleLP_summary$gene %in% XCI_inactivated, "yes", "no")
cosmic_femaleLP_summary$escape <- ifelse(cosmic_femaleLP_summary$gene %in% XCI_escapers, "yes", "no")
# maleLP
cosmic_maleLP_NS$type <- ifelse(cosmic_maleLP_NS$gene %in% c(POG_A,POG_X), "POG", "unknown")
cosmic_maleLP_NS$type <- ifelse(cosmic_maleLP_NS$gene %in% c(TSG_A,TSG_X), "TSG", cosmic_maleLP_NS$type)
cosmic_maleLP_NS$chr <- ifelse(cosmic_maleLP_NS$gene %in% c(TSG_X,POG_X), "X", "unknown")
cosmic_maleLP_NS$chr <- ifelse(cosmic_maleLP_NS$gene %in% c(TSG_A,POG_A), "A", cosmic_maleLP_NS$chr)
cosmic_maleLP_summary <- cosmic_maleLP_NS[(cosmic_maleLP_NS$type != "unknown") | (cosmic_maleLP_NS$chr != "unknown"),]
cosmic_maleLP_summary$Xinactive <- ifelse(cosmic_maleLP_summary$gene %in% XCI_inactivated, "yes", "no")
cosmic_maleLP_summary$escape <- ifelse(cosmic_maleLP_summary$gene %in% XCI_escapers, "yes", "no")

# For Low Prolif
# Minsense for OG, scaled for mutation rate (syn)
# Female
plot_POG_female_LP <- cosmic_femaleLP_summary |>
  filter(type == "POG") |>
  ggplot(aes(x = chr,
             y = LOR)) +
  geom_boxplot() +
  xlab("Chromosome") +
  ylab("Log-odds ratio (nonsynonynous)") +
  theme_bw()

miss_pog_x_s_femaleLP <- cosmic_femaleLP_summary |>  filter((type == "POG") & (chr == "X")) |>
  mutate(x = LOR) |> select(x) |> drop_na() |> unlist() |> as.vector()
miss_pog_a_s_femaleLP <- cosmic_femaleLP_summary |>  filter((type == "POG") & (chr == "A")) |>
  mutate(x = LOR) |> select(x) |> drop_na() |> unlist() |> as.vector()
wilcox.test(miss_pog_x_s_femaleLP, miss_pog_a_s_femaleLP, alternative = "less")
# Male
plot_POG_male_LP <- cosmic_maleLP_summary |>
  filter(type == "POG") |>
  ggplot(aes(x = chr,
             y = LOR)) +
  geom_boxplot() +
  xlab("Chromosome") +
  ylab("Log-odds ratio (nonsynonynous)") +
  theme_bw()

miss_pog_x_s_maleLP <- cosmic_maleLP_summary |>  filter((type == "POG") & (chr == "X")) |>
  mutate(x = LOR) |> select(x) |> drop_na() |> unlist() |> as.vector()
miss_pog_a_s_maleLP <- cosmic_maleLP_summary |>  filter((type == "POG") & (chr == "A")) |>
  mutate(x = LOR) |> select(x) |> drop_na() |> unlist() |> as.vector()
wilcox.test(miss_pog_x_s_maleLP, miss_pog_a_s_maleLP, alternative = "less")

ggsave("plots/LorProlif.png", plot_POG_female_LP + plot_POG_male_LP + plot_annotation(tag_levels = "A"))

# 
# 
# # Two-way, SRH-test, interaction
# pog_all <- c(miss_pog_x_s_femaleLP, 
#              miss_pog_a_s_femaleLP,
#              miss_pog_x_s_maleLP,
#              miss_pog_a_s_maleLP)
# pog_chr <- c(rep("X",length(miss_pog_x_s_femaleLP)), 
#              rep("A",length(miss_pog_a_s_femaleLP)), 
#              rep("X",length(miss_pog_x_s_maleLP)), 
#              rep("A",length(miss_pog_a_s_maleLP)))
# pog_sex <- c(rep("F",length(c(miss_pog_x_s_femaleLP, miss_pog_a_s_femaleLP))),
#              rep("M",length(c(miss_pog_x_s_maleLP, miss_pog_a_s_maleLP))))
# summary(lm(rank(pog_all) ~ pog_sex + pog_chr))
# anova(lm(rank(pog_all) ~ pog_sex + pog_chr))






# Mutations (unique  about arrival)

unique_first_cosmic_female <- cosmic |>
  filter(V5 %in% fsamples) |>
  filter(((V3 == "stop_gained") & (V1 %in% TSG_all$GENE_SYMBOL)) |
           ((V3 == "missense_variant") & (V1 %in% POG_all$GENE_SYMBOL))) |>
  mutate(V6 =  ifelse(V1 %in% POG_all$GENE_SYMBOL, "POG", "unknown")) |>
  mutate(V6 =  ifelse(V1 %in% TSG_all$GENE_SYMBOL, "TSG", V6)) |>
  mutate(V7 = ifelse(V1 %in% c(TSG_X,POG_X), "X", "unknown")) |>
  mutate(V7 = ifelse(V1 %in% c(TSG_A,POG_A), "A", V7)) |>
  add_count(V5, name = "V8") |>
  filter(V8 == 1)
unique_first_cosmic_male <- cosmic |>
  filter(V5 %in% msamples) |>
  filter(((V3 == "stop_gained") & (V1 %in% TSG_all$GENE_SYMBOL)) |
           ((V3 == "missense_variant") & (V1 %in% POG_all$GENE_SYMBOL))) |>
  mutate(V6 =  ifelse(V1 %in% POG_all$GENE_SYMBOL, "POG", "unknown")) |>
  mutate(V6 =  ifelse(V1 %in% TSG_all$GENE_SYMBOL, "TSG", V6)) |>
  mutate(V7 = ifelse(V1 %in% c(TSG_X,POG_X), "X", "unknown")) |>
  mutate(V7 = ifelse(V1 %in% c(TSG_A,POG_A), "A", V7)) |>
  add_count(V5, name = "V8") |>
  filter(V8 == 1)

unique_first_cosmic_female |>
  select(V6, V7) |>
  table()
unique_first_cosmic_male |>
  select(V6, V7) |>
  table()

fisher.test(matrix(c(29,21, 28,31), ncol = 2), alternative = "greater")
fisher.test(matrix(c(2206,300, 2434,350), ncol = 2), alternative = "greater")
# not significant but indicative!

# "mutable" sites
mutable_TSG_A <- cosmic_summary |>
  filter((chr == "A") & (type == "TSG")) |>
  pull(codons) |>
  sum( na.rm=TRUE)
mutable_TSG_X <- cosmic_summary |>
  filter((chr == "X") & (type == "TSG")) |>
  pull(codons) |>
  sum( na.rm=TRUE)
mutable_POG_A <- cosmic_summary |>
  filter((chr == "A") & (type == "POG")) |>
  pull(nonsynSites) |>
  sum( na.rm=TRUE)
mutable_POG_X <- cosmic_summary |>
  filter((chr == "X") & (type == "POG")) |>
  pull(nonsynSites) |>
  sum( na.rm=TRUE)

c(mutable_TSG_A, mutable_TSG_X, mutable_POG_A, mutable_POG_X)

E_TSG_A_f = (177421/(177421+375125))*(300+2206)
E_POG_A_f = (375125/(177421+375125))*(300+2206)
G_A_f <- 2 * ( 300*log(300/E_TSG_A_f) + 2206*log(2206/E_POG_A_f) )

E_TSG_X_f = (12962/(12962+7254))*(21+29)
E_POG_X_f = (7254/(12962+7254))*(21+29)
G_X_f <- 2 * ( 21*log(21/E_TSG_X_f) + 29*log(29/E_POG_X_f) )

E_TSG_A_m = (177421/(177421+375125))*(350+2434)
E_POG_A_m = (375125/(177421+375125))*(350+2434)
G_A_m <- 2 * ( 350*log(350/E_TSG_A_m) + 2434*log(2434/E_POG_A_m) )

E_TSG_X_m = (12962/(12962+7254))*(31+28)
E_POG_X_m = (7254/(12962+7254))*(31+28)
G_X_m <- 2 * ( 31*log(31/E_TSG_X_m) + 28*log(28/E_POG_X_m) )


write.table(data.frame(Chromosome = c("A", "X", "A", "X"),
           Sex = c("female", "female", "male", "male"),
           observed_tsg = c(300, 21, 350, 31),
           expected_tsg = c(E_TSG_A_f, E_TSG_X_f, E_TSG_A_m, E_TSG_X_m),
           observed_pog = c(2206, 29, 2434, 28),
           expected_pog = c(E_POG_A_f, E_POG_X_f, E_POG_A_m, E_POG_X_m),
           p = p.adjust(c(pchisq(G_A_f, df = 1, lower.tail = FALSE), 
                 pchisq(G_X_f, df = 1, lower.tail = FALSE),
                 pchisq(G_A_m, df = 1, lower.tail = FALSE),
                 pchisq(G_X_m, df = 1, lower.tail = FALSE)), method = "bonferroni")
           ),
           file = "results/table_arrival.tab")






# Escapers (not enough datapoints)

# Female no-escapers
cosmic_female_summary |>
  filter(type == "TSG") |>
  ggplot(aes(x = paste(chr, Xinactive),
             y = log((nonsense/codons)/(syn/synSites)))) +
  geom_boxplot()
miss_tsg_x_s_female_noXin <- cosmic_female_summary |>  filter((type == "TSG") & (chr == "X") & (Xinactive == "no")) |>
  mutate(x = log((nonsense/codons)/(syn/synSites))) |> select(x) |> drop_na() |> unlist() |> as.vector()
miss_tsg_x_s_female_Xin <- cosmic_female_summary |>  filter((type == "TSG") & (chr == "X") & (Xinactive == "yes")) |>
  mutate(x = log((nonsense/codons)/(syn/synSites))) |> select(x) |> drop_na() |> unlist() |> as.vector()
# inactive vs nonInactive
wilcox.test(miss_tsg_x_s_female_Xin, miss_tsg_x_s_female_noXin, alternative = "greater")
# x inactive vs A
wilcox.test(miss_tsg_x_s_female_Xin, miss_tsg_a_s_female, alternative = "greater")







q()

