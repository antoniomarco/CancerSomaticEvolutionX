#!/usr/bin/bash


# Evaluate the approximated forms of the change in frequencies equation.
# Generates Supplementary Figure 1, 2 and 3
Rscript scripts/01_evaluate_approximation_change_function.R

# Run the computer model
# Generates Figures 1, 2, 3 and 4
Rscript scripts/02_run_model_somatic_evolution.R

# ## Download and parse COSMIC data. REGISTRATION AND LOG IN REQUIRED!!
# # https://cancer.sanger.ac.uk/cosmic/download#/
# # Beta Data Downloads (release v99, 28th November 2023)
# tar -xvf Cosmic_GenomeScreensMutant_Tsv_v99_GRCh38.tar
# zcat Cosmic_GenomeScreensMutant_v99_GRCh38.tsv.gz | awk -F'\t' '{print $1 "\t" $19 "" $20 "\t" $12 "\t" $4}' | sort -k 1,2 | uniq -c | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $5}' | grep -v GENE_SYMBOL | gzip > datasets/Cosmic_GenomeScreen_parsed.tab.gz
# # Census
# tar -xvf Cosmic_CancerGeneCensus_Tsv_v99_GRCh38.tar
# mv Cosmic_CancerGeneCensus_v99_GRCh38.tsv.gz datasets/
# # Sample information
# tar -xvf Cosmic_Sample_Tsv_v99_GRCh38.tar
# zcat Cosmic_Sample_v99_GRCh38.tsv.gz | awk -F '\t' '$14 == "f" {print $1}' | sort | uniq | gzip > datasets/female_samples.txt.gz
# zcat Cosmic_Sample_v99_GRCh38.tsv.gz | awk -F '\t' '$14 == "m" {print $1}' | sort | uniq | gzip > datasets/male_samples.txt.gz
# # Remove temp files
# rm Cosmic_CancerGeneCensus_Tsv_v99_GRCh38.tar README_Cosmic_CancerGeneCensus_v99_GRCh38.txt README_Cosmic_GenomeScreensMutant_v99_GRCh38.txt Cosmic_GenomeScreensMutant_Tsv_v99_GRCh38.tar README_Cosmic_Sample_v99_GRCh38.txt Cosmic_GenomeScreensMutant_v99_GRCh38.tsv.gz Cosmic_Sample_Tsv_v99_GRCh38.tar Cosmic_Sample_v99_GRCh38.tsv.gz


# Analyze Genomic Mutation Screens
# Generates Figures 5 and 6
Rscript scripts/03_empirial_rates_from_genomic_screens.R


exit 0
