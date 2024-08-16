#POG HAP VS DIP
# Load required packages
library(ggplot2)
library(gridExtra)
library(patchwork)


# Function to simulate genotype frequencies and record allele frequencies over time
simulate_somatic_genotype_frequencies <- function(initial_AA_freq = 1, mutation_rate = 0.01, heterozygosity_coefficient = 1, selection_coefficient = 0.1, num_generations = 1000, ploidy = 2) {
  if(ploidy == 2){
    # Create empty data frame to store genotype and allele frequencies
    genotypes <- data.frame(generation = 0:num_generations,
                            AA = numeric(num_generations + 1),
                            Aa = numeric(num_generations + 1),
                            aa = numeric(num_generations + 1),
                            A_freq = numeric(num_generations + 1),
                            a_freq = numeric(num_generations + 1))
    
    # Set initial genotype frequencies
    genotypes$AA[1] <- initial_AA_freq
    genotypes$Aa[1] <- 1 - initial_AA_freq
    genotypes$aa[1] <- 0
    
    # Set initial allele frequencies
    genotypes$A_freq[1] <- initial_AA_freq + 0.5 * (1 - initial_AA_freq)
    genotypes$a_freq[1] <- 1 - genotypes$A_freq[1]
    
    # Iterate over generations
    for (t in 1:num_generations) {
      # Calculate genotype frequencies after mutation
      genotypes$AA[t + 1] <- genotypes$AA[t] * (1 - 2 * mutation_rate)
      genotypes$Aa[t + 1] <- genotypes$Aa[t] * (1 - mutation_rate) + 2 * mutation_rate * genotypes$AA[t]
      genotypes$aa[t + 1] <- genotypes$aa[t] + mutation_rate * genotypes$Aa[t]
      
      # Apply selection to genotype frequencies
      fitness <- c(1, 1 + heterozygosity_coefficient * selection_coefficient, 1 + selection_coefficient)
      genotypes$AA[t + 1] <- genotypes$AA[t + 1] * fitness[1]
      genotypes$Aa[t + 1] <- genotypes$Aa[t + 1] * fitness[2]
      genotypes$aa[t + 1] <- genotypes$aa[t + 1] * fitness[3]
      
      # Normalize genotype frequencies
      freq_sum <- sum(genotypes$AA[t + 1], genotypes$Aa[t + 1], genotypes$aa[t + 1])
      genotypes$AA[t + 1] <- genotypes$AA[t + 1] / freq_sum
      genotypes$Aa[t + 1] <- genotypes$Aa[t + 1] / freq_sum
      genotypes$aa[t + 1] <- genotypes$aa[t + 1] / freq_sum
      
      # Calculate allele frequencies
      genotypes$A_freq[t + 1] <- genotypes$AA[t + 1] + 0.5 * genotypes$Aa[t + 1]
      genotypes$a_freq[t + 1] <- 1 - genotypes$A_freq[t + 1]
    }
    
    return(genotypes)
    
  } else if (ploidy == 1) {
    # Create empty data frame to store genotype and allele frequencies
    genotypes <- data.frame(generation = 0:num_generations,
                            A = numeric(num_generations + 1),
                            a = numeric(num_generations + 1))
    
    # Set initial genotype frequencies
    genotypes$A[1] <- initial_AA_freq
    genotypes$a[1] <- 1 - initial_AA_freq
    
    # Iterate over generations
    for (t in 1:num_generations) {
      # Calculate genotype frequencies after mutation
      genotypes$A[t + 1] <- genotypes$A[t] * (1 - mutation_rate)
      genotypes$a[t + 1] <- genotypes$a[t] + mutation_rate * genotypes$A[t]
      
      # Apply selection to genotype frequencies
      fitness <- c(1, 1 + selection_coefficient)
      genotypes$A[t + 1] <- genotypes$A[t + 1] * fitness[1]
      genotypes$a[t + 1] <- genotypes$a[t + 1] * fitness[2]
      
      # Normalize genotype frequencies
      freq_sum <- sum(genotypes$A[t + 1], genotypes$a[t + 1])
      genotypes$A[t + 1] <- genotypes$A[t + 1] / freq_sum
      genotypes$a[t + 1] <- genotypes$a[t + 1] / freq_sum
    }
    
    return(genotypes)
    
  } else {
    cat("ERROR: ploidy value not implemented")
    q(save = FALSE)
  }

}



# Set selection strengths
selection_strengths <- c(0.01, 0.05, 0.1)

# FIGURE: POG Dip vs Hap
# Iterate over selection strengths
plots <- list()
for (selection_coeff in selection_strengths) {
  # Simulate genotype and allele frequencies over time for heterozygosity_coefficient = 0
  genotype_data_pog_dip <- simulate_somatic_genotype_frequencies(initial_AA_freq = 1, mutation_rate = 0.01, heterozygosity_coefficient = 1, selection_coefficient = selection_coeff, num_generations = 250, ploidy = 2)
  genotype_data_pog_hap <- simulate_somatic_genotype_frequencies(initial_AA_freq = 1, mutation_rate = 0.01, heterozygosity_coefficient = 1, selection_coefficient = selection_coeff, num_generations = 250, ploidy = 1)
 
  # Create plots for current selection strength
  plot_temp <- ggplot() +
    geom_line(data = genotype_data_pog_dip, aes(x = generation, y = Aa + aa), linewidth = 1, linetype = "solid") +
    geom_line(data = genotype_data_pog_hap, aes(x = generation, y = a), linewidth = 1, linetype = "dashed") +
    labs(x = "Generation", y = "Frequency (c)") +
     scale_color_manual(values = c("Aa + aa h=1" = "black", "a" = "black"),
                        name = "Genotypes",
                        labels = c("Aa + aa h=1" = "a", "a" = "Aa + aa h=1"),
                        guide = guide_legend(override.aes = list(linetype = c("solid", "dashed"))))+
#      theme_gray(base_size = 15) +
    theme(plot.title = element_text(size=14)) +
    theme(panel.background = element_rect(fill = 'white', color = 'black'),
          panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
    labs(title = paste("s = ", selection_coeff, sep = ""), tag = "b") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag = element_text()) +
    theme(legend.position = c(0.5, 0.25),
        legend.background = element_rect(fill = "white", color = "black")) +
    theme(aspect.ratio=3/3)
  
  # Add plot to the list
  plots[[length(plots) + 1]] <- plot_temp
}
ggsave(filename = "plots/Figure_1.png", width = 18, units = "cm", plot  = (plots[[1]] & xlab(NULL)) + (plots[[2]] & ylab(NULL)) + (plots[[3]] & ylab(NULL) & xlab(NULL)) + plot_annotation(tag_levels = "A"))



# FIGURE: TSG Dip vs Hap
# Iterate over selection strengths
plots <- list()
for (selection_coeff in selection_strengths) {
  # Simulate genotype and allele frequencies over time for heterozygosity_coefficient = 0
  genotype_data_tsg_dip <- simulate_somatic_genotype_frequencies(initial_AA_freq = 1, mutation_rate = 0.01, heterozygosity_coefficient = 0, selection_coefficient = selection_coeff, num_generations = 250, ploidy = 2)
  genotype_data_tsg_hap <- simulate_somatic_genotype_frequencies(initial_AA_freq = 1, mutation_rate = 0.01, heterozygosity_coefficient = 0, selection_coefficient = selection_coeff, num_generations = 250, ploidy = 1)
 
  # Create plots for current selection strength
  plot_temp <- ggplot() +
    geom_line(data = genotype_data_tsg_dip, aes(x = generation, y = aa), linewidth = 1, linetype = "solid") +
    geom_line(data = genotype_data_tsg_hap, aes(x = generation, y = a), linewidth = 1, linetype = "dashed") +
    labs(x = "Generation", y = "Frequency (c)") +
     scale_color_manual(values = c("Aa + aa h=1" = "black", "a" = "black"),
                        name = "Genotypes",
                        labels = c("Aa + aa h=1" = "a", "a" = "Aa + aa h=1"),
                        guide = guide_legend(override.aes = list(linetype = c("solid", "dashed"))))+
#      theme_gray(base_size = 15) +
    theme(plot.title = element_text(size=14)) +
    theme(panel.background = element_rect(fill = 'white', color = 'black'),
          panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
    labs(title = paste("s = ", selection_coeff, sep = ""), tag = "b") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag = element_text()) +
    theme(legend.position = c(0.5, 0.25),
        legend.background = element_rect(fill = "white", color = "black")) +
    theme(aspect.ratio=3/3)
  
  # Add plot to the list
  plots[[length(plots) + 1]] <- plot_temp
}
ggsave(filename = "plots/Figure_3.png", width = 18, units = "cm", plot  = (plots[[1]] & xlab(NULL)) + (plots[[2]] & ylab(NULL)) + (plots[[3]] & ylab(NULL) & xlab(NULL)) + plot_annotation(tag_levels = "A"))



# Haf fixation times analysis
selection_strengths_fine <- seq(0,1,0.01)
half_times_selection_scan <- function(mutation_rate = 0.01) {
  mut_rate <- mutation_rate

  thalfs_TSG_dip <- vector()
  thalfs_TSG_hap <- vector()
  thalfs_POG_dip <- vector()
  thalfs_POG_hap <- vector()

  # Iterate over selection strengths
  for (selection_coeff in selection_strengths_fine) {
    # Simulate genotype and allele frequencies over time for heterozygosity_coefficient = 0
    genotype_data_tsg_dip <- simulate_somatic_genotype_frequencies(initial_AA_freq = 1, mutation_rate = mut_rate, heterozygosity_coefficient = 0, selection_coefficient = selection_coeff, num_generations = 250, ploidy = 2)
    genotype_data_tsg_hap <- simulate_somatic_genotype_frequencies(initial_AA_freq = 1, mutation_rate = mut_rate, heterozygosity_coefficient = 0, selection_coefficient = selection_coeff, num_generations = 250, ploidy = 1)
    thalfs_TSG_dip <- c(thalfs_TSG_dip, min(which(genotype_data_tsg_dip$aa > 0.5)))
    thalfs_TSG_hap <- c(thalfs_TSG_hap, min(which(genotype_data_tsg_hap$a > 0.5)))
      
    # Simulate genotype and allele frequencies over time for heterozygosity_coefficient = 1
    genotype_data_pog_dip <- simulate_somatic_genotype_frequencies(initial_AA_freq = 1, mutation_rate = mut_rate, heterozygosity_coefficient = 1, selection_coefficient = selection_coeff, num_generations = 250, ploidy = 2)
    genotype_data_pog_hap <- simulate_somatic_genotype_frequencies(initial_AA_freq = 1, mutation_rate = mut_rate, heterozygosity_coefficient = 1, selection_coefficient = selection_coeff, num_generations = 250, ploidy = 1)  
    thalfs_POG_dip <- c(thalfs_POG_dip, min(which((genotype_data_pog_dip$aa + genotype_data_pog_dip$Aa) > 0.5)))
    thalfs_POG_hap <- c(thalfs_POG_hap, min(which(genotype_data_pog_hap$a > 0.5)))  
  }

    return(list(thalfs_TSG_dip, thalfs_TSG_hap, thalfs_POG_dip, thalfs_POG_hap))

}

thalfs_list_u001 <- half_times_selection_scan(mutation_rate = 0.01)
thalfs_list_u005 <- half_times_selection_scan(mutation_rate = 0.05)
thalfs_list_u01 <- half_times_selection_scan(mutation_rate = 0.1)

# FIGURE sel scan POG
plot_POG_001 <- ggplot() +
  geom_line(aes(x = selection_strengths_fine, y = thalfs_list_u001[[3]]), linewidth = 1, linetype = "solid") +
  geom_line(aes(x = selection_strengths_fine, y = thalfs_list_u001[[4]]), linewidth = 1, linetype = "dashed") +
  labs(x = "Selection coefficient", y = "t(c>=0.5)") +
  theme(plot.title = element_text(size=14)) +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  labs(title = "u = 0.01") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(aspect.ratio=3/3)
plot_POG_005 <- ggplot() +
  geom_line(aes(x = selection_strengths_fine, y = thalfs_list_u005[[3]]), linewidth = 1, linetype = "solid") +
  geom_line(aes(x = selection_strengths_fine, y = thalfs_list_u005[[4]]), linewidth = 1, linetype = "dashed") +
  labs(x = "Selection coefficient", y = "t(c>=0.5)") +
  theme(plot.title = element_text(size=14)) +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  labs(title = "u = 0.05") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(aspect.ratio=3/3)
plot_POG_01 <- ggplot() +
  geom_line(aes(x = selection_strengths_fine, y = thalfs_list_u01[[3]]), linewidth = 1, linetype = "solid") +
  geom_line(aes(x = selection_strengths_fine, y = thalfs_list_u01[[4]]), linewidth = 1, linetype = "dashed") +
  labs(x = "Selection coefficient", y = "t(c>=0.5)") +
  theme(plot.title = element_text(size=14)) +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  labs(title = "u = 0.1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(aspect.ratio=3/3)

ggsave(filename = "plots/Figure_2.png", width = 18, units = "cm", plot  = (plot_POG_001 & xlab(NULL)) + (plot_POG_005 & ylab(NULL)) + (plot_POG_01 & ylab(NULL) & xlab(NULL)) + plot_annotation(tag_levels = "A"))



# FIGURE sel scan TSG
plot_TSG_001 <- ggplot() +
  geom_line(aes(x = selection_strengths_fine, y = thalfs_list_u001[[1]]), linewidth = 1, linetype = "solid") +
  geom_line(aes(x = selection_strengths_fine, y = thalfs_list_u001[[2]]), linewidth = 1, linetype = "dashed") +
  labs(x = "Selection coefficient", y = "t(c>=0.5)") +
  theme(plot.title = element_text(size=14)) +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  labs(title = "u = 0.01") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(aspect.ratio=3/3)
plot_TSG_005 <- ggplot() +
  geom_line(aes(x = selection_strengths_fine, y = thalfs_list_u005[[1]]), linewidth = 1, linetype = "solid") +
  geom_line(aes(x = selection_strengths_fine, y = thalfs_list_u005[[2]]), linewidth = 1, linetype = "dashed") +
  labs(x = "Selection coefficient", y = "t(c>=0.5)") +
  theme(plot.title = element_text(size=14)) +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  labs(title = "u = 0.05") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(aspect.ratio=3/3)
plot_TSG_01 <- ggplot() +
  geom_line(aes(x = selection_strengths_fine, y = thalfs_list_u01[[1]]), linewidth = 1, linetype = "solid") +
  geom_line(aes(x = selection_strengths_fine, y = thalfs_list_u01[[2]]), linewidth = 1, linetype = "dashed") +
  labs(x = "Selection coefficient", y = "t(c>=0.5)") +
  theme(plot.title = element_text(size=14)) +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
  labs(title = "u = 0.1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(aspect.ratio=3/3)

ggsave(filename = "plots/Figure_4.png", width = 18, units = "cm", plot  = (plot_TSG_001 & xlab(NULL)) + (plot_TSG_005 & ylab(NULL)) + (plot_TSG_01 & ylab(NULL) & xlab(NULL)) + plot_annotation(tag_levels = "A"))



quit()
