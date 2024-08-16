# Exact versus approx.
# Load required packages
library(ggplot2)
library(gridExtra)
library(patchwork)


# deltaC haploid

delta_c_h <- function(c=0, u=0.01, s=0.1, approx = FALSE){
  if(approx==FALSE){
    return(((1-c)*u+(1-c^2-c*(1-c)*u)*s)/(1+(c+(1-c)*u)*s))
  }else if(approx==TRUE){
    return(((1-c)*u+(1-c^2)*s)/(1+c*s))
  }else{
    return("Not TRUE/FALSE in 'approx' value.")
  }
}

plots <- list()
for(custom_u in c(0.001, 0.01, 0.1)){
  for(custom_s in c(0.001, 0.01, 0.1)){  
    c_exact_approx <- data.frame(c = seq(0,0.99, 0.01), exact = delta_c_h(c=seq(0,0.99, 0.01), u = custom_u, s = custom_s, approx = FALSE), approx = delta_c_h(c=seq(0,0.99, 0.01), u = custom_u, s = custom_s, approx = TRUE))
    plot_temp <- ggplot() +
    geom_line(data = c_exact_approx, aes(x = c, y = exact), linewidth = 1, linetype = "solid", color = "darkblue") +
    geom_line(data = c_exact_approx, aes(x = c, y = approx), linewidth = 1, linetype = "dashed", color = "darkred") + 
    labs(x = "Frequency (c)", y = "Change in frequency (\u0394 c)") +
     scale_color_manual(values = c("Exact" = "blue", "Approximation" = "red"),
                        name = "Genotypes",
                        labels = c("Aa + aa h=1" = "a", "a" = "Aa + aa h=1"),
                        guide = guide_legend(override.aes = list(linetype = c("solid", "dashed"))))+
#      theme_gray(base_size = 15) +
    theme(plot.title = element_text(size=11)) +
    theme(panel.background = element_rect(fill = 'white', color = 'black'),
          panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
    labs(title = paste("s = ", custom_s, "; u = ", custom_u, sep = ""), tag = "b") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag = element_text()) +
    theme(legend.position = c(0.5, 0.25),
        legend.background = element_rect(fill = "white", color = "black")) +
    theme(aspect.ratio=3/3)
    # Add plot to the list
    plots[[length(plots) + 1]] <- plot_temp
  }
}
ggsave(filename = "plots/Supplementary_Figure_1.png", plot  = (plots[[1]] + plots[[2]] + plots[[3]]) / (plots[[4]] + plots[[5]] + plots[[6]]) / (plots[[7]] + plots[[8]] + plots[[9]]) + plot_annotation(tag_levels = "A"))



# C diploid dominant

delta_c_d <- function(c=0, u=0.01, s=0.1, approx = FALSE){
  if(approx==FALSE){
    return((2*(1-c)*u+(1-c^2-c*2*(1-c)*u)*s)/(1+(c+2*(1-c)*u)*s))
  }else if(approx==TRUE){
    return((2*(1-c)*u+(1-c^2)*s)/(1+c*s))
  }else{
    return("Not TRUE/FALSE in 'approx' value.")
  }
}


plots <- list()
for(custom_u in c(0.001, 0.01, 0.1)){
  for(custom_s in c(0.001, 0.01, 0.1)){
    c_exact_approx <- data.frame(c = seq(0,0.99, 0.01), exact = delta_c_d(c=seq(0,0.99, 0.01), u = custom_u, s = custom_s, approx = FALSE), approx = delta_c_d(c=seq(0,0.99, 0.01), u = custom_u, s = custom_s, approx = TRUE))
    plot_temp <- ggplot() +
    geom_line(data = c_exact_approx, aes(x = c, y = exact), linewidth = 1, linetype = "solid", color = "darkblue") +
    geom_line(data = c_exact_approx, aes(x = c, y = approx), linewidth = 1, linetype = "dashed", color = "darkred") + 
    labs(x = "Frequency (c)", y = "Change in frequency (\u0394 c)") +
     scale_color_manual(values = c("Exact" = "blue", "Approximation" = "red"),
                        name = "Genotypes",
                        labels = c("Aa + aa h=1" = "a", "a" = "Aa + aa h=1"),
                        guide = guide_legend(override.aes = list(linetype = c("solid", "dashed"))))+
#      theme_gray(base_size = 15) +
    theme(plot.title = element_text(size=11)) +
    theme(panel.background = element_rect(fill = 'white', color = 'black'),
          panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
    labs(title = paste("s = ", custom_s, "; u = ", custom_u, sep = ""), tag = "b") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag = element_text()) +
    theme(legend.position = c(0.5, 0.25),
        legend.background = element_rect(fill = "white", color = "black")) +
    theme(aspect.ratio=3/3)
  
  # Add plot to the list
  plots[[length(plots) + 1]] <- plot_temp
    }
}
ggsave(filename = "plots/Supplementary_Figure_2.png", plot  = (plots[[1]] + plots[[2]] + plots[[3]]) / (plots[[4]] + plots[[5]] + plots[[6]]) / (plots[[7]] + plots[[8]] + plots[[9]]) + plot_annotation(tag_levels = "A"))



# C diploid recessive

delta_c_dr <- function(c=0, u=0.01, s=0.1, approx = FALSE){
  if(approx==FALSE){
    return((2*u^2*(1-c)+(c-c^2+2*u^2-4*c*u^2+2*c^2*u^2)*s)/(1+(c+2*u^2*(1-c))*s))
  }else if(approx==TRUE){
    return((c*(1-c)*s)/(1+c*s))
  }else{
    return("Not TRUE/FALSE in 'approx' value.")
  }
}


plots <- list()
for(custom_u in c(0.001, 0.01, 0.1)){
  for(custom_s in c(0.001, 0.01, 0.1)){
    c_exact_approx <- data.frame(c = seq(0,0.99, 0.01), exact = delta_c_dr(c=seq(0,0.99, 0.01), u = custom_u, s = custom_s, approx = FALSE), approx = delta_c_dr(c=seq(0,0.99, 0.01), u = custom_u, s = custom_s, approx = TRUE))
    plot_temp <- ggplot() +
    geom_line(data = c_exact_approx, aes(x = c, y = exact), linewidth = 1, linetype = "solid", color = "darkblue") +
    geom_line(data = c_exact_approx, aes(x = c, y = approx), linewidth = 1, linetype = "dashed", color = "darkred") + 
    labs(x = "Frequency (c)", y = "Change in frequency (\u0394 c)") +
     scale_color_manual(values = c("Exact" = "blue", "Approximation" = "red"),
                        name = "Genotypes",
                        labels = c("Aa + aa h=1" = "a", "a" = "Aa + aa h=1"),
                        guide = guide_legend(override.aes = list(linetype = c("solid", "dashed"))))+
#      theme_gray(base_size = 15) +
    theme(plot.title = element_text(size=11)) +
    theme(panel.background = element_rect(fill = 'white', color = 'black'),
          panel.grid.major = element_line(color = 'grey', linetype = 'dotted'),) + 
    labs(title = paste("s = ", custom_s, "; u = ", custom_u, sep = ""), tag = "b") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag = element_text()) +
    theme(legend.position = c(0.5, 0.25),
        legend.background = element_rect(fill = "white", color = "black")) +
    theme(aspect.ratio=3/3)
  
  # Add plot to the list
  plots[[length(plots) + 1]] <- plot_temp
    }
}
ggsave(filename = "plots/Supplementary_Figure_3.png", plot  = (plots[[1]] + plots[[2]] + plots[[3]]) / (plots[[4]] + plots[[5]] + plots[[6]]) / (plots[[7]] + plots[[8]] + plots[[9]]) + plot_annotation(tag_levels = "A"))







q()
