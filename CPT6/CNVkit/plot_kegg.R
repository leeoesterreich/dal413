library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)

# Read the data
gain_data <- read.delim("gain/KEGG_2021_Human.Human.enrichr.reports.txt")
loss_data <- read.delim("loss/KEGG_2021_Human.Human.enrichr.reports.txt")

# Process gain data
gain_data <- gain_data %>%
  mutate(log10_p_adj = -log10(Adjusted.P.value),
         Genes = gsub(";", "/", Genes)) %>%
  select(Term, log10_p_adj, Genes) %>%
  filter(log10_p_adj > 2) %>%  # Filter for significant pathways
  arrange(desc(log10_p_adj)) %>%
  slice_head(n = 10)  # Take top 10

# Process loss data
loss_data <- loss_data %>%
  mutate(log10_p_adj = -log10(Adjusted.P.value),
         Genes = gsub(";", "/", Genes)) %>%
  select(Term, log10_p_adj, Genes) %>%
  filter(log10_p_adj > 2) %>%  # Filter for significant pathways
  arrange(desc(log10_p_adj)) %>%
  slice_head(n = 10)  # Take top 10

# Find max x-axis limit for both plots
max_x <- max(c(gain_data$log10_p_adj, loss_data$log10_p_adj)) * 1.5

# Custom theme without gridlines
clean_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 8, margin = margin(r = 20)),  # Add margin for pathway names
    axis.text.x = element_text(size = 10),
    plot.margin = margin(5, 20, 5, 100),  # Adjusted margins
    panel.background = element_rect(fill = "white", color = NA),
    axis.line.x = element_line(color = "black"),  # Add x-axis line for better alignment
    plot.title = element_text(size = 12, hjust = 0, margin = margin(b = 10))
  )

# Function to create plot
create_plot <- function(data, type, color, show_xlab = FALSE) {
  # Calculate position for gene labels
  data$label_x <- data$log10_p_adj + 0.1
  
  p <- ggplot(data, aes(x = log10_p_adj, y = reorder(Term, log10_p_adj))) +
    geom_bar(stat = "identity", fill = color) +
    geom_text(aes(x = label_x, label = Genes), hjust = 0, size = 2, color = "gray40") +
    clean_theme +
    scale_x_continuous(limits = c(0, max_x), expand = c(0, 0), position = "bottom") +  # Remove expansion
    scale_y_discrete(labels = function(x) paste0("    ", x))  # Add space for labels
  
  # Add x-axis label only to bottom plot
  if(show_xlab) {
    p <- p + xlab("-log10(p.adjust)")
  } else {
    p <- p + xlab(NULL)
  }
  
  # Add type label with matching color
  p <- p + theme(plot.title = element_text(color = color)) +
    ggtitle(type)
  
  return(p)
}

# Create plots
p1 <- create_plot(gain_data, "Gain", "#FDA4AF", FALSE)
p2 <- create_plot(loss_data, "Loss", "#7DD3FC", TRUE)

# Combine plots with equal heights and aligned axes
pdf("kegg_enrichment.pdf", width = 12, height = 10)
grid.arrange(p1, p2, 
             heights = c(1, 1),
             padding = unit(-0.5, "lines"))  # Negative padding to ensure alignment
dev.off() 