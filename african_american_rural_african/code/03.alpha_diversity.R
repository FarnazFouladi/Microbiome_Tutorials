# Create output directory
output_dir <- "Outputs/Alpha_Diversity"
dir.create(output_dir, recursive = T, showWarnings = F)

# Calculate alpha diversity metrics
set.seed(13654)
otu_rarefied <- vegan::rrarefy(t(otu_table), min(colSums(otu_table)))
metadata$Shannon_Index <- vegan::diversity(otu_rarefied, index = "shannon")
metadata$Richness <- rowSums(otu_rarefied > 0)
metadata$Chao1 <- estimateR(otu_rarefied)[2, ]

#### Analysis 1: AAM vs ARF at baseline

# Subset to baseline
meta_sub1 <- metadata %>% filter(diet == "ED1")

# Alpha diversity metrics
alpha_metrics <- c("Shannon_Index", "Richness", "Chao1")

# For each alpha diversity metric, perform a wilcoxon test and create a boxplot
alpha_plots1 <- lapply(alpha_metrics, function(alpha) {
  stat <- meta_sub1 %>%
    rstatix::wilcox_test(as.formula(paste0(alpha, "~ nationality"))) %>%
    add_xy_position()

  # plot
  alpha_fig <- ggplot(meta_sub1, aes(x = nationality, y = .data[[alpha]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.5, aes(color = nationality)) +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    stat_pvalue_manual(stat, label = "wilcoxon test p = {p}") +
    labs(title = alpha)
})

# Save plots
pdf(file.path(output_dir, "alpha_div_nationality.pdf"), 15, 5)
print(ggpubr::ggarrange(plotlist = alpha_plots1, common.legend = T, legend = "right", nrow = 1))
dev.off()

#### Analysis 2: Compare each timepoint to the baseline

alpha_plots2 <- lapply(alpha_metrics, function(alpha) {
  # Subset to paired data
  paired_metadata <- metadata %>%
    group_by(subject) %>%
    filter(n() == 6) %>%
    ungroup()

  stat <- paired_metadata %>%
    group_by(nationality) %>%
    rstatix::pairwise_wilcox_test(as.formula(paste0(alpha, "~ diet")), paired = T, ref.group = "ED1") %>%
    add_xy_position()


  alpha_fig <- ggplot(paired_metadata, aes(x = diet, y = .data[[alpha]])) +
    geom_violin(trim = TRUE, alpha = 0.5, aes(fill = diet)) +
    geom_jitter(position = position_jitter(width = 0.1), size = 1, alpha = 0.5, color = "black") +
    geom_boxplot(outlier.shape = NA, width = 0.1, alpha = 0.5) +
    geom_line(aes(group = subject), alpha = 0.3, color = "gray") +
    scale_fill_manual(values = mycolors)+
    theme_classic() +
    facet_wrap(. ~ nationality) +
    labs(title = alpha) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    stat_pvalue_manual(stat, hide.ns = T)
})


# Save plot
pdf(file.path(output_dir, "alpha_div_diet.pdf"), 7, 10)
print(ggpubr::ggarrange(plotlist = alpha_plots2, ncol = 1))
dev.off()
