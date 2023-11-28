# Create output directory
output_dir <- "Outputs/Beta_Diversity"
dir.create(output_dir, recursive = T, showWarnings = F)


# Create a relative abundance table
otu_relab <- sweep(otu_table, 2, colSums(otu_table), "/")
colSums(otu_relab)

# Bray-Curtis diversity
bc <- as.data.frame(as.matrix(vegan::vegdist(t(otu_relab), method = "bray")))

#### Analysis 1:

# Subset to baseline
bc_sub1 <- bc[meta_sub1$sample, meta_sub1$sample]

# PCOA
pco <- ape::pcoa(bc_sub1) #Another option:cmdscale
pco_df <- data.frame(PC1 = pco$vectors[, 1], PC2 = pco$vectors[, 2])
eig_vals <- round(pco$values[, "Relative_eig"][1:2] * 100, 2)
names(eig_vals) <- c("PC1", "PC2")

stopifnot(all.equal(rownames(pco_df), meta_sub1$sample))
pco_df <- cbind(pco_df, meta_sub1)

# PERMANOVA
fit <- vegan::adonis2(bc_sub1 ~ nationality + sex + bmi_group, data = meta_sub1, by = "margin")

# Plot
plot_beta <- ggplot(data = pco_df, aes(PC1, PC2, color = nationality)) +
  geom_point(alpha = 0.5, size = 3) +
  theme_classic() +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = paste0("PC1 ", "(", eig_vals["PC1"], "%", ")"), y = paste0("PC2 ", "(", eig_vals["PC2"], "%", ")"),
    title = "Bray Curtis Dissimilarity", caption = paste("PERMANOVA R2 =", round(fit$R2[1], 3), "p =", round(fit$`Pr(>F)`[1], 3))
  )

pdf(file.path(output_dir, "beta_div_nationality.pdf"), 7, 5)
print(plot_beta)
dev.off()

#### Analysis 2:

# Lists to save plots and dataframes
beta_plots <- list()
baseline_diff_list <- list()


for (nt in c("AAM", "AFR")) {
  # Subset to nationality
  meta_sub2 <- metadata %>% filter(nationality == nt)
  bc_sub2 <- bc[meta_sub2$sample, meta_sub2$sample]

  # PCOA
  pco <- ape::pcoa(bc_sub2)
  pco_df <- data.frame(PC1 = pco$vectors[, 1], PC2 = pco$vectors[, 2])
  eig_vals <- round(pco$values[, "Relative_eig"][1:2] * 100, 2)
  names(eig_vals) <- c("PC1", "PC2")

  # Add metadata
  stopifnot(all.equal(rownames(pco_df), meta_sub2$sample))
  pco_df <- cbind(pco_df, meta_sub2)

  # PERMANOVA
  fit <- vegan::adonis2(bc_sub2 ~ diet + sex + bmi_group, data = meta_sub2,by = "margin",strata =meta_sub2$subject )
  
  # Try LDM::permanovaFL()
  #fit_fl <- LDM::permanovaFL(bc_sub2 | (sex + bmi_group) ~ diet ,data = meta_sub2,cluster.id = "subject" )
  
  # Plot
  plot_beta <- ggplot(data = pco_df, aes(PC1, PC2, color = diet)) +
    geom_point(alpha = 0.7, size = 3) +
    theme_classic() +
    scale_color_manual(values = mycolors) +
    labs(
      x = paste0("PC1 ", "(", eig_vals["PC1"], "%", ")"), y = paste0("PC2 ", "(", eig_vals["PC2"], "%", ")"),
      title = nt, caption = paste("PERMANOVA R2 =", round(fit$R2[1], 3), "p =", round(fit$`Pr(>F)`[1], 3))
    )

  beta_plots[[nt]] <- plot_beta

  # Difference from the baseline

  # Loop through each subject and extract the bc dissimilarity of each timepoint from baseline
  for (id in unique(meta_sub2$subject)) {
    # Subset to subject
    meta_sub2_tmp <- meta_sub2 %>% filter(subject == id)

    # Row: baseline
    # Columns: Timepoints after baseline
    baseline_diff_tmp <- bc_sub2[meta_sub2_tmp$sample[meta_sub2_tmp$diet == "ED1"], meta_sub2_tmp$sample[meta_sub2_tmp$diet != "ED1"]]

    if (nrow(baseline_diff_tmp) != 0) {
      df_tmp <- data.frame(
        subject = id, diet = meta_sub2_tmp$diet[meta_sub2_tmp$diet != "ED1"],
        baseline_diff = as.numeric(baseline_diff_tmp[1, ]),
        nationality = nt
      )
      baseline_diff_list[[paste0(nt, "_", id)]] <- df_tmp
    } else {
      cat("No baseline sample for subject", id, "\n")
    }
  }
}

# Save PCOA plots
pdf(file.path(output_dir, "beta_diversity_diet.pdf"), 10, 5)
p <- ggpubr::ggarrange(plotlist = beta_plots, nrow = 1, common.legend = T, legend = "right")
annotate_figure(p, top = text_grob("Bray-Curtis Dissimilarity", face = "bold", size = 14))
dev.off()


# Convert the list of baseline differences to a dataframe
baseline_diff <- baseline_diff_list %>% bind_rows()

# Pairwise wilcoxon test for paired data
stat <- baseline_diff %>%
  group_by(subject) %>%
  filter(n() == 5) %>%
  ungroup() %>%
  group_by(nationality) %>%
  rstatix::pairwise_wilcox_test(baseline_diff ~ diet, paired = T) %>%
  add_xy_position()

baseline_diff_fig <- ggplot(baseline_diff, aes(x = diet, y = baseline_diff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.5, aes(color = diet)) +
  # geom_line(aes(group = subject), alpha = 0.3, color = "gray")+
  scale_color_manual(values = mycolors[-1]) +
  theme_classic() +
  facet_wrap(. ~ nationality) +
  labs(title = "Bray Curtis Dissimilariry from Baseline", y = "Dissimilarity from baseline") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_pvalue_manual(stat, hide.ns = T)

pdf(file.path(output_dir, "bc_from_baseline.pdf"), 7, 5)
print(baseline_diff_fig)
dev.off()

## Try CLME

# African
baseline_diff_afr <- baseline_diff %>% filter(nationality == "AFR")
clme_output_afr <- clme(baseline_diff ~ diet + (1|subject), data = baseline_diff_afr,
                    constraints = const, seed=42, nsim=100)
clme_summary_afr <- summary(clme_output_afr)
plot(clme_summary_afr, ci = TRUE, legendx = "bottomright", inset = 0.08)

# American
baseline_diff_aam <- baseline_diff %>% filter(nationality == "AAM")
clme_output_aam <- clme(baseline_diff ~ diet + (1|subject), data = baseline_diff_aam,
                        constraints = const, seed=42, nsim=100)
clme_summary_aam <- summary(clme_output_aam)
plot(clme_summary_aam, ci = TRUE, legendx = "bottomright", inset = 0.08)




