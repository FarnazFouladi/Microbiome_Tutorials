# Create output directory
output_dir <- "Outputs/Differential_Abundance"
dir.create(output_dir, recursive = T, showWarnings = F)

# Analysis 1

# ANCOM-BC2 Requires data to be a phyloseq TreeSummarizedExperiment
# See https://bioconductor.org/packages/release/bioc/vignettes/TreeSummarizedExperiment/inst/doc/Introduction_to_treeSummarizedExperiment.html
# https://www.bioconductor.org/packages/devel/bioc/vignettes/mia/inst/doc/mia.html

print(dietswap)

#### Analysis 1:

# Subset to baseline
phy_sub1 <- subset_samples(dietswap, diet == "ED1")

diff_out1 <- ancombc2(
  data = phy_sub1, assay_name = "counts", tax_level = "Genus",
  fix_formula = "nationality + sex + bmi_group", rand_formula = NULL,
  p_adj_method = "BH", pseudo_sens = TRUE,
  prv_cut = 0.15, lib_cut = 1000, s0_perc = 0.05,
  group = "nationality", struc_zero = TRUE, neg_lb = TRUE,
  alpha = 0.05, n_cl = 2, verbose = TRUE
)


# Bias corrected abundances
bias_corrected_df <- diff_out1$bias_correct_log_table

# Replace NAs with the minimum plus half of the minimum
bias_corrected_df <- apply(bias_corrected_df, 2, function(x) {
  x[is.na(x)] <- min(x, na.rm = T) + (min(x, na.rm = T) / 2)
  return(x)
})

# Structural zeros
zero_df <- diff_out1$zero_ind
sig_taxa_with_st_zero <- zero_df$taxon[zero_df %>%
  column_to_rownames("taxon") %>%
  rowSums(.) == 1]

# Extract results
res1 <- diff_out1$res

# Subset to significant taxa and then create new columns for plotting
res1_sig <- res1 %>%
  filter(diff_nationalityAFR) %>%
  arrange(lfc_nationalityAFR) %>%
  mutate(
    taxon = factor(taxon, levels = unique(taxon)),
    change = ifelse(lfc_nationalityAFR < 0, "negative_lfc", "positive_lfc"),
    p_sig = ifelse(q_nationalityAFR < 0.001, "***", ifelse(q_nationalityAFR < 0.01, "**", "*")),
    y_p = ifelse(lfc_nationalityAFR > 0, lfc_nationalityAFR + 0.2, lfc_nationalityAFR - 0.2),
    ss_test = ifelse(passed_ss_nationalityAFR, "aquamarine3", "black")
  )

# lfc Plot
res1_fig <- ggplot(res1_sig, aes(x = taxon, y = lfc_nationalityAFR, fill = change)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = res1_sig$ss_test,size = 5)) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_text(aes(y = y_p, label = p_sig)) +
  labs(x = "")

pdf(file.path(output_dir,"sig_lfc_nationality.pdf"),8,5)
print(res1_fig)
dev.off()


# Boxplots for significant taxa

# Merge the corrected bias abundances and merge with metadata
taxa_meta_merged <- t(bias_corrected_df) %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(meta(phy_sub1))

boxplots <- lapply(as.character(res1_sig$taxon), function(taxon) {
  # Create a dataframe to include p-values
  pval_tmp <- data.frame(
    group1 = "AAM", group2 = "AFR", p = signif(res1_sig[res1_sig$taxon == taxon, ]$q_nationalityAFR, 3),
    y.position = 0.3 + max(taxa_meta_merged[, as.character(taxon)])
  )


  plot <- ggplot(taxa_meta_merged, aes(y = .data[[as.character(taxon)]], x = .data[["nationality"]])) +
    geom_boxplot() +
    geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.5, aes(color = nationality)) +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    labs(title = taxon) +
    stat_pvalue_manual(pval_tmp, label = "p")
})

ggsave(
  filename = file.path(output_dir, "differential_abundance_nationality.pdf"),
  plot = marrangeGrob(boxplots, nrow = 1, ncol = 1),
  width = 5, height = 5
)

#### Analysis 2:

outputs_list <- list()

for (nt in c("AAM", "AFR")) {
  # Subset to nationality
  phy_sub <- subset_samples(dietswap, nationality == nt)

  diff_out2 <- ancombc2(
    data = phy_sub, assay_name = "counts", tax_level = "Genus",
    fix_formula = "diet + sex + bmi_group", rand_formula = "(1|subject)",
    p_adj_method = "BH", pseudo_sens = TRUE,
    prv_cut = 0.15, lib_cut = 1000, s0_perc = 0.05,
    group = "diet", struc_zero = TRUE, neg_lb = TRUE,
    global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
    alpha = 0.05, n_cl = 2, verbose = TRUE,
    trend_control = list(
      contrast = list(
        matrix(c(1, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, -1, 1), # Increasing trend
          nrow = 5,
          byrow = TRUE
        ),
        matrix(c(-1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1), # decreasing trend
          nrow = 5,
          byrow = TRUE
        )
      ),
      node = list(5, 5),
      solver = "ECOS",
      B = 10
    )
  )

  outputs_list[[nt]] <- diff_out2
}

## Extract results:

# Bias corrected abundance
bias_corrected_df <- outputs_list[["AAM"]]$bias_correct_log_table %>%
  rownames_to_column("taxon") %>%
  full_join(outputs_list[["AFR"]]$bias_correct_log_table %>% rownames_to_column("taxon")) %>%
  column_to_rownames("taxon")

# Replace NAs with minimum + minimum/2
bias_corrected_df <- apply(bias_corrected_df, 2, function(x) {
  x[is.na(x)] <- min(x, na.rm = T) + (min(x, na.rm = T) / 2)
  return(x)
})

# Primary results
res1 <- rbind(outputs_list[["AAM"]]$res, outputs_list[["AFR"]]$res) %>%
  mutate(nationality = c(rep("AAM", nrow(outputs_list[["AAM"]]$res)), rep("AFR", nrow(outputs_list[["AFR"]]$res))))

# Dunnet results
dunn_results <- rbind(outputs_list[["AAM"]]$res_dunn, outputs_list[["AFR"]]$res_dunn) %>%
  mutate(nationality = c(rep("AAM", nrow(outputs_list[["AAM"]]$res_dunn)), rep("AFR", nrow(outputs_list[["AFR"]]$res_dunn))))

# Trend results
trend_results <- rbind(outputs_list[["AAM"]]$res_trend, outputs_list[["AFR"]]$res_trend) %>%
  mutate(nationality = c(rep("AAM", nrow(outputs_list[["AAM"]]$res_trend)), rep("AFR", nrow(outputs_list[["AFR"]]$res_trend))))

# Subset to significant taxa and then create new columns for plotting
dunn_sig <- dunn_results %>%
  filter(diff_dietHE1 | diff_dietHE2 | diff_dietDI1 | diff_dietDI2 | diff_dietED2) %>%
  pivot_longer(
    cols = !c(taxon, nationality),
    cols_vary = "slowest",
    names_to = c(".value", "Diet"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  mutate(p_sig = ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", "ns")))) %>%
  mutate_if(is.numeric, function(x) round(x, 2)) %>%
  mutate(
    label = ifelse(p_sig != "ns", paste0(lfc, "\n", p_sig), lfc),
    Diet = factor(gsub("diet", "", Diet), levels = c("HE1", "HE2", "DI1", "DI2", "ED2"))
  )

# Heatmap of lfc
dunn_fig <- ggplot(dunn_sig, aes(y = taxon, x = Diet, fill = lfc)) +
  geom_tile(color = "black") +
  geom_text(aes(label = label), color = "black", size = 3) +
  labs(x = NULL, y = NULL, title = "Log fold changes from baseline") +
  theme_minimal() +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    na.value = "white", midpoint = 0, limit = c(-max(abs(dunn_sig$lfc)), max(abs(dunn_sig$lfc))),
    name = "LFC"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  facet_wrap(~nationality, scales = "free")

pdf(file.path(output_dir, "dunnett_heatmap.pdf"), 15, 10)
print(dunn_fig)
dev.off()

# Trend test

# Subset to significant taxa and then create new columns for plotting
trend_sig <- trend_results %>%
  filter(diff_abn) %>%
  pivot_longer(
    cols = lfc_dietHE1:se_dietED2,
    cols_vary = "slowest",
    names_to = c(".value", "Diet"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  mutate_if(is.numeric, function(x) round(x, 2)) %>%
  mutate(Diet = factor(gsub("diet", "", Diet), levels = c("HE1", "HE2", "DI1", "DI2", "ED2")))

trend_fig <- ggplot(trend_sig, aes(y = taxon, x = Diet, fill = lfc)) +
  geom_tile(color = "black") +
  geom_text(aes(label = lfc), color = "black", size = 3) +
  labs(x = NULL, y = NULL, title = "Significant taxa monotonically decreased or increased") +
  theme_minimal() +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    na.value = "white", midpoint = 0, limit = c(-max(abs(trend_sig$lfc)), max(abs(trend_sig$lfc))),
    name = "LFC"
  ) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  facet_wrap(~nationality, scales = "free")

pdf(file.path(output_dir, "trend_heatmap.pdf"), 15, 10)
print(trend_fig)
dev.off()


# Boxplots

# Merge the corrected bias abundances and merge with metadata
taxa_meta_merged <- t(bias_corrected_df) %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(meta(dietswap))

boxplots <- lapply(as.character(dunn_sig$taxon), function(taxon) {
  # Subset to taxon
  taxa_meta_sub <- data.frame(taxon = taxa_meta_merged[, taxon], diet = taxa_meta_merged$diet, nationality = taxa_meta_merged$nationality)
  dunn_sig_sub <- dunn_sig[dunn_sig$taxon == taxon, ] %>%
    mutate(group1 = "ED1", group2 = Diet)

  # Create a dataframe of pvalues
  pval_tmp <- taxa_meta_sub %>%
    group_by(nationality) %>%
    get_y_position(formula = taxon ~ diet, ref.group = "ED1") %>%
    left_join(dunn_sig_sub) %>%
    filter(!is.na(p_sig) & p_sig != "ns")


  plot <- ggplot(taxa_meta_sub, aes(y = taxon, x = diet)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.5, aes(color = diet)) +
    # geom_line(aes(group = subject), alpha = 0.5, color = "gray")+
    scale_color_manual(values = mycolors) +
    theme_classic() +
    labs(title = taxon, y = "bias-corrected log abundance") +
    facet_wrap(~nationality) +
    stat_pvalue_manual(pval_tmp, label = "p_sig", hide.ns = T)
})



ggsave(
  filename = file.path(output_dir, "dunnet_boxplots.pdf"),
  plot = marrangeGrob(boxplots, nrow = 1, ncol = 1),
  width = 10, height = 5
)
