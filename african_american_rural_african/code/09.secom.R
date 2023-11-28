# Create output directory
output_dir <- "Outputs/Correlation"
dir.create(output_dir, recursive = T, showWarnings = F)


#### Analysis 1:

plot_list <- list()

for (nt in c("AAM", "AFR")) {
  # Subset to nationality
  phy_sub <- subset_samples(dietswap, nationality == nt & Diet == "ED1")

  set.seed(123)
  ## Linear relationships
  res_linear <- secom_linear(
    data = list(phy_sub), assay_name = "counts",
    tax_level = "Phylum", pseudo = 0,
    prv_cut = 0.5, lib_cut = 1000, corr_cut = 0.5,
    wins_quant = c(0.05, 0.95), method = "pearson",
    soft = FALSE, thresh_len = 20, n_cv = 10,
    thresh_hard = 0.3, max_p = 0.005, n_cl = 2
  )

  ## Nonlinear relationships
  res_dist <- secom_dist(
    data = list(phy_sub), assay_name = "counts",
    tax_level = "Phylum", pseudo = 0,
    prv_cut = 0.5, lib_cut = 1000, corr_cut = 0.5,
    wins_quant = c(0.05, 0.95), R = 1000,
    thresh_hard = 0.3, max_p = 0.005, n_cl = 2
  )



  ## Plotting for linear relationship
  corr_linear <- res_linear$corr_fl
  cooccur_linear <- res_linear$mat_cooccur

  # Filter by co-occurrence
  overlap <- 10
  corr_linear[cooccur_linear < overlap] <- 0

  # Convert to a long dataframe
  df_linear <- data.frame(get_upper_tri(corr_linear)) %>%
    rownames_to_column("var1") %>%
    pivot_longer(cols = -var1, names_to = "var2", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(value = round(value, 2))

  tax_name <- sort(union(df_linear$var1, df_linear$var2))
  df_linear$var1 <- factor(df_linear$var1, levels = tax_name)
  df_linear$var2 <- factor(df_linear$var2, levels = tax_name)

  # Plot
  heat_linear_fl <- df_linear %>%
    ggplot(aes(var2, var1, fill = value)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white", na.value = "grey",
      midpoint = 0, limit = c(-1, 1), space = "Lab",
      name = NULL
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    geom_text(aes(var2, var1, label = value), color = "black", size = 4) +
    labs(x = NULL, y = NULL, title = nt, subtitle = "Pearson (Filtering)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 45, vjust = 1, size = 12, hjust = 1,
        face = "italic"
      ),
      axis.text.y = element_text(size = 12, face = "italic"),
      strip.text.x = element_text(size = 14),
      strip.text.y = element_text(size = 14),
      legend.text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 15),
      panel.grid.major = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    coord_fixed()

  plot_list[[paste0(nt, "_linear")]] <- heat_linear_fl

  ## Plotting for non-linear relationship
  corr_dist <- res_dist$dcorr_fl
  cooccur_dist <- res_dist$mat_cooccur

  # Filter by co-occurrence
  overlap <- 10
  corr_dist[cooccur_dist < overlap] <- 0

  df_dist <- data.frame(get_upper_tri(corr_dist)) %>%
    rownames_to_column("var1") %>%
    pivot_longer(cols = -var1, names_to = "var2", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(value = round(value, 2))

  tax_name <- sort(union(df_dist$var1, df_dist$var2))
  df_dist$var1 <- factor(df_dist$var1, levels = tax_name)
  df_dist$var2 <- factor(df_dist$var2, levels = tax_name)

  heat_dist_fl <- df_dist %>%
    ggplot(aes(var2, var1, fill = value)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white", na.value = "grey",
      midpoint = 0, limit = c(-1, 1), space = "Lab",
      name = NULL
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    geom_text(aes(var2, var1, label = value), color = "black", size = 4) +
    labs(x = NULL, y = NULL, title = nt, subtitle = "Distance (Filtering)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 45, vjust = 1, size = 12, hjust = 1,
        face = "italic"
      ),
      axis.text.y = element_text(size = 12, face = "italic"),
      strip.text.x = element_text(size = 14),
      strip.text.y = element_text(size = 14),
      legend.text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 15),
      panel.grid.major = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    ) +
    coord_fixed()

  plot_list[[paste0(nt, "_non_linear")]] <- heat_dist_fl
}

# Save plot
pdf(file.path(output_dir, "nationality_secom.pdf"), 10, 10)
ggpubr::ggarrange(plotlist = plot_list, nrow = 2, ncol = 2)
dev.off()



#### Analysis 2:

plot_list <- list()

for (nt in c("AAM", "AFR")) {
  # Subset to nationality
  phy_sub <- subset_samples(dietswap, nationality == nt)

  # subset to different timepoints:
  phy_sub_ED1 <- subset_samples(phy_sub, diet == "ED1")
  phy_sub_HE1 <- subset_samples(phy_sub, diet == "HE1")
  phy_sub_HE2 <- subset_samples(phy_sub, diet == "HE2")
  phy_sub_DI1 <- subset_samples(phy_sub, diet == "DI1")
  phy_sub_DI2 <- subset_samples(phy_sub, diet == "DI2")
  phy_sub_ED2 <- subset_samples(phy_sub, diet == "ED2")

  # Change the sample names to subject names. There should be overalpping in names across dataframes
  sample_names(phy_sub_ED1) <- meta(phy_sub_ED1)$subject
  sample_names(phy_sub_HE1) <- meta(phy_sub_HE1)$subject
  sample_names(phy_sub_HE2) <- meta(phy_sub_HE2)$subject
  sample_names(phy_sub_DI1) <- meta(phy_sub_DI1)$subject
  sample_names(phy_sub_DI2) <- meta(phy_sub_DI2)$subject
  sample_names(phy_sub_ED2) <- meta(phy_sub_ED2)$subject

  # Create TreeSummarizedExperiment objects
  tse1 <- makeTreeSummarizedExperimentFromPhyloseq(phy_sub_ED1)
  tse2 <- makeTreeSummarizedExperimentFromPhyloseq(phy_sub_HE1)
  tse3 <- makeTreeSummarizedExperimentFromPhyloseq(phy_sub_HE2)
  tse4 <- makeTreeSummarizedExperimentFromPhyloseq(phy_sub_DI1)
  tse5 <- makeTreeSummarizedExperimentFromPhyloseq(phy_sub_DI2)
  tse6 <- makeTreeSummarizedExperimentFromPhyloseq(phy_sub_ED2)

  set.seed(123)
  ## Linear relationships
  res_linear <- secom_linear(
    data = list(
      ED1 = tse1, HE1 = tse2, HE2 = tse3,
      DI1 = tse4, DI2 = tse5, ED2 = tse6
    ),
    assay_name = c(rep("counts", 6)),
    tax_level = c(rep("Phylum", 6)), pseudo = 0,
    prv_cut = 0.5, lib_cut = 1000, corr_cut = 0.5,
    wins_quant = c(0.05, 0.95), method = "pearson",
    soft = FALSE, thresh_len = 20, n_cv = 10,
    thresh_hard = 0.3, max_p = 0.005, n_cl = 2
  )


  ## Nonlinear relationships
  res_dist <- secom_dist(
    data = list(
      ED1 = tse1, HE1 = tse2, HE2 = tse3,
      DI1 = tse4, DI2 = tse5, ED2 = tse6
    ),
    assay_name = c(rep("counts", 6)),
    tax_level = c(rep("Phylum", 6)),
    pseudo = 0,
    prv_cut = 0.2, lib_cut = 1000, corr_cut = 0.5,
    wins_quant = c(0.05, 0.95), R = 1000,
    thresh_hard = 0.3, max_p = 0.005, n_cl = 2
  )



  # Loop through each measurement and create a heatmap
  for (measure in c("Pearson", "Distance")) {
    if (measure == "Pearson") {
      corr_mat <- res_linear$corr_fl
      cooccur_mat <- res_linear$mat_cooccur
    } else {
      corr_mat <- res_dist$dcorr_fl
      cooccur_mat <- res_dist$mat_cooccur
    }

    # Filter by co-occurrence
    overlap <- 10
    corr_mat[cooccur_mat < overlap] <- 0

    # Change the order of taxa based on timepoints
    tax_name <- colnames(corr_mat)
    tax_name <- tax_name[c(grep("ED1", tax_name), grep("HE1", tax_name), grep("HE2", tax_name), grep("DI1", tax_name), grep("DI2", tax_name), grep("ED2", tax_name))]
    corr_mat <- corr_mat[tax_name, tax_name]

    # Convert to long format
    df <- data.frame(get_upper_tri(corr_mat)) %>%
      rownames_to_column("var1") %>%
      pivot_longer(cols = -var1, names_to = "var2", values_to = "value") %>%
      filter(!is.na(value)) %>%
      mutate(
        var2 = gsub("\\...", " - ", var2),
        value = round(value, 2)
      )

    df$var1 <- factor(df$var1, levels = tax_name)
    df$var2 <- factor(df$var2, levels = tax_name)
    colors <- RColorBrewer::brewer.pal(6, "Dark2")
    names(colors) <- c("ED1", "HE1", "HE2", "DI1", "DI2", "ED2")
    txt_color <- sapply(stringi::stri_extract(tax_name, regex = "[A-Z]*[1-9]"), function(x) as.character(colors[x]))

    # Heatmap
    heat_fl <- df %>%
      ggplot(aes(x = var1, y = var2, fill = value)) +
      geom_tile(color = "black") +
      scale_fill_gradient2(
        low = "blue", high = "red", mid = "white",
        na.value = "grey", midpoint = 0, limit = c(-1, 1),
        space = "Lab", name = NULL
      ) +
      scale_x_discrete(drop = FALSE) +
      scale_y_discrete(drop = FALSE) +
      geom_text(aes(x = var1, y = var2, label = value), color = "black", size = 2) +
      labs(x = NULL, y = NULL, title = nt, subtitle = ifelse(measure == "Pearson", "Pearson (Filtering)", "Distance (Filtering)")) +
      theme_bw() +
      geom_vline(xintercept = 6.5, color = "blue", linetype = "dashed") +
      geom_hline(yintercept = 6.5, color = "blue", linetype = "dashed") +
      geom_vline(xintercept = 12.5, color = "blue", linetype = "dashed") +
      geom_hline(yintercept = 12.5, color = "blue", linetype = "dashed") +
      geom_vline(xintercept = 18.5, color = "blue", linetype = "dashed") +
      geom_hline(yintercept = 18.5, color = "blue", linetype = "dashed") +
      geom_vline(xintercept = 24.5, color = "blue", linetype = "dashed") +
      geom_hline(yintercept = 24.5, color = "blue", linetype = "dashed") +
      geom_vline(xintercept = 30.5, color = "blue", linetype = "dashed") +
      geom_hline(yintercept = 30.5, color = "blue", linetype = "dashed") +
      theme(
        axis.text.x = element_text(
          angle = 45, vjust = 1, size = 12, hjust = 1,
          face = "italic", color = txt_color
        ),
        axis.text.y = element_text(
          size = 12, face = "italic",
          color = txt_color
        ),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 15),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none"
      ) +
      coord_fixed()

    plot_list[[paste0(nt, "_", measure)]] <- heat_fl
  }
}

# Save plot
pdf(file.path(output_dir, "AAM_secom.pdf"), 20, 20)
ggpubr::ggarrange(plot_list[[1]], plot_list[[2]], nrow = 1, ncol = 2)
dev.off()



pdf(file.path(output_dir, "AFR_secom.pdf"), 20, 20)
ggpubr::ggarrange(plot_list[[3]], plot_list[[4]], nrow = 1, ncol = 2)
dev.off()
