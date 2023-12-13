# Compare different beta diversity metrics
# Compare NMDS versus MDS

# Output directory
output_dir <- "Outputs/Beta_Diversity"


library(mia)
# https://www.bioconductor.org/packages/devel/bioc/vignettes/mia/inst/doc/mia.html#merging-and-agglomeration-based-on-taxonomic-information.
data("GlobalPatterns", package = "mia")
GlobalPatterns
tse <- GlobalPatterns

# Retrieve data
otu_df <- as.data.frame(assay(tse), check.names = FALSE)
tax_df <- data.frame(rowData(tse), check.names = FALSE)
meta <- data.frame(colData(tse), check.names = FALSE)
tree_phy <- rowTree(tse)

meta$SampleID <- as.character(meta$X.SampleID)
all.equal(meta$SampleID, colnames(otu_df))
all.equal(rownames(tax_df), rownames(otu_df))

# Convert to relative abundance
otu_relab <- sweep(otu_df, 2, colSums(otu_df), "/")

# CLR transformation
otu_relab1 <- otu_relab
otu_relab1 <- otu_relab1 + min(otu_relab1[otu_relab1 != 0]) / 2
otu_clr <- vegan::decostand(otu_relab1, MARGIN = 2, method = "clr")


# colors
cols_sampletypes <- c("#FF7F00","#A6CEE3","#1F78B4","grey", "#023858","#B30000","#6A3D9A","#B15928","#FB9A99")

################
# Beta diversity
################

# Functions to perform MDS (Principal Coordinate Analysis)
# Difference between pcoa and pco? Distance matrix versus covariance matrix
pcoa_function <- function(dist_df, meta_df = meta, group = "SampleType", distance_metric = "Bray Curtis Dissimilarity") {
  # Principal coordinate analysis
  pco <- ape::pcoa(dist_df)
  pco_df <- data.frame(PC1 = pco$vectors[, 1], PC2 = pco$vectors[, 2])
  eig_vals <- round(pco$values[, "Relative_eig"][1:2] * 100, 2)
  names(eig_vals) <- c("PC1", "PC2")
  stopifnot(all.equal(rownames(pco_df), meta_df$SampleID))
  pco_df <- cbind(pco_df, meta_df)

  # Plot
  plot_beta <- ggplot(data = pco_df, aes(PC1, PC2, color = .data[[group]])) +
    geom_point(alpha = 0.8, size = 3) +
    theme_classic() +
    # scale_color_brewer(palette = "Dark2") +
    labs(
      x = paste0("PC1 ", "(", eig_vals["PC1"], "%", ")"), y = paste0("PC2 ", "(", eig_vals["PC2"], "%", ")"),
      title = distance_metric
    )+
    scale_color_manual(values = cols_sampletypes)

  return(list(plot_beta, pco))
}


# Function to perform NMDS
# No assumption about linear relationship
# An iterative process
# Arranges points to maximize rank-order correlation between real-world distance and ordination space distance
# Goodness of fit is measured by stress (the disagreement between observed distances and fitted distances in the reduced dimension)
#
nmds_function <- function(dist_df , meta_df = meta, group = "SampleType", distance_metric = "Bray Curtis Dissimilarity") {
  # nonmetric multidimensional scaling
  pco <- vegan::metaMDS(dist_df, k = 2, autotransform = F)
  pco_df <- data.frame(PC1 = pco$points[, 1], PC2 = pco$points[, 2])
  stopifnot(all.equal(rownames(pco_df), meta_df$SampleID))
  pco_df <- cbind(pco_df, meta_df)

  # Plot
  plot_beta <- ggplot(data = pco_df, aes(PC1, PC2, color = .data[[group]])) +
    geom_point(alpha = 0.8, size = 3) +
    theme_classic() +
    # scale_color_brewer(palette = "Dark2") +
    labs(
      x = "NMDS-1", y = "NMDS-2",
      title = distance_metric
    )+
    scale_color_manual(values = cols_sampletypes)

  return(list(plot_beta, pco))
}

# Goodness of fit for MDS
# Obtained from: https://microbiome.github.io/OMA/community-similarity.html#unsupervised-ordination
# For NMDS, monotonic regression is used (rank order)
calculate_stress <- function(d0 = as.matrix(bc), dp = pcoa_bc_dist) {
  # Calculate stress i.e. relative difference in the original and
  # projected dissimilarities
  stress <- sum((dp - d0)^2) / sum(d0^2)

  ord <- order(as.vector(d0))
  df <- data.frame(
    d0 = as.vector(d0)[ord],
    dmds = as.vector(dp)[ord]
  )

  plot <- ggplot(df, aes(x = d0, y = dmds)) +
    geom_smooth() +
    geom_point() +
    labs(
      title = "Shepard plot",
      x = "Original distance",
      y = "MDS distance",
      subtitle = paste("Stress:", round(stress, 2))
    ) +
    theme_bw()

  return(plot)
}


# 1. Bray-Curtis distance
bc <- as.data.frame(as.matrix(vegan::vegdist(t(otu_relab), method = "bray")))
pcoa_bc <- pcoa_function(bc, meta, "SampleType")
nmds_bc <- nmds_function(bc, meta, "SampleType")

pcoa_bc_p <- pcoa_bc[[1]]
nmds_bc_p <- nmds_bc[[1]]

pcoa_bc_df <- pcoa_bc[[2]]$vectors[, 1:2]
nmds_bc_df <- nmds_bc[[2]]$points[, 1:2]

pcoa_bc_dist <- as.matrix(dist(pcoa_bc_df))
nmds_bc_dist <- as.matrix(dist(nmds_bc_df))

calculate_stress(as.matrix(bc), pcoa_bc_dist)
calculate_stress(as.matrix(bc), nmds_bc_dist) # Note: This stress function is not used for nmds
stressplot(nmds_bc[[2]]) # R^ = 1-Stress^2
nmds_bc[[2]]$stress # 0.1660493



# 2. Aitchison distance
aitch <- as.data.frame(as.matrix(vegan::vegdist(t(otu_clr), method = "euclidean")))

pcoa_aitch <- pcoa_function(aitch, meta, "SampleType", "Aitchison Distance")
nmds_aitch <- nmds_function(aitch, meta, "SampleType", "Aitchison Distance")

pcoa_aitch_p <- pcoa_aitch[[1]]
nmds_aitch_p <- nmds_aitch[[1]]

pcoa_aitch_df <- pcoa_aitch[[2]]$vectors[, 1:2]
nmds_aitch_df <- nmds_aitch[[2]]$points[, 1:2]

pcoa_aitch_dist <- as.matrix(dist(pcoa_aitch_df))
nmds_aitch_dist <- as.matrix(dist(nmds_aitch_df))
calculate_stress(as.matrix(aitch), pcoa_aitch_dist)
calculate_stress(as.matrix(aitch), nmds_aitch_dist) # Note: This stress function is not used for nmds
stressplot(nmds_aitch[[2]]) # R^ = 1-Stress^2
nmds_aitch[[2]]$stress # 0.1027627



# 3. Weighted Unifrac
unifrac <- GUniFrac::GUniFrac(t(otu_relab), tree_phy, size.factor = NULL, alpha = c(0, 0.5, 1), verbose = TRUE)
w_unifrac <- unifrac$unifracs[, , "d_1"]
# uw_unifrac <- unifrac$unifracs[,,"d_0"] # Unweighted unifrac

pcoa_unifrac <- pcoa_function(w_unifrac, meta, "SampleType", "Unifrac Distance")
nmds_unifrac <- nmds_function(w_unifrac, meta, "SampleType", "Unifrac Distance")

pcoa_unifrac_p <- pcoa_unifrac[[1]]
nmds_unifrac_p <- nmds_unifrac[[1]]

pcoa_unifrac_df <- pcoa_unifrac[[2]]$vectors[, 1:2]
nmds_unifrac_df <- nmds_unifrac[[2]]$points[, 1:2]

pcoa_unifrac_dist <- as.matrix(dist(pcoa_unifrac_df))
nmds_unifrac_dist <- as.matrix(dist(nmds_unifrac_df))
calculate_stress(as.matrix(w_unifrac), pcoa_unifrac_dist)
calculate_stress(as.matrix(w_unifrac), nmds_unifrac_dist) # Note: This stress function is not used for nmds
stressplot(nmds_unifrac[[2]]) # R^ = 1-Stress^2
nmds_unifrac[[2]]$stress # 0.1244809


# 4. Rao's quadratic entropy (Explore)

library(picante)
ultra_tree <- phytools::force.ultrametric(tree_phy, method = "extend")
# rao_qa <- raoD(t(otu_df),ultra_tree)
rao_qa <- readRDS(file.path(output_dir, "raoD_outputs.rds"))
names(rao_qa)
dist_btw <- round(rao_qa$H, digits = 8)

pcoa_rao <- pcoa_function(dist_btw, meta, "SampleType", "Rao's based distance")
nmds_rao <- nmds_function(dist_btw, meta, "SampleType", "Rao's based distance")

pcoa_rao_p <- pcoa_rao[[1]]
nmds_rao_p <- nmds_rao[[1]]

pcoa_rao_df <- pcoa_rao[[2]]$vectors[, 1:2]
nmds_rao_df <- nmds_rao[[2]]$points[, 1:2]

pcoa_rao_dist <- as.matrix(dist(pcoa_rao_df))
nmds_rao_dist <- as.matrix(dist(nmds_rao_df))
calculate_stress(as.matrix(dist_btw), pcoa_rao_dist)
calculate_stress(as.matrix(dist_btw), nmds_rao_dist)
stressplot(nmds_rao[[2]]) # R^ = 1-Stress^2
nmds_rao[[2]]$stress # 0.1190026


# Compare distances
plot(as.vector(dist_btw), as.vector(w_unifrac), xlab = "Rao's quadratic entropy", ylab = "Weighted Unifrac")
plot(as.vector(dist_btw), as.vector(as.matrix(bc)), xlab = "Rao's quadratic entropy", ylab = "Bray-Curtis")
plot(as.vector(dist_btw), as.vector(as.matrix(aitch)), xlab = "Rao's quadratic entropy", ylab = "Aitchison's Distance")
plot(as.vector(as.matrix(bc)), as.vector(as.matrix(aitch)), xlab = "Bray Curtis Distance", ylab = "Aitchison's Distance")


# Save plots
pdf(file.path(output_dir, "pcoa_nmds_plots.pdf"), 15, 15)
ggarrange(pcoa_bc_p, nmds_bc_p,
  pcoa_aitch_p, nmds_aitch_p,
  pcoa_unifrac_p, nmds_unifrac_p,
  pcoa_rao_p, nmds_rao_p,
  nrow = 4, ncol = 2
)

dev.off()
