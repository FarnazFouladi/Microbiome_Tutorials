# In this script, I will use linear models/linear mixed effects models for comparing 
# alpha diversity metrics between groups

# Create output directory
output_dir <- "Outputs/Alpha_Diversity_lm"
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

# For each alpha diversity metric, perform a linear model and create a boxplot
alpha_plots1 <- lapply(alpha_metrics, function(alpha) {
  fit <- summary(lm(as.formula(paste0(alpha, "~ nationality + sex + bmi_group")), data = meta_sub1))

  stat <- get_y_position(meta_sub1, as.formula(paste0(alpha, "~ nationality"))) %>%
    mutate(p = p_format(fit$coefficients[2, 4],3))

  # plot
  alpha_fig <- ggplot(meta_sub1, aes(x = nationality, y = .data[[alpha]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.5, aes(color = nationality)) +
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    stat_pvalue_manual(stat, label = "linear model p = {p}") +
    labs(title = alpha)
})

# Save plots
pdf(file.path(output_dir, "alpha_div_nationality_lm.pdf"), 15, 5)
print(ggpubr::ggarrange(plotlist = alpha_plots1, common.legend = T, legend = "right", nrow = 1))
dev.off()

#### Analysis 2: Compare each timepoints to the baseline

alpha_plots <- list()
for (nt in c("AAM", "AFR")) {
  # Subset to nationality
  meta_sub2 <- metadata %>% dplyr::filter(nationality == nt)

  alpha_plots2 <- lapply(alpha_metrics, function(alpha) {
    # Mixed linear model
    fit <- lmerTest::lmer(as.formula(paste0(alpha, "~ diet + sex + bmi_group + (1|subject)")), data = meta_sub2)
    fit_summ <- summary(fit)

    # Note:
    # We can do pairwise comparison using emmeans or glht
    # emmeans::emmeans(fit,specs = pairwise ~ diet )
    # summary(multcomp::glht(fit, linfct = mcp(diet = 'Tukey')))
    # summary(multcomp::glht(fit, linfct = mcp(diet= matrix( c(0,0,-1,1,0,0),1,dimnames = list("HE2 - DI1")))))
    # summary(multcomp::glht(fit, linfct = matrix(c(0,0,-1,1,0,0,0,0,0),1)))


    # Create a dataframe including stats results for labeling the figure
    stat <- meta_sub2 %>%
      get_y_position(as.formula(paste0(alpha, "~ diet")), ref.group = "ED1") %>%
      mutate(
        p = fit_summ$coefficients[2:6, "Pr(>|t|)"], p.adj = p.adjust(p, "BH"),
        p.adj.signif = ifelse(p.adj < 0.001, "***", ifelse(p.adj < 0.01, "**", ifelse(p.adj < 0.05, "*", "ns")))
      )

    # plot
    alpha_fig <- ggplot(meta_sub2, aes(x = diet, y = .data[[alpha]])) +
      geom_violin(trim = TRUE, alpha = 0.5, aes(fill = diet)) +
      geom_jitter(position = position_jitter(width = 0.1), size = 1, alpha = 0.5, color = "black") +
      geom_boxplot(outlier.shape = NA, width = 0.1, alpha = 0.5) +
      geom_line(aes(group = subject), alpha = 0.3, color = "gray") +
      scale_fill_manual(values = mycolors) +
      theme_classic() +
      facet_wrap(. ~ nationality) +
      labs(title = alpha) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "none"
      ) +
      stat_pvalue_manual(stat, hide.ns = T)
  })

  alpha_plots[[nt]] <- alpha_plots2
}

# Save plot
pdf(file.path(output_dir, "alpha_div_diet_lm.pdf"), 7, 10)
print(ggpubr::ggarrange(alpha_plots[[1]][[1]], alpha_plots[[2]][[1]],
  alpha_plots[[1]][[2]], alpha_plots[[2]][[2]],
  alpha_plots[[1]][[3]], alpha_plots[[2]][[3]],
  ncol = 2, nrow = 3
))
dev.off()

# Additional note:
# Time as a numeric value:
# Not able to do pairwise comparisons, however, with spline function we can see if the slope changes at a specific timepoint (knots)
# https://cran.r-project.org/web/packages/lspline/vignettes/lspline.html
# Example:
alpha = "Shannon_Index"
fit <- lmerTest::lmer(Shannon_Index ~ lspline::lspline(timepoint,knots = c(4), marginal = T) + sex + bmi_group + (1|subject),data = meta_sub2)

# Another option: splinectomeR
# https://www.frontiersin.org/articles/10.3389/fmicb.2018.00785/full
# It is suitable when we want to see if two groups (for example, treatment vs. control) differ significantly overtime
# Example
library(splinectomeR)
permu_result <- permuspliner(data = metadata, xvar = "timepoint", yvar = "Shannon_Index",cases = "subject",category = "nationality",
                             perms = 999, retain_perm = T)
permuspliner.plot.permdistance(permu_result, xlabel="timepoint")
permuspliner.plot.permsplines(permu_result, xvar="timepoint", yvar="Shannon_Index")

# Trend over time for each group
trend_results <- trendyspliner(data = metadata, xvar = "timepoint", yvar = "Shannon_Index",cases = "subject",category = "nationality",group = "AAM",
                               perms = 999, retain_perm = T,mean_center=F)
trendyspliner.plot.perms(trend_results, xlabel = 'timepoint', ylabel = 'Shannon_Index')
trend_results <- trendyspliner(data = metadata, xvar = "timepoint", yvar = "Shannon_Index",cases = "subject",category = "nationality",group = "AFR",
                               perms = 999, retain_perm = T,mean_center=F)
trendyspliner.plot.perms(trend_results, xlabel = 'timepoint', ylabel = 'Shannon_Index')


