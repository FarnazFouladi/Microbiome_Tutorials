############## Does alpha diversity change overtime in AFR

# Subset to nationality
meta_sub2 <- metadata %>% dplyr::filter(nationality == "AFR")

# Mean at each timepoint
meta_mean <- meta_sub2 %>% group_by(diet) %>% summarise(mean = mean(Shannon_Index), nationality = unique(nationality))

plot <- ggplot(meta_sub2, aes(x = diet, y = Shannon_Index)) +
  geom_jitter(position = position_jitter(width = 0.1), size = 1, alpha = 0.5, color = "black") +
  geom_violin(trim = TRUE, alpha = 0.5) +
  geom_line(data = meta_mean, aes(x = diet,y = mean,group = nationality), color = "red",size = 1 )+
  geom_point(data = meta_mean, aes(x = diet,y = mean,group = nationality), color = "red",size = 3 )+
  scale_color_manual(values = mycolors) +
  theme_classic() +
  labs(title = "Shannon Index") +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )


###1. Linear Mixed Effect Models
lmer_out <- lmerTest::lmer(Shannon_Index ~ diet + sex + bmi_group + (1|subject),data = meta_sub2)
summary(lmer_out)
# Disadvantages:
# * No information about the order of timepoints/diet is used 
# * It compares only to baseline 


###2. Spline 
library(splinectomeR)
trend_results <- trendyspliner(data = meta_sub2 , xvar = "timepoint", yvar = "Shannon_Index",cases = "subject",category = "nationality",group = "AFR",
                               perms = 999, retain_perm = T)
#p-value = 0.002
trendyspliner.plot.perms(trend_results, xlabel = 'timepoint', ylabel = 'Shannon_Index')
# Disadvantages:
# * It does not specify the trend (increasing, decreasing,...)
# * It does not provide stats at different timepoints

###3. Linear Mixed Effects Models under Inequality Constraints
library(CLME)
const <- list(order = "simple", node = 1, decreasing = FALSE)
clme_output <- clme(Shannon_Index ~ diet + sex + bmi_group + (1|subject), data = meta_sub2,
             constraints = const, seed=42, nsim=100)
clme_summary <- summary(clme_output)
plot(clme_summary, ci = TRUE, legendx = "bottomright", inset = 0.08)
clme_output$theta


