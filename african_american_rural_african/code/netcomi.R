library(NetCoMi)
data("amgut1.filt")
data("amgut2.filt.phy")


net_spring <- netConstruct(amgut1.filt,
                           filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 50),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           measure = "spring",
                           measurePar = list(nlambda=10, 
                                             rep.num=10),
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 2,
                           seed = 123456)


props_spring <- netAnalyze(net_spring, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)
#?summary.microNetProps
summary(props_spring, numbNodes = 5L)


plotHeat(mat = props_spring$graphletLCC$gcm1,
         pmat = props_spring$graphletLCC$pAdjust1,
         type = "mixed",
         title = "GCM", 
         colorLim = c(-1, 1),
         mar = c(2, 0, 2, 0))

# Add rectangles highlighting the four types of orbits
graphics::rect(xleft   = c( 0.5,  1.5, 4.5,  7.5),
               ybottom = c(11.5,  7.5, 4.5,  0.5),
               xright  = c( 1.5,  4.5, 7.5, 11.5),
               ytop    = c(10.5, 10.5, 7.5,  4.5),
               lwd = 2, xpd = NA)

text(6, -0.2, xpd = NA, 
     "Significance codes:  ***: 0.001;  **: 0.01;  *: 0.05")


p <- plot(props_spring, 
          nodeColor = "cluster", 
          nodeSize = "eigenvector",
          title1 = "Network on OTU level with SPRING associations", 
          showTitle = TRUE,
          cexTitle = 2.3)

legend(0.7, 1.1, cex = 2.2, title = "estimated association:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)

p$q1$Arguments$cut

###################Spearman
net_pears <- netConstruct(amgut2.filt.phy,  
                          measure = "pearson",
                          normMethod = "clr",
                          zeroMethod = "multRepl",
                          sparsMethod = "threshold",
                          thresh = 0.3,
                          verbose = 3)

props_pears <- netAnalyze(net_pears, 
                          clustMethod = "cluster_fast_greedy")

plot(props_pears,
     edgeInvisFilter = "threshold",
     edgeInvisPar = 0.4,
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     repulsion = 0.8,
     rmSingles = TRUE,
     labelScale = FALSE,
     cexLabels = 1.6,
     nodeSizeSpread = 3,
     cexNodes = 2,
     hubBorderCol = "darkgray",
     title1 = "Network on OTU level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 2.3)

legend(0.7, 1.1, cex = 2.2, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"),
       bty = "n", horiz = TRUE)



net_pears_unsigned <- netConstruct(data = net_pears$assoEst1,
                                   dataType = "correlation", 
                                   sparsMethod = "threshold",
                                   thresh = 0.3,
                                   dissFunc = "unsigned",
                                   verbose = 3)

hist(net_pears$assoMat1, 100, xlim = c(-1, 1), ylim = c(0, 400),
     xlab = "Estimated correlation", 
     main = "Estimated correlations after sparsification")

hist(net_pears$adjaMat1, 100, ylim = c(0, 400),
     xlab = "Adjacency values", 
     main = "Transformed adjacencies (using the \"signed\" method)")

hist(net_pears_unsigned$adjaMat1, 100, ylim = c(0, 400),
     xlab = "Adjacency values", 
     main = "Transformed adjacencies (using the \"unsigned\" method)")


library(phyloseq)
data("amgut2.filt.phy")

# Agglomerate to genus level
amgut_genus <- tax_glom(amgut2.filt.phy, taxrank = "Rank6")

# Taxonomic table
taxtab <- as(tax_table(amgut_genus), "matrix")

# Rename taxonomic table and make Rank6 (genus) unique
amgut_genus_renamed <- renameTaxa(amgut_genus, 
                                  pat = "<name>", 
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Rank6")


# Network construction and analysis
net_genus <- netConstruct(amgut_genus_renamed,
                          taxRank = "Rank6",
                          measure = "pearson",
                          zeroMethod = "multRepl",
                          normMethod = "clr",
                          sparsMethod = "threshold",
                          thresh = 0.3,
                          verbose = 3)

props_genus <- netAnalyze(net_genus, clustMethod = "cluster_fast_greedy")

# Compute layout
graph3 <- igraph::graph_from_adjacency_matrix(net_genus$adjaMat1, 
                                              weighted = TRUE)
set.seed(123456)
lay_fr <- igraph::layout_with_fr(graph3)

# Row names of the layout matrix must match the node names
rownames(lay_fr) <- rownames(net_genus$adjaMat1)

plot(props_genus,
     layout = lay_fr,
     shortenLabels = "intelligent",
     labelLength = 10,
     labelPattern = c(5, "'", 3, "'", 3),
     nodeSize = "fix",
     nodeColor = "gray",
     cexNodes = 0.8,
     cexHubs = 1.1,
     cexLabels = 1.2,
     title1 = "Network on genus level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 2.3)

legend(0.7, 1.1, cex = 2.2, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)


# Compare two network:
# Split the phyloseq object into two groups
amgut_season_yes <- phyloseq::subset_samples(amgut2.filt.phy, 
                                             SEASONAL_ALLERGIES == "yes")
amgut_season_no <- phyloseq::subset_samples(amgut2.filt.phy, 
                                            SEASONAL_ALLERGIES == "no")

amgut_season_yes


n_yes <- phyloseq::nsamples(amgut_season_yes)

# Network construction
net_season <- netConstruct(data = amgut_season_no, 
                           data2 = amgut_season_yes,  
                           filtTax = "highestVar",
                           filtTaxPar = list(highestVar = 50),
                           filtSamp = "highestFreq",
                           filtSampPar = list(highestFreq = n_yes),
                           measure = "spring",
                           measurePar = list(nlambda=10, 
                                             rep.num=10),
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 2,
                           seed = 123456)

props_season <- netAnalyze(net_season, 
                           centrLCC = FALSE,
                           avDissIgnoreInf = TRUE,
                           sPathNorm = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = c("degree", "eigenvector"),
                           hubQuant = 0.9,
                           lnormFit = TRUE,
                           normDeg = FALSE,
                           normBetw = FALSE,
                           normClose = FALSE,
                           normEigen = FALSE)

summary(props_season)


plot(props_season, 
     sameLayout = FALSE, 
     nodeColor = "cluster",
     nodeSize = "mclr",
     labelScale = FALSE,
     cexNodes = 1.5, 
     cexLabels = 2.5,
     cexHubLabels = 3,
     cexTitle = 3.7,
     groupNames = c("No seasonal allergies", "Seasonal allergies"),
     hubBorderCol  = "gray40")

legend("bottom", title = "estimated association:", legend = c("+","-"), 
       col = c("#009900","red"), inset = 0.02, cex = 4, lty = 1, lwd = 4, 
       bty = "n", horiz = TRUE)


comp_season <- netCompare(props_season, 
                          permTest = FALSE, 
                          verbose = FALSE,
                          seed = 123456)

summary(comp_season, 
        groupNames = c("No allergies", "Allergies"),
        showCentr = c("degree", "between", "closeness"), 
        numbNodes = 5)


net_season_pears <- netConstruct(data = amgut_season_no, 
                                 data2 = amgut_season_yes, 
                                 filtTax = "highestVar",
                                 filtTaxPar = list(highestVar = 50),
                                 measure = "pearson", 
                                 normMethod = "clr",
                                 sparsMethod = "none", 
                                 thresh = 0.2,
                                 verbose = 3)

# Differential network construction
diff_season <- diffnet(net_season_pears,
                       diffMethod = "fisherTest", 
                       adjust = "lfdr")
