#### Analysis 1:

# Subset to AAM and AFC

# Subset to nationality
phy_sub_aam <- subset_samples(dietswap, nationality == "AAM" & diet == "ED1" )
phy_sub_afr <- subset_samples(dietswap, nationality == "AFR" & diet == "ED1")

print(phy_sub_aam)
print(phy_sub_afr)


# Sample size is different? If the difference is large, we might want to subsample the dataset with higher number of samples
n_sample <- nsamples(phy_sub_afr)

# Network construction
net_nationaity <- netConstruct(data = phy_sub_aam, 
                           data2 = phy_sub_afr,  
                           filtTax = "highestVar",
                           filtTaxPar = list(highestVar = 50),
                           #filtSamp = "highestFreq",
                           #filtSampPar = list(highestFreq = n_sample),
                           taxRank= "Genus",
                           measure = "spieceasi", #methods: "pearson", "spearman", "bicor", "sparcc", "cclasso", "ccrepe", "spieceasi" (default), "spring", "gcoda" and "propr" 
                           measurePar = list(method='mb', lambda.min.ratio=1e-2, nlambda=15),
                           #normMethod = "none", 
                           #zeroMethod = "none",
                           #sparsMethod = "none", 
                           dissFunc = "signed",
                           seed = 123456)


# Network comparison
props_nationality <- netAnalyze(net_nationaity, 
                           centrLCC = TRUE,
                           sPathNorm = FALSE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = c("degree", "eigenvector"),
                           hubQuant = 0.9)

summary(props_nationality)


pdf(file.path(output_dir,"spieceasi_network.pdf"),50,40)
plot(props_nationality, 
     sameLayout = TRUE, 
     repulsion = 0.95,
     layoutGroup = "union",
     rmSingles = "inboth", 
     nodeSize = "mclr", 
     labelScale = FALSE,
     cexNodes = 1.5, 
     cexLabels = 2.5,
     cexHubLabels = 3,
     cexTitle = 3.8,
     groupNames = c("AAM", "AFR"),
     hubBorderCol  = "gray40")

legend("bottom", title = "estimated association:", legend = c("+","-"), 
       col = c("#009900","red"), inset = 0.02, cex = 4, lty = 1, lwd = 4, 
       bty = "n", horiz = TRUE)

dev.off()


# Quantitative network comparison without permutation
comp_nationality <- netCompare(props_nationality, 
                          permTest = FALSE, 
                          seed = 123456)

summary(comp_nationality, 
        groupNames = c("AAM", "AFR"),
        showCentr = c("degree", "between", "closeness"), 
        numbNodes = 5)

# Differential network construction
diff_nationality <- diffnet(net_nationaity,
                       diffMethod = "fisherTest",alpha = 0.9,
                       adjust = "lfdr")


plot(diffnet_nationaity)

# For unadjusted p-value
plot(diff_nationality, adjusted = FALSE,
     mar = c(2, 2, 5, 15), legendPos = c(1.2,1.2),
     legendArgs = list(bty = "n"),
     legendGroupnames = c("AAM", "AFR"),
     legendTitle = "Correlations:")











# Quantitative network comparison with permutation
comp_nationality <- netCompare(props_nationality, 
                          permTest = TRUE,
                          nPerm = 5, 
                          storeAssoPerm = TRUE,
                          fileStoreAssoPerm = "Outputs/Correlation/assoPerm_comp",
                          storeCountsPerm = FALSE, 
                          seed = 123456)

summary(comp_nationality)



# Construct differential network
diffnet_nationaity <- diffnet(net_nationaity, diffMethod = "permute", nPerm = 5, 
                         fileLoadAssoPerm = "Outputs/Correlation/assoPerm_comp",alpha = 0.5,
                         storeCountsPerm = FALSE)
plot(diffnet_nationaity)

# For unadjusted p-value
plot(diffnet_nationaity, adjusted = FALSE,
     mar = c(2, 2, 5, 15), legendPos = c(1.2,1.2),
     legendArgs = list(bty = "n"),
     legendGroupnames = c("AAM", "AFR"),
     legendTitle = "Correlations:")



assoPerm <- filematrix::fm.open(filenamebase = "Outputs/Correlation/assoPerm_comp" , readonly = TRUE)
dim(as.matrix(assoPerm))

filematrix::close("Outputs/Correlation/assoPerm_comp")
