# Running Spiec-Easi for the example dataset and Constructing networks
library(NetCoMi)

#### Analysis 1:
# Comparing the microbial networks between AAM and AFR

# Creating an output directory
output_dir <- "Outputs/Network_Analysis"
dir.create(file.path(output_dir))

# Subset to AAM and AFC

# Subset to each nationality
phy_sub_aam <- subset_samples(dietswap, nationality == "AAM" & diet == "ED1" )
phy_sub_afr <- subset_samples(dietswap, nationality == "AFR" & diet == "ED1")

print(phy_sub_aam)
print(phy_sub_afr)


# Sample size is different? If the difference is large, we might want to subsample the dataset with higher number of samples
n_sample <- nsamples(phy_sub_afr)

# Network construction
net_nationaity <- netConstruct(data = phy_sub_aam, 
                           data2 = phy_sub_afr,  
                           #filtTax = "highestVar", # Keep the 50 taxa with highest variance
                           #filtTaxPar = list(highestVar = 50),
                           #filtSamp = "totalReads",
                           #filtSampPar = list(filtSamp = 1000),
                           taxRank= "Genus",
                           measure = "spieceasi", 
                           measurePar = list(method='mb'), # lambda.min.ratio=1e-2, nlambda=15),
                           #measure = "spring",
                           #measurePar = list(nlambda=10, 
                           #                 rep.num=10),
                           #
                           #normMethod = "none", 
                           #zeroMethod = "none",
                           #sparsMethod = "none", 
                           dissFunc = "signed",
                           seed = 123456)


# Network comparison
props_nationality <- netAnalyze(net_nationaity, 
                           centrLCC = TRUE,
                           sPathNorm = FALSE, #the shortest path is the minimum sum of dissimilarities between two nodes.
                           clustMethod = "cluster_fast_greedy",
                           hubPar = c("degree", "eigenvector"),
                           hubQuant = 0.9)

summary(props_nationality)


pdf(file.path(output_dir,"spieceasi_network.pdf"),60,40)
plot(props_nationality, 
     sameLayout = TRUE, 
     repulsion = 0.95,
     shortenLabels = "simple",
     labelLength = 15,
     charToRm = "et rel.",
     layoutGroup = "union",
     rmSingles = "inboth", 
     nodeSize = "mclr", 
     #nodeSizeSpread = 3,
     nodeColor = "cluster",
     sameClustCol = TRUE,
     labelScale = FALSE,
     labelFont=0.8,
     cexNodes = 1.5, 
     cexLabels = 2.5,
     cexHubLabels = 3,
     cexTitle = 3.8,
     #nodeFilter = "clustMin", 
     #nodeFilterPar = 20, 
     #edgeInvisFilter = "threshold",
     #edgeInvisPar = 0.1,
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

# Differential association
diff_nationality <- diffnet(net_nationaity,
                       diffMethod = "discordant",#lfdrThresh=0.2,
                       #diffMethod = "fisherTest",
                       adjust = "lfdr")


plot(diffnet_nationaity)

# For unadjusted p-value
plot(diff_nationality, adjusted = FALSE,
     mar = c(2, 2, 5, 15), legendPos = c(1.2,1.2),
     labelFont=1.5,
     cexLabels = 0.4,
     labelScale = FALSE,
     legendArgs = list(bty = "n"),
     legendGroupnames = c("AAM", "AFR"),
     legendTitle = "Correlations:")


# Show only taxa with association difference more than 0.2
diff_mat <- diff_nationality$diffMat
diff_taxa <- rownames(diff_mat)[rowSums(diff_mat>0.2)>0]

pdf(file.path(output_dir,"spieceasi_network_diff.pdf"),60,40)
plot(props_nationality, 
     nodeFilter = "names",
     nodeFilterPar = diff_taxa,
     nodeColor = "gray",
     cexNodes = 1.5, 
     cexLabels = 3,
     cexTitle = 3.8,
     sameLayout = TRUE, 
     layoutGroup = "union",
     rmSingles = FALSE, 
     nodeSize = "clr",
     edgeTranspHigh = 20,
     labelScale = FALSE,
     groupNames = c("AAM", "AFR"),
     hubBorderCol  = "gray40")
dev.off()
     
  

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
     labelFont=1.5,
     cexLabels = 0.4,
     labelScale = FALSE,
     legendArgs = list(bty = "n"),
     legendGroupnames = c("AAM", "AFR"),
     legendTitle = "Correlations:")



assoPerm <- filematrix::fm.open(filenamebase = "Outputs/Correlation/assoPerm_comp" , readonly = TRUE)
dim(as.matrix(assoPerm))

filematrix::close("Outputs/Correlation/assoPerm_comp")
