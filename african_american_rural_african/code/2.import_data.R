# Setting working directory
#setwd("~/Documents/codes/african_american_rural_african")

# Get data
# See https://joey711.github.io/phyloseq/
data(dietswap, package = "microbiome")
print(dietswap)

# Relevant publication for this data:
#https://www.nature.com/articles/ncomms7342

# Add diet group to the metadata
sample_data(dietswap)$diet <- paste0(sample_data(dietswap)$group, sample_data(dietswap)$timepoint.within.group)
sample_data(dietswap)$diet <- factor(sample_data(dietswap)$diet,
  levels = c("ED1", "HE1", "HE2", "DI1", "DI2", "ED2")
)

# Retrieve data from the phyloseq object
otu_table <- as.data.frame(abundances(dietswap))
tax_table <- as.data.frame(tax_table(dietswap))
metadata <- meta(dietswap)

# Sanity check
all.equal(colnames(otu_table),metadata$sample)
all.equal(rownames(otu_table),rownames(tax_table))

# Create a phylum level
phylum_table <- otu_table %>%
  mutate(phylum = tax_table$Phylum) %>%
  remove_rownames() %>%
  group_by(phylum) %>%
  summarise_all(list(sum)) %>%
  column_to_rownames("phylum")

# Alternative method
phylum_phy <- aggregate_taxa(dietswap, level = "Phylum")

all.equal(phylum_table,as.data.frame(abundances(phylum_phy)))

# Check out number of reads
read_number <- colSums(otu_table)
hist(read_number)
range(read_number)

# Alternative method
read_number <- readcount(dietswap)

# Set colors
mycolors <- c("#919191", "#FFA000", "#FF6F00", "#C5E1A5", "#7CB342", "#1B5E20")
