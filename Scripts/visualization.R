# Step 7: Visualizing and interpreting results in R
# Import the ggplot2 library for plotting
library(ggplot2)

# Read the annotated variant data
variants <- read.table("b1_annotated.hg19_multianno.txt", header=TRUE, sep="\t")

# Generate a bar plot to visualize variant distribution by functional category
ggplot(variants, aes(x=Func.refGene)) +
  geom_bar() +
  theme_minimal() +
  xlab("Gene Function") +
  ylab("Number of Variants") +
  ggtitle("Distribution of Variants by Gene Function")
