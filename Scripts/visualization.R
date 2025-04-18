# Step 7: Visualizing and interpreting results in R
library(ggplot2)

variants <- read.table("b1_annotated.hg19_multianno.txt", header=TRUE, sep="\t")

ggplot(variants, aes(x=Func.refGene)) +
  geom_bar() +
  theme_minimal() +
  xlab("Gene Function") +
  ylab("Number of Variants") +
  ggtitle("Distribution of Variants by Gene Function")
