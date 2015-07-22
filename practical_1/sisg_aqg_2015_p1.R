# ===============================================================
# Summer Institute in Statistical Genetics, 2015, Seattle
# Advanced Quantitative Genetics
# Practical 1 
# Date:
# Author:
# ===============================================================
#setwd("~/Desktop/SISG_AQG_2015") # Mac and Linux

setwd("C:/Users/vonholle/Documents/grad_school/misc/training/sisg-2015/module23/sisg-2015-module23") # Windows

args(rnorm)
height <- rnorm(1000, 170, 2.8)
mean(height)
sd(height)
plot(density(height, adjust = 3), main = '')

y.means <- c(5, 15, 20)
y.sd <- 5

# Generate a SNP. Sample from the binomial distribution with success probability 0.4
# Success probability corresponds to the MAF
snp <- rbinom(1000, 2, 0.4)
# Generate Y given the SNP and bind to data frame
y <- rnorm(1000, y.means[factor(snp)], y.sd)
snp.data <- data.frame(cbind(y, snp))


# Plot the phenotype as a factored boxplot of the underlying genotype counts
# Load ggplot2
library(ggplot2)
p <- ggplot(snp.data, aes(x = factor(snp), y = y))
p + geom_boxplot() + xlab("Reference allele count")
