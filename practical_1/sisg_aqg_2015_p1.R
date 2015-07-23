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


twin.data <- read.table("practical_1/data/twin_height_bmi.txt", header = T)

dim(twin.data)
head(twin.data)

twin.mz <- twin.data[twin.data$twin == "MZ", ]
cor(twin.mz $ht_t1, twin.mz $ht_t2)

p <- ggplot(twin.mz, aes(x = ht_t1, y = ht_t2))
p + geom_point() + geom_smooth(method = "lm")
# Dizygotic twins
twin.dz <- twin.data[twin.data$twin == "DZ", ]
cor(twin.dz$ht_t1, twin.dz$ht_t2)
p <- ggplot(twin.dz, aes(x = ht_t1, y = ht_t2))
p + geom_point() + geom_smooth(method = "lm")

ped <- read.table("practical_1/data/height_bmi.ped", header = F)
dim(ped)

names(ped) <- c("FAMILY", "PERSON", "FATHER", "MOTHER", "SEX")

colnames(ped)[6:11] <- c("TWIN_STATUS", "HEIGHT","NI_1", "BMI",
                          "NI_2", "NI_3")
head(ped)


dat <- read.table("practical_1/data/height_bmi.dat")

dim(dat)
head(dat)


install.packages("kinship2")
library("kinship2")

ped.sub <- ped[1:12, ]

p <- pedigree(ped.sub$PERSON, ped.sub$FATHER,
              ped.sub$MOTHER, ped.sub$SEX)


plot.pedigree(p)


# from practical_1_final_plot.R, sent by Maciej Trzaskowski, 2015/07/22
# ............................

twin.data <- read.table("practical_1/data/twin_height_bmi.txt", header = T)

twin.data$rel <- ifelse(twin.data$twin=="MZ",1,2)

mztmp=which(twin.data$twin=="MZ")
dztmp=which(twin.data$twin=="DZ")


rb=cor(twin.data$bmi_t1, twin.data$bmi_t2,use="pairwise.complete.obs")
rbmz=cor(twin.data$bmi_t1[mztmp], twin.data$bmi_t2[mztmp],use="pairwise.complete.obs")
rbdz=cor(twin.data$bmi_t1[dztmp], twin.data$bmi_t2[dztmp],use="pairwise.complete.obs")
(h2b=2*(rbmz-rbdz))
rh=cor(twin.data$ht_t1, twin.data$ht_t2,use="pairwise.complete.obs")
rhmz=cor(twin.data$ht_t1[mztmp], twin.data$ht_t2[mztmp],use="pairwise.complete.obs")
rhdz=cor(twin.data$ht_t1[dztmp], twin.data$ht_t2[dztmp],use="pairwise.complete.obs")
(h2h=2*(rhmz-rhdz))    

pdf('raw_twin_cor.pdf')      
par(mfrow=c(2,2))
plot(twin.data$ht_t1[mztmp],twin.data$ht_t2[mztmp],pch=20,col='firebrick',main="MZ",xlab="Twin1",ylab="Twin2")
legend('bottomright',legend=sprintf("Height r=%s",round(rhmz,2)))        
plot(twin.data$ht_t1[dztmp],twin.data$ht_t2[dztmp],pch=20,col='firebrick',main="DZ",xlab="Twin1",ylab="Twin2")
legend('bottomright',legend=sprintf("Height r=%s",round(rhdz,2))) 
plot(twin.data$bmi_t1[mztmp],twin.data$bmi_t2[mztmp],pch=20,col='firebrick',main="MZ",xlab="Twin1",ylab="Twin2")
legend('bottomright',legend=sprintf("BMI r=%s",round(rbmz,2)))        
plot(twin.data$bmi_t1[dztmp],twin.data$bmi_t2[dztmp],pch=20,col='firebrick',main="DZ",xlab="Twin1",ylab="Twin2")
legend('bottomright',legend=sprintf("BMI r=%s",round(rbdz,2))) 
dev.off()

mzmtmp=which(twin.data$twin=="MZ" & twin.data$sex==1)
dzmtmp=which(twin.data$twin=="DZ" & twin.data$sex==1)

mzftmp=which(twin.data$twin=="MZ" & twin.data$sex==1)
dzftmp=which(twin.data$twin=="DZ" & twin.data$sex==1)

twin.data$zht_t1 <- NA
twin.data$zht_t2 <- NA
twin.data$zbmi_t1 <- NA
twin.data$zbmi_t2 <- NA

twin.data$zht_t1[mzmtmp] = (twin.data$ht_t1[mzmtmp]-mean(twin.data$ht_t1[mzmtmp],na.rm=T)) / sd(twin.data$ht_t1[mzmtmp])
twin.data$zht_t1[mzftmp] = (twin.data$ht_t1[mzftmp]-mean(twin.data$ht_t1[mzftmp],na.rm=T)) / sd(twin.data$ht_t1[mzftmp])
twin.data$zht_t1[dzmtmp] = (twin.data$ht_t1[dzmtmp]-mean(twin.data$ht_t1[dzmtmp],na.rm=T)) / sd(twin.data$ht_t1[dzmtmp])
twin.data$zht_t1[dzftmp] = (twin.data$ht_t1[dzftmp]-mean(twin.data$ht_t1[dzftmp],na.rm=T)) / sd(twin.data$ht_t1[dzftmp])

twin.data$zht_t2[mzmtmp] = (twin.data$ht_t2[mzmtmp]-mean(twin.data$ht_t2[mzmtmp],na.rm=T)) / sd(twin.data$ht_t2[mzmtmp])
twin.data$zht_t2[mzftmp] = (twin.data$ht_t2[mzftmp]-mean(twin.data$ht_t2[mzftmp],na.rm=T)) / sd(twin.data$ht_t2[mzftmp])
twin.data$zht_t2[dzmtmp] = (twin.data$ht_t2[dzmtmp]-mean(twin.data$ht_t2[dzmtmp],na.rm=T)) / sd(twin.data$ht_t2[dzmtmp])
twin.data$zht_t2[dzftmp] = (twin.data$ht_t2[dzftmp]-mean(twin.data$ht_t2[dzftmp],na.rm=T)) / sd(twin.data$ht_t2[dzftmp])


twin.data$zbmi_t1[mzmtmp] = (twin.data$bmi_t1[mzmtmp]-mean(twin.data$bmi_t1[mzmtmp],na.rm=T)) / sd(twin.data$bmi_t1[mzmtmp])
twin.data$zbmi_t1[mzftmp] = (twin.data$bmi_t1[mzftmp]-mean(twin.data$bmi_t1[mzftmp],na.rm=T)) / sd(twin.data$bmi_t1[mzftmp])
twin.data$zbmi_t1[dzmtmp] = (twin.data$bmi_t1[dzmtmp]-mean(twin.data$bmi_t1[dzmtmp],na.rm=T)) / sd(twin.data$bmi_t1[dzmtmp])
twin.data$zbmi_t1[dzftmp] = (twin.data$bmi_t1[dzftmp]-mean(twin.data$bmi_t1[dzftmp],na.rm=T)) / sd(twin.data$bmi_t1[dzftmp])

twin.data$zbmi_t2[mzmtmp] = (twin.data$bmi_t2[mzmtmp]-mean(twin.data$bmi_t2[mzmtmp],na.rm=T)) / sd(twin.data$bmi_t2[mzmtmp])
twin.data$zbmi_t2[mzftmp] = (twin.data$bmi_t2[mzftmp]-mean(twin.data$bmi_t2[mzftmp],na.rm=T)) / sd(twin.data$bmi_t2[mzftmp])
twin.data$zbmi_t2[dzmtmp] = (twin.data$bmi_t2[dzmtmp]-mean(twin.data$bmi_t2[dzmtmp],na.rm=T)) / sd(twin.data$bmi_t2[dzmtmp])
twin.data$zbmi_t2[dzftmp] = (twin.data$bmi_t2[dzftmp]-mean(twin.data$bmi_t2[dzftmp],na.rm=T)) / sd(twin.data$bmi_t2[dzftmp])


rb=cor(twin.data$zbmi_t1, twin.data$zbmi_t2,use="pairwise.complete.obs")
rbmz=cor(twin.data$zbmi_t1[mztmp], twin.data$zbmi_t2[mztmp],use="pairwise.complete.obs")
rbdz=cor(twin.data$zbmi_t1[dztmp], twin.data$zbmi_t2[dztmp],use="pairwise.complete.obs")
(h2b=2*(rbmz-rbdz))
rh=cor(twin.data$zht_t1, twin.data$zht_t2,use="pairwise.complete.obs")
rhmz=cor(twin.data$zht_t1[mztmp], twin.data$zht_t2[mztmp],use="pairwise.complete.obs")
rhdz=cor(twin.data$zht_t1[dztmp], twin.data$zht_t2[dztmp],use="pairwise.complete.obs")
(h2h=2*(rhmz-rhdz))    

pdf('zscore_sex_twin_cor.pdf')
par(mfrow=c(2,2))
plot(twin.data$zht_t1[mztmp],twin.data$zht_t2[mztmp],pch=20,col='firebrick',main="MZ",xlab="Twin1",ylab="Twin2")
legend('bottomright',legend=sprintf("Height r=%s",round(rhmz,2)))        
plot(twin.data$zht_t1[dztmp],twin.data$zht_t2[dztmp],pch=20,col='firebrick',main="DZ",xlab="Twin1",ylab="Twin2")
legend('bottomright',legend=sprintf("Height r=%s",round(rhdz,2))) 
plot(twin.data$bmi_t1[mztmp],twin.data$bmi_t2[mztmp],pch=20,col='firebrick',main="MZ",xlab="Twin1",ylab="Twin2")
legend('bottomright',legend=sprintf("BMI r=%s",round(rbmz,2)))        
plot(twin.data$bmi_t1[dztmp],twin.data$bmi_t2[dztmp],pch=20,col='firebrick',main="DZ",xlab="Twin1",ylab="Twin2")
legend('bottomright',legend=sprintf("BMI r=%s",round(rbdz,2))) 
dev.off()

