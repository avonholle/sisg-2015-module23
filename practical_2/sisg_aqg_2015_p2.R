# ===============================================================
# Summer Institute in Statistical Genetics, 2015, Seattle
# Advanced Quantitative Genetics
# Practical 2 
# Date:
# Author:
# ===============================================================
# setwd("~/Desktop/SISG_AQG_2015") # Mac and Linux
# setwd("c:/path/to/my/directory")

setwd("C:/Users/vonholle/Documents/grad_school/misc/training/sisg-2015/module23/sisg-2015-module23") # Windows


# Note: I borrowed the next few lines of code from howrigan_snp.cleaning.pdf (author info: https://www.atgu.mgh.harvard.edu/people/daniel_howrigan)

# Once you've set up your PLINK directory for this analysis and have set R to that directory, do this:
system("./plink",intern=TRUE)
# If is says PLINK! in the first few lines, you're ready to go!

# THE 'system()' COMMAND:
# #The system command is a way to use your computer terminal from within R. Within the system() function,
# whatever you put in quotes will be read as if that's what you typed in your terminal window. I'm doing all
# my analyses using a UNIX shell (bash), and so in order to use the PLINK program, I have to put a './'
# before I write 'plink'. Windows users only need to write 'plink'. In addition, I use 'intern=TRUE' so that
# the output will be displayed in the R console.
# #NOTE: this won't look as good in R as in terminal unless you have your width setup accordingly
options(width=100) #works well


#Let's see if these files will work in PLINK
# Once you've set up your binary directory for this analysis and have set R to that directory, do this:
#system("./plink --file test.snp --make-bed --out test2",intern=TRUE)

data <- read.table("practical_2/data/data.txt", header = T)
dim(data)
data[1:2, 1:24]

# Investigate the pi_hat_d values for the autosomes and X chromosome
data[1:2, 25:47]
# Take a look at the rest of the data matrix
data[1:2, 48:59]


# Row means over a matrix can be calculated with the below function
rowMeans(data[, 2:23])
# A similar function does not exist for row standard deviations and thus we
# must appeal to the apply function
apply(data[, 2:23], 1, sd)
pihat.mean <- rowMeans(data[, 2:23])
pihat.sd <- apply(data[, 2:23], 1, sd)
ibd.mean <- rowMeans(data[, 25:46])
ibd.sd <- apply(data[, 25:46], 1, sd)
plot(pihat.mean, ibd.mean, pch = 20,
col = 4, ylab = "Dominance relationship (mean IBD2 sharing)", xlab = "Additive relationship (mean IBD
          sharing)")
# Regress IBD2 on pi-hat. Use the lm function un R. Try ?lm if you are interested
reg <- lm(ibd.mean ~ pihat.mean)
abline(reg, lwd = 2.5)


# Regress IBD on pi-hat. Use the lm function un R. Try ?lm if you are interested
summary(lm(ibd.mean ~ pihat.mean))


# Read in .ped and .dat files
ped <- read.table("practical_2/data/qtdt.ped")
dat <- read.table("practical_2/data/qtdt.dat")
ibd <- read.table("practical_2/data/qtdt.ibd")

# The head of the .ped file. Two columns have been cut to fit in the listing
head(ped)
head(dat)
head(ibd)

# Let's see if these files will work in QTDT
# Once you've set up your binary directory for this analysis and have set R to that directory, do this:
system("./qtdt qtdt_mac -d practical_2/data/qtdt.dat -p practical_2/data/qtdt.ped -a- -we -veg")
system("./qtdt qtdt_mac -d practical_2/data/qtdt.dat -p practical_2/data/qtdt.ped -i practical_2/data/qtdt.ibd 3 -a- -weg -vega")


