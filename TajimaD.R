########## r2vcftools installation --------------------------------------------------
## Install VCFtools
# GitHub: https://github.com/vcftools/vcftools
# Manual: https://vcftools.github.io/man_latest.html

## Install the r2vcftools library from GitHub
install.packages("devtools")
devtools::install_github("bcm-uga/LEA")
devtools::install_github("nspope/r2vcftools", force=T)

### Load r2vcftools
library(r2vcftools)

## Load VCFsummary function
VCFsummary <- function(snps){
  Nid <- nrow(snps@meta)
  Nsnps <- length(snps@site_id)
  cat(paste(Nid, "individuals and", Nsnps, "SNPs."), sep="\n")
}

########## Load VCF file -----------------------------------------------------------
## Download example VCF file "Imaurandioides.vcf" from figshare: https://ndownloader.figshare.com/files/10990757
url <- "https://ndownloader.figshare.com/files/10990757"
download.file(url, destfile	= "Imaurandioides.vcf")

## Load vcf file into vcfLink object              
snps <-  vcfLink("Imaurandioides.vcf", overwriteID=T)
snps
VCFsummary(snps) 

########## Filter to a single SNP per contig ---------------------------------------
## This thins SNPs to a given distance in bp from one another.
##  Setting the distance higher than the length of the contig ensures that you'll have a single SNP per contig.
snps_thin <- Filter(snps, filterOptions(thin=300)) 
VCFsummary(snps_thin) 


########## Subset by genetic clusters/populations ----------------------------------
## Create dummy population ID and insert into VCF meta data
snps_thin@meta$Pop <- rep(1:3, each = 10, length.out=nrow(snps_thin@meta))
head(snps_thin@meta)

##Subset by populations
P1 <- snps_thin@meta[snps_thin@meta$Pop==1, "sample_name"]
P2 <- snps_thin@meta[snps_thin@meta$Pop==2, "sample_name"]
P3 <- snps_thin@meta[snps_thin@meta$Pop==3, "sample_name"]

snps_p1 <- Subset(snps_thin, samples=P1)
snps_p2 <- Subset(snps_thin, samples=P2)
snps_p3 <- Subset(snps_thin, samples=P3)

VCFsummary(snps_p1) 
VCFsummary(snps_p2) 
VCFsummary(snps_p3) 

########## Estimate Tajima's D -----------------------------------------------------
## Estimate Tajima's D, bias-corrected for MAF
?TajimaD

tajd_p1 <- TajimaD(snps_p1, nboot=10000, maf=0.03, use_vcftools_D=FALSE)
tajd_p1$results
str(tajd_p1$simulations)

tajd_2 <- TajimaD(snps_p2, nboot=10000, maf=0.03, use_vcftools_D=FALSE)
tajd_p2$results
str(tajd_p2$simulations)

tajd_3 <- TajimaD(snps_p3, nboot=10000, maf=0.03, use_vcftools_D=FALSE)
tajd_p3$results
str(tajd_p3$simulations)

########## Saving and loading Tajima's D results -----------------------------------
## Save
save(tajd_p1, file = "tajd_p1.Rdata")
save(tajd_p2, file = "tajd_p2.Rdata")
save(tajd_p3, file = "tajd_p3.Rdata")

## Load 
load("tajd_p1.Rdata")
load("tajd_p2.Rdata")
load("tajd_p3.Rdata")

########## Plotting observed Tajima's D against the null distribution -------------------
## The null distribution (histogram) is shown next to the observed Tajima's D value (red line)
library(ggplot2)

ggplot(data.frame(x=tajd_p1$simulations$'Null, bias-corrected')) + geom_histogram(aes(x=x), binwidth=0.01) + geom_vline(xintercept=mean(tajd_p1$simulations$'Bootstrap, bias-corrected'), lty=2, col="red") + geom_vline(xintercept=0) + theme_bw()
ggplot(data.frame(x=tajd_p2$simulations$'Null, bias-corrected')) + geom_histogram(aes(x=x), binwidth=0.01) + geom_vline(xintercept=mean(tajd_p2$simulations$'Bootstrap, bias-corrected'), lty=2, col="red") + geom_vline(xintercept=0) + theme_bw()
ggplot(data.frame(x=tajd_p3$simulations$'Null, bias-corrected')) + geom_histogram(aes(x=x), binwidth=0.01) + geom_vline(xintercept=mean(tajd_p3$simulations$'Bootstrap, bias-corrected'), lty=2, col="red") + geom_vline(xintercept=0) + theme_bw()

