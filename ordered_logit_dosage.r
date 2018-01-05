args <- commandArgs(trailingOnly = TRUE)

genotype_file <- args[1]
fam_file <- args[2]
pheno_file <- args[3]
covar_file <- args[4]
event_name <- args[5]
time_name <- args[6]
cov_list <- args[7]
outfile <- args[8]

# event_name="ftremit_event"
# time_name="ftremit_time"
# cov_list=c("C1","C2","C3")

library(data.table)
library(survival)


genotype <- fread(paste('gunzip -c',genotype_file) ,header=F,data.table=F)
fam <- read.table(fam_file,stringsAsFactors=F,header=F)
pheno <- read.table(pheno_file,stringsAsFactors=F,header=T)
covar <- read.table(covar_file,stringsAsFactors=F,header=T)

# genotype <- fread('gunzip -c /home/caroline/Desktop/pgbd/qc/imputation/dasuqc1_bip_pgbd_mix_am-qc.hg19.ch.fl/dosecnt/dos_bip_pgbd_mix_am-qc.hg19.ch.fl.chr10_027_030.out.dosage.gz.doscnt.gz',header=T,data.table=F)
# fam <- read.table('/home/caroline/Desktop/pgbd/qc/imputation/dasuqc1_bip_pgbd_mix_am-qc.hg19.ch.fl/qc1/dos_bip_pgbd_mix_am-qc.hg19.ch.fl.chr5_123_126.out.dosage.fam',stringsAsFactors=F,header=F)
# pheno <- read.table('/home/caroline/Desktop/pgbd/qc/imputation/ftremit.pheno.txt',stringsAsFactors=F,header=T)
# covar <- read.table('/home/caroline/Desktop/pgbd/qc/imputation/bip_pgbd_mix_am-qc-eur_pca.menv.mds_cov',stringsAsFactors=F,header=T)
# time_name="ftremit_time"
# event_name="ftremit_event"


covar$FID_IID <- paste(covar[,1],covar[,2],sep="_")
pheno$FID_IID <- paste(pheno[,1],pheno[,2],sep="_")
phenocov <- merge(pheno,covar,by="FID_IID")
phenocov$ev <- Surv(time=phenocov[,time_name],event=phenocov[,event_name])
#phenocov <- phenocov[!is.na(phenocov[,event_name]),] #Only take subjects with events

#Subset data frame to relevant subjects only

names(genotype)[-c(1:3)] <- paste(fam[,1],fam[,2],sep="_")

#Get marker names 
markers <- genotype[,1]

#Transpose genotype
genotype2 <- data.frame(t(genotype[-c(1:3)]))
names(genotype2) <- c(markers)


#make a dataframe of subject names that goes in the order they appear in the file
dat_order <- as.data.frame(rownames(genotype2))

#name the subjects the same as found in the matching column in the batch file
names(dat_order) <- c('FID_IID')

#establish a note of what the ordering is
dat_order$order <- c(1:(dim(genotype2)[1]))


#merge batch data with subject name order data
covariates_ordered0 <- merge(phenocov,dat_order,by='FID_IID')

#sort the data by order so the phenotype now lines up
covariates_ordered <- covariates_ordered0[order(covariates_ordered0$order),]


#Filter down genotypes to ones with covariates
genotype_3a <- genotype2[row.names(genotype2) %in% covariates_ordered$FID_IID,]
row.names(genotype_3a) == covariates_ordered$FID_IID #Should be true

#Remove monomorphic variants
 monocheck <- function(x,tol)
 {
  abs(max(x,na.rm=T) - min(x,na.rm=T)) <= tol
 }

monmarkers <- apply(genotype_3a,2,monocheck,tol=0.000000001)
genotype_3 <- genotype_3a[,!monmarkers]


#Define a matrix of covariates
covs_use <- strsplit(cov_list,",")[[1]]
print(paste("Using covariates", cov_list))
covmat <- as.matrix(covariates_ordered[,covs_use])
#covmat <-  as.matrix(covariates_ordered[,c("C1","C2","C3")])



#Define a null model for the lrtest
nullmod <- coxph(phenocov$ev ~ covmat )


regfun <- function(genotype,event,covset,nullmod)
{    

    try(
        gee.fit <- coxph(event ~ covset + genotype)
        , silent = TRUE)
                
                      
    if(exists("gee.fit"))
    {
        if(class(gee.fit) !=  "try-error")
        {
         lrt <- anova(gee.fit)[3,4] # Likelihood ratio test p-value
         r <- c(coef(summary(gee.fit))["genotype",],lrt)
      
        } 
    } else {
                r <- rep(NA,6) #If the model didn't fit, just supply impossible values
        }
    return(r)
}


results <- t(apply(genotype_3,2,regfun,event=phenocov$ev,covset=covmat,nullmod=nullmod))
#results <- t(apply(genotype_3[,1:100],2,regfun,event=phenocov$ev,covset=covmat,nullmod=nullmod))
#coxph(

write.table(results,outfile)

        
