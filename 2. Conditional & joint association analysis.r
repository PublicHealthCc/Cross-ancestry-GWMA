# Conditional 1

# EAS-LC
R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())
people <- 'EAS'
type <- 'LC'
sample_size <- 247989
gwas_tag <- '/data1/Trans_meta/meta_dataset_fixed/EAS_meta/EAS_meta_fixed.txt'
gwas <- fread(gwas_tag)
gwas$chr_bp <- paste0(gwas$CHR,':',gwas$BP)
data <- fread('trans_meta_sig_SNP.csv',fill = T)
data <- subset(data,Subtype==type)
data$chr_bp <- paste0(data$CHR,':',data$pos)
data1 <- subset(data,!no==people)
data2 <- subset(data,is.na(no))
dataa <- rbind(data1,data2)
for (rsID in dataa$rsid) {
  system(paste0('mkdir -p ',people,'/',type,'/',rsID))
  system(paste0('plink --noweb --bfile /data/Public/1000Genome/1kg_new/chr_all_1kgv3_2015_EAS --snp ',rsID,' --window 1000 --make-bed --out ',people,'/',type,'/',rsID,'/',rsID))
  dat <- fread(paste0(people,'/',type,'/',rsID,'/',rsID,'.bim'))
  names(dat) <- c('chr','SNP','CM','BP','ALT','REF')
  system(paste0('plink --noweb --bfile ',people,'/',type,'/',rsID,'/',rsID,' --freq --out ',people,'/',type,'/',rsID,'/',rsID))
  freq <- read.table(paste0(people,'/',type,'/',rsID,'/',rsID,'.frq'),h=T)
  dat1 <- merge(dat,freq,by='SNP')
  dat1$chr_bp <- paste0(dat1$chr,':',dat1$BP)
  gwas1 <- subset(gwas,chr_bp%in%dat1$chr_bp)
  names(dat1) <- paste0('ref_',names(dat1))
  names(dat1)[ncol(dat1)] <- 'chr_bp'
  dat2 <- merge(dat1, gwas1 ,by='chr_bp')
  a <- subset(dat2,ref_ALT==toupper(Allele1))
  b <- subset(dat2,ref_ALT!=toupper(Allele1))
  b$Effect <- -1*b$Effect
  dat3 <- rbind(a,b)
  re <- dat3[,c('ref_SNP','ref_A1','ref_A2','ref_MAF','Effect','StdErr','P')]
  colnames(re) <- c('SNP','A1','A2','freq','b','se','p')
  re$N <- sample_size
  write.table(re,file=paste0(people,'/',type,'/',rsID,'/',rsID,'.ma'),row.names=F,quote=F,sep='\t')
  write.table(rsID,file=paste0(people,'/',type,'/',rsID,'/cond.snplist.',rsID),row.names=F,quote=F,col.names=F,sep='\t')
  system(paste0('gcta  --bfile ',people,'/',type,'/',rsID,'/',rsID,' --cojo-wind 10000 --cojo-file ',people,'/',type,'/',rsID,'/',rsID,'.ma --cojo-cond ',people,'/',type,'/',rsID,'/cond.snplist.',rsID,' --out ',people,'/',type,'/',rsID,'/',rsID,'_condition.top'))
  con <- read.table(paste0(people,'/',type,'/',rsID,'/',rsID,'_condition.top.cma.cojo'),h=T)
  dat <- merge(re, con[,c('SNP','bC','bC_se','pC')], by='SNP',all.x=T)
  nrow(re)
  nrow(dat)
  head(dat)
  write.table(dat,file=paste0(people,'/',type,'/',rsID,'/',rsID,'_condition.top-sQTL.txt'),row.names=F,col.names=T,quote=F,sep='\t')
 
}

# 整理pC<5e-08的结果到EAS_LC_conditional_res.csv


# EAS-LUAD
R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())
people <- 'EAS'
type <- 'LUAD'
sample_size <- 29633
gwas_tag <- '/data1/Trans_meta/meta_dataset_fixed/EAS_meta/EAS_AD_meta_fixed.txt'
gwas <- fread(gwas_tag)
gwas$chr_bp <- paste0(gwas$CHR,':',gwas$BP)
data <- fread('trans_meta_sig_SNP.csv',fill = T)
data <- subset(data,Subtype==type)
data$chr_bp <- paste0(data$CHR,':',data$pos)
data1 <- subset(data,!no==people)
data2 <- subset(data,is.na(no))
dataa <- rbind(data1,data2)
for (rsID in dataa$rsid) {
  system(paste0('mkdir -p ',people,'/',type,'/',rsID))
  system(paste0('plink --noweb --bfile /data/Public/1000Genome/1kg_new/chr_all_1kgv3_2015_EAS --snp ',rsID,' --window 1000 --make-bed --out ',people,'/',type,'/',rsID,'/',rsID))
  dat <- fread(paste0(people,'/',type,'/',rsID,'/',rsID,'.bim'))
  names(dat) <- c('chr','SNP','CM','BP','ALT','REF')
  system(paste0('plink --noweb --bfile ',people,'/',type,'/',rsID,'/',rsID,' --freq --out ',people,'/',type,'/',rsID,'/',rsID))
  freq <- read.table(paste0(people,'/',type,'/',rsID,'/',rsID,'.frq'),h=T)
  dat1 <- merge(dat,freq,by='SNP')
  dat1$chr_bp <- paste0(dat1$chr,':',dat1$BP)
  gwas1 <- subset(gwas,chr_bp%in%dat1$chr_bp)
  names(dat1) <- paste0('ref_',names(dat1))
  names(dat1)[ncol(dat1)] <- 'chr_bp'
  dat2 <- merge(dat1, gwas1 ,by='chr_bp')
  a <- subset(dat2,ref_ALT==toupper(Allele1))
  b <- subset(dat2,ref_ALT!=toupper(Allele1))
  b$Effect <- -1*b$Effect
  dat3 <- rbind(a,b)
  re <- dat3[,c('ref_SNP','ref_A1','ref_A2','ref_MAF','Effect','StdErr','P')]
  colnames(re) <- c('SNP','A1','A2','freq','b','se','p')
  re$N <- sample_size
  write.table(re,file=paste0(people,'/',type,'/',rsID,'/',rsID,'.ma'),row.names=F,quote=F,sep='\t')
  write.table(rsID,file=paste0(people,'/',type,'/',rsID,'/cond.snplist.',rsID),row.names=F,quote=F,col.names=F,sep='\t')
  system(paste0('gcta  --bfile ',people,'/',type,'/',rsID,'/',rsID,' --cojo-wind 10000 --cojo-file ',people,'/',type,'/',rsID,'/',rsID,'.ma --cojo-cond ',people,'/',type,'/',rsID,'/cond.snplist.',rsID,' --out ',people,'/',type,'/',rsID,'/',rsID,'_condition.top'))
  con <- read.table(paste0(people,'/',type,'/',rsID,'/',rsID,'_condition.top.cma.cojo'),h=T)
  dat <- merge(re, con[,c('SNP','bC','bC_se','pC')], by='SNP',all.x=T)
  nrow(re)
  nrow(dat)
  head(dat)
  write.table(dat,file=paste0(people,'/',type,'/',rsID,'/',rsID,'_condition.top-sQTL.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}


# 整理pC<5e-08的结果到EAS_LUAD_conditional_res.csv


# EAS-LUSC
R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())
people <- 'EAS'
type <- 'LUSC'
sample_size <- 21800
gwas_tag <- '/data1/Trans_meta/meta_dataset_fixed/EAS_meta/EAS_SC_meta_fixed.txt'
gwas <- fread(gwas_tag)
gwas$chr_bp <- paste0(gwas$CHR,':',gwas$BP)
data <- fread('trans_meta_sig_SNP.csv',fill = T)
data <- subset(data,Subtype==type)
data$chr_bp <- paste0(data$CHR,':',data$pos)
data1 <- subset(data,!no==people)
data2 <- subset(data,is.na(no))
dataa <- rbind(data1,data2)
data <- fread('trans_meta_sig_SNP.csv',fill = T)
data <- subset(data,Subtype==type)
data$chr_bp <- paste0(data$CHR,':',data$pos)
data1 <- subset(data,!no==people)
data2 <- subset(data,is.na(no))
dataa <- rbind(data1,data2)

for (rsID in dataa$rsid) {
  system(paste0('mkdir -p ',people,'/',type,'/',rsID))
  system(paste0('plink --noweb --bfile /data/Public/1000Genome/1kg_new/chr_all_1kgv3_2015_EAS --snp ',rsID,' --window 1000 --make-bed --out ',people,'/',type,'/',rsID,'/',rsID))
  dat <- fread(paste0(people,'/',type,'/',rsID,'/',rsID,'.bim'))
  names(dat) <- c('chr','SNP','CM','BP','ALT','REF')
  system(paste0('plink --noweb --bfile ',people,'/',type,'/',rsID,'/',rsID,' --freq --out ',people,'/',type,'/',rsID,'/',rsID))
  freq <- read.table(paste0(people,'/',type,'/',rsID,'/',rsID,'.frq'),h=T)
  dat1 <- merge(dat,freq,by='SNP')
  dat1$chr_bp <- paste0(dat1$chr,':',dat1$BP)
  gwas1 <- subset(gwas,chr_bp%in%dat1$chr_bp)
  names(dat1) <- paste0('ref_',names(dat1))
  names(dat1)[ncol(dat1)] <- 'chr_bp'
  dat2 <- merge(dat1, gwas1 ,by='chr_bp')
  a <- subset(dat2,ref_ALT==toupper(Allele1))
  b <- subset(dat2,ref_ALT!=toupper(Allele1))
  b$Effect <- -1*b$Effect
  dat3 <- rbind(a,b)
  re <- dat3[,c('ref_SNP','ref_A1','ref_A2','ref_MAF','Effect','StdErr','P')]
  colnames(re) <- c('SNP','A1','A2','freq','b','se','p')
  re$N <- sample_size
  write.table(re,file=paste0(people,'/',type,'/',rsID,'/',rsID,'.ma'),row.names=F,quote=F,sep='\t')
  write.table(rsID,file=paste0(people,'/',type,'/',rsID,'/cond.snplist.',rsID),row.names=F,quote=F,col.names=F,sep='\t')
  system(paste0('gcta  --bfile ',people,'/',type,'/',rsID,'/',rsID,' --cojo-wind 10000 --cojo-file ',people,'/',type,'/',rsID,'/',rsID,'.ma --cojo-cond ',people,'/',type,'/',rsID,'/cond.snplist.',rsID,' --out ',people,'/',type,'/',rsID,'/',rsID,'_condition.top'))
  con <- read.table(paste0(people,'/',type,'/',rsID,'/',rsID,'_condition.top.cma.cojo'),h=T)
  dat <- merge(re, con[,c('SNP','bC','bC_se','pC')], by='SNP',all.x=T)
  nrow(re)
  nrow(dat)
  head(dat)
  write.table(dat,file=paste0(people,'/',type,'/',rsID,'/',rsID,'_condition.top-sQTL.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}


# 整理pC<5e-08的结果到EAS_LUSC_conditional_res.csv



# EUR-LC
R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())
people <- 'EUR'
type <- 'LC'
sample_size <- 291871
gwas_tag <- '/data1/Trans_meta/meta_dataset_fixed/EUR_meta/EUR_meta_fixed.txt'
gwas <- fread(gwas_tag)
gwas$chr_bp <- paste0(gwas$CHR,':',gwas$BP)
data <- fread('trans_meta_sig_SNP.csv',fill = T)
data <- subset(data,Subtype==type)
data$chr_bp <- paste0(data$CHR,':',data$pos)
data1 <- subset(data,!no==people)
data2 <- subset(data,is.na(no))
dataa <- rbind(data1,data2)

for (rsID in dataa$rsid) {
  system(paste0('mkdir -p ',people,'/',type,'/',rsID))
  system(paste0('plink --noweb --bfile /data/Public/1000Genome/1kg_new/chr_all_1kgv3_2015_EAS --snp ',rsID,' --window 1000 --make-bed --out ',people,'/',type,'/',rsID,'/',rsID))
  dat <- fread(paste0(people,'/',type,'/',rsID,'/',rsID,'.bim'))
  names(dat) <- c('chr','SNP','CM','BP','ALT','REF')
  system(paste0('plink --noweb --bfile ',people,'/',type,'/',rsID,'/',rsID,' --freq --out ',people,'/',type,'/',rsID,'/',rsID))
  freq <- read.table(paste0(people,'/',type,'/',rsID,'/',rsID,'.frq'),h=T)
  dat1 <- merge(dat,freq,by='SNP')
  dat1$chr_bp <- paste0(dat1$chr,':',dat1$BP)
  gwas1 <- subset(gwas,chr_bp%in%dat1$chr_bp)
  names(dat1) <- paste0('ref_',names(dat1))
  names(dat1)[ncol(dat1)] <- 'chr_bp'
  dat2 <- merge(dat1, gwas1 ,by='chr_bp')
  a <- subset(dat2,ref_ALT==toupper(Allele1))
  b <- subset(dat2,ref_ALT!=toupper(Allele1))
  b$Effect <- -1*b$Effect
  dat3 <- rbind(a,b)
  re <- dat3[,c('ref_SNP','ref_A1','ref_A2','ref_MAF','Effect','StdErr','P')]
  colnames(re) <- c('SNP','A1','A2','freq','b','se','p')
  re$N <- sample_size
  write.table(re,file=paste0(people,'/',type,'/',rsID,'/',rsID,'.ma'),row.names=F,quote=F,sep='\t')
  write.table(rsID,file=paste0(people,'/',type,'/',rsID,'/cond.snplist.',rsID),row.names=F,quote=F,col.names=F,sep='\t')
  system(paste0('gcta  --bfile ',people,'/',type,'/',rsID,'/',rsID,' --cojo-wind 10000 --cojo-file ',people,'/',type,'/',rsID,'/',rsID,'.ma --cojo-cond ',people,'/',type,'/',rsID,'/cond.snplist.',rsID,' --out ',people,'/',type,'/',rsID,'/',rsID,'_condition.top'))
  con <- read.table(paste0(people,'/',type,'/',rsID,'/',rsID,'_condition.top.cma.cojo'),h=T)
  dat <- merge(re, con[,c('SNP','bC','bC_se','pC')], by='SNP',all.x=T)
  nrow(re)
  nrow(dat)
  head(dat)
  write.table(dat,file=paste0(people,'/',type,'/',rsID,'/',rsID,'_condition.top-sQTL.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}


# 整理pC<5e-08的结果到EUR_LC_conditional_res.csv



# EUR-LUAD
R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())
people <- 'EUR'
type <- 'LUAD'
sample_size <- 266072
gwas_tag <- '/data1/Trans_meta/meta_dataset_fixed/EUR_meta/EUR_AD_meta_fixed.txt'
gwas <- fread(gwas_tag)
gwas$chr_bp <- paste0(gwas$CHR,':',gwas$BP)
data <- fread('trans_meta_sig_SNP.csv',fill = T)
data <- subset(data,Subtype==type)
data$chr_bp <- paste0(data$CHR,':',data$pos)
data1 <- subset(data,!no==people)
data2 <- subset(data,is.na(no))
dataa <- rbind(data1,data2)

for (rsID in dataa$rsid) {
  system(paste0('mkdir -p ',people,'/',type,'/',rsID))
  system(paste0('plink --noweb --bfile /data/Public/1000Genome/1kg_new/chr_all_1kgv3_2015_EUR --snp ',rsID,' --window 1000 --make-bed --out ',people,'/',type,'/',rsID,'/',rsID))
  dat <- fread(paste0(people,'/',type,'/',rsID,'/',rsID,'.bim'))
  names(dat) <- c('chr','SNP','CM','BP','ALT','REF')
  system(paste0('plink --noweb --bfile ',people,'/',type,'/',rsID,'/',rsID,' --freq --out ',people,'/',type,'/',rsID,'/',rsID))
  freq <- read.table(paste0(people,'/',type,'/',rsID,'/',rsID,'.frq'),h=T)
  dat1 <- merge(dat,freq,by='SNP')
  dat1$chr_bp <- paste0(dat1$chr,':',dat1$BP)
  gwas1 <- subset(gwas,chr_bp%in%dat1$chr_bp)
  names(dat1) <- paste0('ref_',names(dat1))
  names(dat1)[ncol(dat1)] <- 'chr_bp'
  dat2 <- merge(dat1, gwas1 ,by='chr_bp')
  a <- subset(dat2,ref_ALT==toupper(Allele1))
  b <- subset(dat2,ref_ALT!=toupper(Allele1))
  b$Effect <- -1*b$Effect
  dat3 <- rbind(a,b)
  re <- dat3[,c('ref_SNP','ref_A1','ref_A2','ref_MAF','Effect','StdErr','P')]
  colnames(re) <- c('SNP','A1','A2','freq','b','se','p')
  re$N <- sample_size
  write.table(re,file=paste0(people,'/',type,'/',rsID,'/',rsID,'.ma'),row.names=F,quote=F,sep='\t')
  write.table(rsID,file=paste0(people,'/',type,'/',rsID,'/cond.snplist.',rsID),row.names=F,quote=F,col.names=F,sep='\t')
  system(paste0('gcta  --bfile ',people,'/',type,'/',rsID,'/',rsID,' --cojo-wind 10000 --cojo-file ',people,'/',type,'/',rsID,'/',rsID,'.ma --cojo-cond ',people,'/',type,'/',rsID,'/cond.snplist.',rsID,' --out ',people,'/',type,'/',rsID,'/',rsID,'_condition.top'))
  con <- read.table(paste0(people,'/',type,'/',rsID,'/',rsID,'_condition.top.cma.cojo'),h=T)
  dat <- merge(re, con[,c('SNP','bC','bC_se','pC')], by='SNP',all.x=T)
  nrow(re)
  nrow(dat)
  head(dat)
  write.table(dat,file=paste0(people,'/',type,'/',rsID,'/',rsID,'_condition.top-sQTL.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}


# 整理pC<5e-08的结果到EUR_LUAD_conditional_res.csv


# EUR-LUSC
R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())
people <- 'EUR'
type <- 'LUSC'
sample_size <- 261951
gwas_tag <- '/data1/Trans_meta/meta_dataset_fixed/EUR_meta/EUR_SC_meta_fixed.txt'
gwas <- fread(gwas_tag)
gwas$chr_bp <- paste0(gwas$CHR,':',gwas$BP)
data <- fread('trans_meta_sig_SNP.csv',fill = T)
data <- subset(data,Subtype==type)
data$chr_bp <- paste0(data$CHR,':',data$pos)
data1 <- subset(data,!no==people)
data2 <- subset(data,is.na(no))
dataa <- rbind(data1,data2)

for (rsID in dataa$rsid) {
  system(paste0('mkdir -p ',people,'/',type,'/',rsID))
  system(paste0('plink --noweb --bfile /data/Public/1000Genome/1kg_new/chr_all_1kgv3_2015_EUR --snp ',rsID,' --window 1000 --make-bed --out ',people,'/',type,'/',rsID,'/',rsID))
  dat <- fread(paste0(people,'/',type,'/',rsID,'/',rsID,'.bim'))
  names(dat) <- c('chr','SNP','CM','BP','ALT','REF')
  system(paste0('plink --noweb --bfile ',people,'/',type,'/',rsID,'/',rsID,' --freq --out ',people,'/',type,'/',rsID,'/',rsID))
  freq <- read.table(paste0(people,'/',type,'/',rsID,'/',rsID,'.frq'),h=T)
  dat1 <- merge(dat,freq,by='SNP')
  dat1$chr_bp <- paste0(dat1$chr,':',dat1$BP)
  gwas1 <- subset(gwas,chr_bp%in%dat1$chr_bp)
  names(dat1) <- paste0('ref_',names(dat1))
  names(dat1)[ncol(dat1)] <- 'chr_bp'
  dat2 <- merge(dat1, gwas1 ,by='chr_bp')
  a <- subset(dat2,ref_ALT==toupper(Allele1))
  b <- subset(dat2,ref_ALT!=toupper(Allele1))
  b$Effect <- -1*b$Effect
  dat3 <- rbind(a,b)
  re <- dat3[,c('ref_SNP','ref_A1','ref_A2','ref_MAF','Effect','StdErr','P')]
  colnames(re) <- c('SNP','A1','A2','freq','b','se','p')
  re$N <- sample_size
  write.table(re,file=paste0(people,'/',type,'/',rsID,'/',rsID,'.ma'),row.names=F,quote=F,sep='\t')
  write.table(rsID,file=paste0(people,'/',type,'/',rsID,'/cond.snplist.',rsID),row.names=F,quote=F,col.names=F,sep='\t')
  system(paste0('gcta  --bfile ',people,'/',type,'/',rsID,'/',rsID,' --cojo-wind 10000 --cojo-file ',people,'/',type,'/',rsID,'/',rsID,'.ma --cojo-cond ',people,'/',type,'/',rsID,'/cond.snplist.',rsID,' --out ',people,'/',type,'/',rsID,'/',rsID,'_condition.top'))
  con <- read.table(paste0(people,'/',type,'/',rsID,'/',rsID,'_condition.top.cma.cojo'),h=T)
  dat <- merge(re, con[,c('SNP','bC','bC_se','pC')], by='SNP',all.x=T)
  nrow(re)
  nrow(dat)
  head(dat)
  write.table(dat,file=paste0(people,'/',type,'/',rsID,'/',rsID,'_condition.top-sQTL.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}


# 整理pC<5e-08的结果到EUR_LUSC_conditional_res.csv


###_______________________________________________________________________________________________________________________________________
# Conditional 2

# 根据*_conditional_res.csv的结果整理成trans_meta_sig_SNP_con1.csv


R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())

# LC
list <- read.csv('trans_meta_sig_SNP_con1.csv')
list1 <- subset(list,Subtype=='LC')
gwas_tag <- '/data1/Trans_meta/meta_dataset_fixed/Trans_meta/EAS_EUR_meta_fixed.txt'
gwas <- fread(gwas_tag)
gwas$chr_bp <- paste0(gwas$CHR,':',gwas$BP)
dat <- data.frame()
for (i in 1:nrow(list1)) {
  if(!is.na(list1$EAS_Conditional_1_sig[i]) & !is.na(list1$EUR_Conditional_1_sig[i])){
  dat1 <- fread(paste0('./con1/EAS/LC/',list1$rsid[i],'/',list1$rsid[i],'_condition.top.cma.cojo'))
  dat1 <- subset(dat1,pC<5e-08)
  dat2 <- fread(paste0('./con1/EUR/LC/',list1$rsid[i],'/',list1$rsid[i],'_condition.top.cma.cojo'))
  dat2 <- subset(dat2,pC<5e-08)
  dat3 <- rbind(dat1,dat2)
  dat3$chr_bp <- paste0(dat3$Chr,':',dat3$bp)
  write.table(dat3$SNP,file = paste0('./con1/LC/',list1$rsid[i],'_extract_list'),row.names = F,col.names = F,quote = F,sep = '\t')
  gwas1 <- subset(gwas,chr_bp%in%dat3$chr_bp)
  names(gwas1) <- paste0('meta_',names(gwas1))
  names(gwas1)[ncol(gwas1)] <- 'chr_bp'
  gwas1_sig <- subset(gwas1,meta_P<5e-08)
  dat4 <- merge(dat3,gwas1_sig,by = 'chr_bp')
  dat4 <- subset(dat4,!chr_bp%in%paste0(list1$CHR[i],':',list1$pos[i]))
  dat4 <- dat4[order(dat4$pC),]
  dat4 <- subset(dat4,!duplicated(SNP))
  write.csv(dat4,file = paste0('./con1/LC/',list1$rsid[i],'_con1_meta_sig.csv'),row.names = F)
  dat_out <- subset(dat4,meta_P==min(meta_P))
  dat_out$tag <- list1$rsid[i]
  dat <- rbind(dat,dat_out)
  } else if (is.na(list1$EAS_Conditional_1_sig[i])){
    dat2 <- fread(paste0('./con1/EUR/LC/',list1$rsid[i],'/',list1$rsid[i],'_condition.top.cma.cojo'))
    dat2 <- subset(dat2,pC<5e-08)
    dat3 <- dat2
    dat3$chr_bp <- paste0(dat3$Chr,':',dat3$bp)
    write.table(dat3$SNP,file = paste0('./con1/LC/',list1$rsid[i],'_extract_list'),row.names = F,col.names = F,quote = F,sep = '\t')
    gwas1 <- subset(gwas,chr_bp%in%dat3$chr_bp)
    names(gwas1) <- paste0('meta_',names(gwas1))
    names(gwas1)[ncol(gwas1)] <- 'chr_bp'
    gwas1_sig <- subset(gwas1,meta_P<5e-08)
    dat4 <- merge(dat3,gwas1_sig,by = 'chr_bp')
    dat4 <- subset(dat4,!chr_bp%in%paste0(list1$CHR[i],':',list1$pos[i]))
    dat4 <- dat4[order(dat4$pC),]
    dat4 <- subset(dat4,!duplicated(SNP))
    write.csv(dat4,file = paste0('./con1/LC/',list1$rsid[i],'_con1_meta_sig.csv'),row.names = F)
    dat_out <- subset(dat4,meta_P==min(meta_P))
    dat_out$tag <- list1$rsid[i]
    dat <- rbind(dat,dat_out)
  } else{
    dat1 <- fread(paste0('./con1/EAS/LC/',list1$rsid[i],'/',list1$rsid[i],'_condition.top.cma.cojo'))
    dat1 <- subset(dat1,pC<5e-08)
    dat3 <- dat1
    dat3$chr_bp <- paste0(dat3$Chr,':',dat3$bp)
    write.table(dat3$SNP,file = paste0('./con1/LC/',list1$rsid[i],'_extract_list'),row.names = F,col.names = F,quote = F,sep = '\t')
    gwas1 <- subset(gwas,chr_bp%in%dat3$chr_bp)
    names(gwas1) <- paste0('meta_',names(gwas1))
    names(gwas1)[ncol(gwas1)] <- 'chr_bp'
    gwas1_sig <- subset(gwas1,meta_P<5e-08)
    dat4 <- merge(dat3,gwas1_sig,by = 'chr_bp')
    dat4 <- subset(dat4,!chr_bp%in%paste0(list1$CHR[i],':',list1$pos[i]))
    dat4 <- dat4[order(dat4$pC),]
    dat4 <- subset(dat4,!duplicated(SNP))
    write.csv(dat4,file = paste0('./con1/LC/',list1$rsid[i],'_con1_meta_sig.csv'),row.names = F)
    dat_out <- subset(dat4,meta_P==min(meta_P))
    dat_out$tag <- list1$rsid[i]
    dat <- rbind(dat,dat_out)
  }
}


write.csv(dat,file = './con1/LC/con2_tag_SNP.csv',row.names = F)



# 整理成新的表格：conditional2_res_extract.csv



# EAS-LC
R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())
people <- 'EAS'
type <- 'LC'
sample_size <- 247989
gwas_tag <- '/data1/Trans_meta/meta_dataset_fixed/EAS_meta/EAS_meta_fixed.txt'
gwas <- fread(gwas_tag)
gwas$chr_bp <- paste0(gwas$CHR,':',gwas$BP)
data <- fread('conditional2_res_extract.csv',fill = T)
data <- subset(data,!Conditional_2_rsid==0)
dataa <- subset(data,Subtype==type)

for (i in 1:nrow(dataa)) {
  system(paste0('mkdir -p con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i]))
  dat <- fread(paste0('con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i],'.bim'))
  names(dat) <- c('chr','SNP','CM','BP','ALT','REF')
  freq <- read.table(paste0('con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i],'.frq'),h=T)
  dat1 <- merge(dat,freq,by='SNP')
  dat1$chr_bp <- paste0(dat1$chr,':',dat1$BP)
  gwas1 <- subset(gwas,chr_bp%in%dat1$chr_bp)
  names(dat1) <- paste0('ref_',names(dat1))
  names(dat1)[ncol(dat1)] <- 'chr_bp'
  dat2 <- merge(dat1, gwas1 ,by='chr_bp')
  a <- subset(dat2,ref_ALT==toupper(Allele1))
  b <- subset(dat2,ref_ALT!=toupper(Allele1))
  b$Effect <- -1*b$Effect
  dat3 <- rbind(a,b)
  dat3 <- subset(dat3,!ref_SNP==dataa$rsid[i])
  in_list <- fread(paste0('./con1/',type,'/',dataa$rsid[i],'_extract_list'),h = F)
  dat3 <- subset(dat3,ref_SNP%in%in_list$V1)
  re <- dat3[,c('ref_SNP','ref_A1','ref_A2','ref_MAF','Effect','StdErr','P')]
  colnames(re) <- c('SNP','A1','A2','freq','b','se','p')
  re$N <- sample_size
  write.table(re,file=paste0('con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i],'/',dataa$Conditional_2_rsid[i],'.ma'),row.names=F,quote=F,sep='\t')
  write.table(dataa$Conditional_2_rsid[i],file=paste0('con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i],'/cond.snplist.',dataa$Conditional_2_rsid[i]),row.names=F,quote=F,col.names=F,sep='\t')
  system(paste0('plink --bfile con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i], ' --extract ./con1/',type,'/',dataa$rsid[i],'_extract_list --make-bed --out con2/', people,'/',type,'/',dataa$Conditional_2_rsid[i],'/',dataa$Conditional_2_rsid[i]))
  system(paste0('gcta  --bfile con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i],'/',dataa$Conditional_2_rsid[i],' --cojo-wind 10000 --cojo-file con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i],'/',dataa$Conditional_2_rsid[i],'.ma --cojo-cond con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i],'/cond.snplist.',dataa$Conditional_2_rsid[i],' --out con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i],'/',dataa$Conditional_2_rsid[i],'_condition.top'))
  con <- read.table(paste0('con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i],'/',dataa$Conditional_2_rsid[i],'_condition.top.cma.cojo'),h=T)
  dat <- merge(re, con[,c('SNP','bC','bC_se','pC')], by='SNP',all.x=T)
  nrow(re)
  nrow(dat)
  head(dat)
  write.table(dat,file=paste0('con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i],'/',dataa$Conditional_2_rsid[i],'_condition.top-sQTL.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}

# 整理pC<5e-08的结果到EAS_LC_conditional_res2.csv



# EUR-LC
R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())
people <- 'EUR'
type <- 'LC'
sample_size <- 291871
gwas_tag <- '/data1/Trans_meta/meta_dataset_fixed/EUR_meta/EUR_meta_fixed.txt'
gwas <- fread(gwas_tag)
gwas$chr_bp <- paste0(gwas$CHR,':',gwas$BP)
data <- fread('conditional2_res_extract.csv',fill = T)
data <- subset(data,!Conditional_2_rsid==0)
dataa <- subset(data,Subtype==type)

for (i in 1:nrow(dataa)) {
  system(paste0('mkdir -p con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i]))
  dat <- fread(paste0('con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i],'.bim'))
  names(dat) <- c('chr','SNP','CM','BP','ALT','REF')
  freq <- read.table(paste0('con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i],'.frq'),h=T)
  dat1 <- merge(dat,freq,by='SNP')
  dat1$chr_bp <- paste0(dat1$chr,':',dat1$BP)
  gwas1 <- subset(gwas,chr_bp%in%dat1$chr_bp)
  names(dat1) <- paste0('ref_',names(dat1))
  names(dat1)[ncol(dat1)] <- 'chr_bp'
  dat2 <- merge(dat1, gwas1 ,by='chr_bp')
  a <- subset(dat2,ref_ALT==toupper(Allele1))
  b <- subset(dat2,ref_ALT!=toupper(Allele1))
  b$Effect <- -1*b$Effect
  dat3 <- rbind(a,b)
  dat3 <- subset(dat3,!ref_SNP==dataa$rsid[i])
  in_list <- fread(paste0('./con1/',type,'/',dataa$rsid[i],'_extract_list'),h = F)
  dat3 <- subset(dat3,ref_SNP%in%in_list$V1)
  re <- dat3[,c('ref_SNP','ref_A1','ref_A2','ref_MAF','Effect','StdErr','P')]
  colnames(re) <- c('SNP','A1','A2','freq','b','se','p')
  re$N <- sample_size
  write.table(re,file=paste0('con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i],'/',dataa$Conditional_2_rsid[i],'.ma'),row.names=F,quote=F,sep='\t')
  write.table(dataa$Conditional_2_rsid[i],file=paste0('con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i],'/cond.snplist.',dataa$Conditional_2_rsid[i]),row.names=F,quote=F,col.names=F,sep='\t')
  system(paste0('plink --bfile con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i], ' --extract ./con1/',type,'/',dataa$rsid[i],'_extract_list --make-bed --out con2/', people,'/',type,'/',dataa$Conditional_2_rsid[i],'/',dataa$Conditional_2_rsid[i]))
  system(paste0('gcta  --bfile con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i],'/',dataa$Conditional_2_rsid[i],' --cojo-wind 10000 --cojo-file con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i],'/',dataa$Conditional_2_rsid[i],'.ma --cojo-cond con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i],'/cond.snplist.',dataa$Conditional_2_rsid[i],' --out con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i],'/',dataa$Conditional_2_rsid[i],'_condition.top'))
  con <- read.table(paste0('con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i],'/',dataa$Conditional_2_rsid[i],'_condition.top.cma.cojo'),h=T)
  dat <- merge(re, con[,c('SNP','bC','bC_se','pC')], by='SNP',all.x=T)
  nrow(re)
  nrow(dat)
  head(dat)
  write.table(dat,file=paste0('con2/',people,'/',type,'/',dataa$Conditional_2_rsid[i],'/',dataa$Conditional_2_rsid[i],'_condition.top-sQTL.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}

# 整理pC<5e-08的结果到EUR_LC_conditional_res2.csv



###____________________________________________________________________________________________________________________________________________________________
# Conditional 3

# 根据*_conditional_res2.csv的结果整理成trans_meta_sig_SNP_con2.csv

R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())

# LC
list <- read.csv('trans_meta_sig_SNP_con2.csv')
list <- subset(list,!Conditional_2_rsid==0)
list1 <- list[-2,]
gwas_tag <- '/data1/Trans_meta/meta_dataset_fixed/Trans_meta/EAS_EUR_meta_fixed.txt'
gwas <- fread(gwas_tag)
gwas$chr_bp <- paste0(gwas$CHR,':',gwas$BP)

dat <- data.frame()
for (i in 1:nrow(list1)) {
  dat1 <- fread(paste0('/data1/qtl/test/eQTL_for_ZCC/trans_meta/con2/EAS/LC/',list1$Conditional_2_rsid[i],'/',list1$Conditional_2_rsid[i],'_condition.top.cma.cojo'))
  dat1 <- subset(dat1,pC<5e-08)
  dat2 <- fread(paste0('/data1/qtl/test/eQTL_for_ZCC/trans_meta/con2/EUR/LC/',list1$Conditional_2_rsid[i],'/',list1$Conditional_2_rsid[i],'_condition.top.cma.cojo'))
  dat2 <- subset(dat2,pC<5e-08)
  dat3 <- rbind(dat1,dat2)
  dat3$chr_bp <- paste0(dat3$Chr,':',dat3$bp)
  write.table(dat3$SNP,file = paste0('/data1/qtl/test/eQTL_for_ZCC/trans_meta/con2/LC/',list1$Conditional_2_rsid[i],'_extract_list'),row.names = F,col.names = F,quote = F,sep = '\t')
  gwas1 <- subset(gwas,chr_bp%in%dat3$chr_bp)
  names(gwas1) <- paste0('meta_',names(gwas1))
  names(gwas1)[ncol(gwas1)] <- 'chr_bp'
  gwas1_sig <- subset(gwas1,meta_P<5e-08)
  dat4 <- merge(dat3,gwas1_sig,by = 'chr_bp')
  dat4 <- subset(dat4,!chr_bp%in%list1$Conditional_2_chr_bp[i])
  dat4 <- dat4[order(dat4$pC),]
  dat4 <- subset(dat4,!duplicated(SNP))
  write.csv(dat4,file = paste0('/data1/qtl/test/eQTL_for_ZCC/trans_meta/con2/LC/',list1$Conditional_2_rsid[i],'_con2_meta_sig.csv'),row.names = F)
  dat_out <- subset(dat4,meta_P==min(meta_P))
  dat_out$tag <- list1$Conditional_2_rsid[i]
  dat <- rbind(dat,dat_out)
} 
write.csv(dat,file = '/data1/qtl/test/eQTL_for_ZCC/trans_meta/con2/LC/con3_tag_SNP1.csv',row.names = F)

# 整理成新的表格：conditional3_res_extract.csv

# EAS-LC
R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())
people <- 'EAS'
type <- 'LC'
sample_size <- 247989
gwas_tag <- '/data1/Trans_meta/meta_dataset_fixed/EAS_meta/EAS_meta_fixed.txt'
gwas <- fread(gwas_tag)
gwas$chr_bp <- paste0(gwas$CHR,':',gwas$BP)
data <- fread('conditional3_res_extract.csv',fill = T)
data <- subset(data,!Conditional_3_rsid==0)
dataa <- subset(data,Subtype==type)

for (i in 1:nrow(dataa)) {
  system(paste0('mkdir -p con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i]))
  dat <- fread(paste0('con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i],'.bim'))
  names(dat) <- c('chr','SNP','CM','BP','ALT','REF')
  freq <- read.table(paste0('con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i],'.frq'),h=T)
  dat1 <- merge(dat,freq,by='SNP')
  dat1$chr_bp <- paste0(dat1$chr,':',dat1$BP)
  gwas1 <- subset(gwas,chr_bp%in%dat1$chr_bp)
  names(dat1) <- paste0('ref_',names(dat1))
  names(dat1)[ncol(dat1)] <- 'chr_bp'
  dat2 <- merge(dat1, gwas1 ,by='chr_bp')
  a <- subset(dat2,ref_ALT==toupper(Allele1))
  b <- subset(dat2,ref_ALT!=toupper(Allele1))
  b$Effect <- -1*b$Effect
  dat3 <- rbind(a,b)
  dat3 <- subset(dat3,!ref_SNP==dataa$rsid[i])
  dat3 <- subset(dat3,!ref_SNP==dataa$Conditional_2_rsid[i])
  in_list <- fread(paste0('./con2/',type,'/',dataa$Conditional_2_rsid[i],'_extract_list'),h = F)
  dat3 <- subset(dat3,ref_SNP%in%in_list$V1)
  re <- dat3[,c('ref_SNP','ref_A1','ref_A2','ref_MAF','Effect','StdErr','P')]
  colnames(re) <- c('SNP','A1','A2','freq','b','se','p')
  re$N <- sample_size
  write.table(re,file=paste0('con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i],'/',dataa$Conditional_3_rsid[i],'.ma'),row.names=F,quote=F,sep='\t')
  write.table(dataa$Conditional_3_rsid[i],file=paste0('con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i],'/cond.snplist.',dataa$Conditional_3_rsid[i]),row.names=F,quote=F,col.names=F,sep='\t')
  system(paste0('plink --bfile con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i], ' --extract ./con2/',type,'/',dataa$Conditional_2_rsid[i],'_extract_list --make-bed --out con3/', people,'/',type,'/',dataa$Conditional_3_rsid[i],'/',dataa$Conditional_3_rsid[i]))
  system(paste0('gcta  --bfile con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i],'/',dataa$Conditional_3_rsid[i],' --cojo-wind 10000 --cojo-file con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i],'/',dataa$Conditional_3_rsid[i],'.ma --cojo-cond con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i],'/cond.snplist.',dataa$Conditional_3_rsid[i],' --out con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i],'/',dataa$Conditional_3_rsid[i],'_condition.top'))
  con <- read.table(paste0('con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i],'/',dataa$Conditional_3_rsid[i],'_condition.top.cma.cojo'),h=T)
  dat <- merge(re, con[,c('SNP','bC','bC_se','pC')], by='SNP',all.x=T)
  nrow(re)
  nrow(dat)
  head(dat)
  write.table(dat,file=paste0('con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i],'/',dataa$Conditional_3_rsid[i],'_condition.top-sQTL.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}

# 整理pC<5e-08的结果到EAS_LC_conditional_res3.csv


# EUR-LC
R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())
people <- 'EUR'
type <- 'LC'
sample_size <- 291871
gwas_tag <- '/data1/Trans_meta/meta_dataset_fixed/EUR_meta/EUR_meta_fixed.txt'
gwas <- fread(gwas_tag)
gwas$chr_bp <- paste0(gwas$CHR,':',gwas$BP)

data <- fread('conditional3_res_extract.csv',fill = T)
data <- subset(data,!Conditional_3_rsid==0)
dataa <- subset(data,Subtype==type)

for (i in 1:nrow(dataa)) {
  system(paste0('mkdir -p con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i]))
  dat <- fread(paste0('con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i],'.bim'))
  names(dat) <- c('chr','SNP','CM','BP','ALT','REF')
  freq <- read.table(paste0('con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i],'.frq'),h=T)
  dat1 <- merge(dat,freq,by='SNP')
  dat1$chr_bp <- paste0(dat1$chr,':',dat1$BP)
  gwas1 <- subset(gwas,chr_bp%in%dat1$chr_bp)
  names(dat1) <- paste0('ref_',names(dat1))
  names(dat1)[ncol(dat1)] <- 'chr_bp'
  dat2 <- merge(dat1, gwas1 ,by='chr_bp')
  a <- subset(dat2,ref_ALT==toupper(Allele1))
  b <- subset(dat2,ref_ALT!=toupper(Allele1))
  b$Effect <- -1*b$Effect
  dat3 <- rbind(a,b)
  dat3 <- subset(dat3,!ref_SNP==dataa$rsid[i])
  dat3 <- subset(dat3,!ref_SNP==dataa$Conditional_2_rsid[i])
  in_list <- fread(paste0('/./con2/',type,'/',dataa$Conditional_2_rsid[i],'_extract_list'),h = F)
  dat3 <- subset(dat3,ref_SNP%in%in_list$V1)
  re <- dat3[,c('ref_SNP','ref_A1','ref_A2','ref_MAF','Effect','StdErr','P')]
  colnames(re) <- c('SNP','A1','A2','freq','b','se','p')
  re$N <- sample_size
  write.table(re,file=paste0('con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i],'/',dataa$Conditional_3_rsid[i],'.ma'),row.names=F,quote=F,sep='\t')
  write.table(dataa$Conditional_3_rsid[i],file=paste0('con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i],'/cond.snplist.',dataa$Conditional_3_rsid[i]),row.names=F,quote=F,col.names=F,sep='\t')
  system(paste0('plink --bfile con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i], ' --extract ./con2/',type,'/',dataa$Conditional_2_rsid[i],'_extract_list --make-bed --out con3/', people,'/',type,'/',dataa$Conditional_3_rsid[i],'/',dataa$Conditional_3_rsid[i]))
  system(paste0('gcta  --bfile con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i],'/',dataa$Conditional_3_rsid[i],' --cojo-wind 10000 --cojo-file con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i],'/',dataa$Conditional_3_rsid[i],'.ma --cojo-cond con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i],'/cond.snplist.',dataa$Conditional_3_rsid[i],' --out con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i],'/',dataa$Conditional_3_rsid[i],'_condition.top'))
  con <- read.table(paste0('con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i],'/',dataa$Conditional_3_rsid[i],'_condition.top.cma.cojo'),h=T)
  dat <- merge(re, con[,c('SNP','bC','bC_se','pC')], by='SNP',all.x=T)
  nrow(re)
  nrow(dat)
  head(dat)
  write.table(dat,file=paste0('con3/',people,'/',type,'/',dataa$Conditional_3_rsid[i],'/',dataa$Conditional_3_rsid[i],'_condition.top-sQTL.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}

# 整理pC<5e-08的结果到EUR_LC_conditional_res3.csv



###________________________________________________________________________________________________________________________________
# Conditional 4

# 根据*_conditional_res5.csv的结果整理成trans_meta_sig_SNP_con4.csv
R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())

# LC
list <- read.csv('trans_meta_sig_SNP_con3.csv')
list <- subset(list,!Conditional_3_rsid==0)
list1 <- subset(list,Subtype=='LC')
list1 <- list1[1,]

gwas_tag <- '/data1/Trans_meta/meta_dataset_fixed/Trans_meta/EAS_EUR_meta_fixed.txt'
gwas <- fread(gwas_tag)
gwas$chr_bp <- paste0(gwas$CHR,':',gwas$BP)

dat <- data.frame()
for (i in 1:nrow(list1)) {
  dat1 <- fread(paste0('./con3/EAS/LC/',list1$Conditional_3_rsid[i],'/',list1$Conditional_3_rsid[i],'_condition.top.cma.cojo'))
  dat1 <- subset(dat1,pC<5e-08)
  dat2 <- fread(paste0('./con3/EUR/LC/',list1$Conditional_3_rsid[i],'/',list1$Conditional_3_rsid[i],'_condition.top.cma.cojo'))
  dat2 <- subset(dat2,pC<5e-08)
  dat3 <- rbind(dat1,dat2)
  dat3$chr_bp <- paste0(dat3$Chr,':',dat3$bp)
  write.table(dat3$SNP,file = paste0('./con3/LC/',list1$Conditional_3_rsid[i],'_extract_list'),row.names = F,col.names = F,quote = F,sep = '\t')
  gwas1 <- subset(gwas,chr_bp%in%dat3$chr_bp)
  names(gwas1) <- paste0('meta_',names(gwas1))
  names(gwas1)[ncol(gwas1)] <- 'chr_bp'
  gwas1_sig <- subset(gwas1,meta_P<5e-08)
  dat4 <- merge(dat3,gwas1_sig,by = 'chr_bp')
  dat4 <- subset(dat4,!chr_bp%in%list1$Conditional_3_chr_bp[i])
  dat4 <- dat4[order(dat4$pC),]
  dat4 <- subset(dat4,!duplicated(SNP))
  write.csv(dat4,file = paste0('./con3/LC/',list1$Conditional_3_rsid[i],'_con3_meta_sig.csv'),row.names = F)
  dat_out <- subset(dat4,meta_P==min(meta_P))
  dat_out$tag <- list1$Conditional_3_rsid[i]
  dat <- rbind(dat,dat_out)
}

write.csv(dat,file = './con3/LC/con4_tag_SNP.csv',row.names = F)

# 整理成新的表格：conditional4_res_extract.csv

# EAS-LC
R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())
people <- 'EAS'
type <- 'LC'
sample_size <- 247989
gwas_tag <- '/data1/Trans_meta/meta_dataset_fixed/EAS_meta/EAS_meta_fixed.txt'
gwas <- fread(gwas_tag)
gwas$chr_bp <- paste0(gwas$CHR,':',gwas$BP)
data <- fread('conditional4_res_extract.csv',fill = T)
data <- subset(data,!Conditional_4_rsid==0)
dataa <- subset(data,Subtype==type)

for (i in 1:nrow(dataa)) {
  system(paste0('mkdir -p con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i]))
  dat <- fread(paste0('con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i],'.bim'))
  names(dat) <- c('chr','SNP','CM','BP','ALT','REF')
  freq <- read.table(paste0('con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i],'.frq'),h=T)
  dat1 <- merge(dat,freq,by='SNP')
  dat1$chr_bp <- paste0(dat1$chr,':',dat1$BP)
  gwas1 <- subset(gwas,chr_bp%in%dat1$chr_bp)
  names(dat1) <- paste0('ref_',names(dat1))
  names(dat1)[ncol(dat1)] <- 'chr_bp'
  dat2 <- merge(dat1, gwas1 ,by='chr_bp')
  a <- subset(dat2,ref_ALT==toupper(Allele1))
  b <- subset(dat2,ref_ALT!=toupper(Allele1))
  b$Effect <- -1*b$Effect
  dat3 <- rbind(a,b)
  dat3 <- subset(dat3,!ref_SNP==dataa$rsid[i])
  dat3 <- subset(dat3,!ref_SNP==dataa$Conditional_2_rsid[i])
  dat3 <- subset(dat3,!ref_SNP==dataa$Conditional_3_rsid[i])
  in_list <- fread(paste0('./con3/',type,'/',dataa$Conditional_3_rsid[i],'_extract_list'),h = F)
  dat3 <- subset(dat3,ref_SNP%in%in_list$V1)
  re <- dat3[,c('ref_SNP','ref_A1','ref_A2','ref_MAF','Effect','StdErr','P')]
  colnames(re) <- c('SNP','A1','A2','freq','b','se','p')
  re$N <- sample_size
  write.table(re,file=paste0('con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i],'/',dataa$Conditional_4_rsid[i],'.ma'),row.names=F,quote=F,sep='\t')
  write.table(dataa$Conditional_4_rsid[i],file=paste0('con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i],'/cond.snplist.',dataa$Conditional_4_rsid[i]),row.names=F,quote=F,col.names=F,sep='\t')
  system(paste0('plink --bfile con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i], ' --extract ./con3/',type,'/',dataa$Conditional_3_rsid[i],'_extract_list --make-bed --out con4/', people,'/',type,'/',dataa$Conditional_4_rsid[i],'/',dataa$Conditional_4_rsid[i]))
  system(paste0('gcta  --bfile con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i],'/',dataa$Conditional_4_rsid[i],' --cojo-wind 10000 --cojo-file con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i],'/',dataa$Conditional_4_rsid[i],'.ma --cojo-cond con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i],'/cond.snplist.',dataa$Conditional_4_rsid[i],' --out con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i],'/',dataa$Conditional_4_rsid[i],'_condition.top'))
  con <- read.table(paste0('con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i],'/',dataa$Conditional_4_rsid[i],'_condition.top.cma.cojo'),h=T)
  dat <- merge(re, con[,c('SNP','bC','bC_se','pC')], by='SNP',all.x=T)
  nrow(re)
  nrow(dat)
  head(dat)
  write.table(dat,file=paste0('con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i],'/',dataa$Conditional_4_rsid[i],'_condition.top-sQTL.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}


# 整理pC<5e-08的结果到EAS_LC_conditional_res4.csv


# EUR-LC
R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())
people <- 'EUR'
type <- 'LC'
sample_size <- 291871
gwas_tag <- '/data1/Trans_meta/meta_dataset_fixed/EUR_meta/EUR_meta_fixed.txt'
gwas <- fread(gwas_tag)
gwas$chr_bp <- paste0(gwas$CHR,':',gwas$BP)
data <- fread('conditional4_res_extract.csv',fill = T)
data <- subset(data,!Conditional_4_rsid==0)
dataa <- subset(data,Subtype==type)

for (i in 1:nrow(dataa)) {
  system(paste0('mkdir -p con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i]))
  dat <- fread(paste0('con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i],'.bim'))
  names(dat) <- c('chr','SNP','CM','BP','ALT','REF')
  freq <- read.table(paste0('con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i],'.frq'),h=T)
  dat1 <- merge(dat,freq,by='SNP')
  dat1$chr_bp <- paste0(dat1$chr,':',dat1$BP)
  gwas1 <- subset(gwas,chr_bp%in%dat1$chr_bp)
  names(dat1) <- paste0('ref_',names(dat1))
  names(dat1)[ncol(dat1)] <- 'chr_bp'
  dat2 <- merge(dat1, gwas1 ,by='chr_bp')
  a <- subset(dat2,ref_ALT==toupper(Allele1))
  b <- subset(dat2,ref_ALT!=toupper(Allele1))
  b$Effect <- -1*b$Effect
  dat3 <- rbind(a,b)
  dat3 <- subset(dat3,!ref_SNP==dataa$rsid[i])
  dat3 <- subset(dat3,!ref_SNP==dataa$Conditional_2_rsid[i])
  dat3 <- subset(dat3,!ref_SNP==dataa$Conditional_3_rsid[i])
  in_list <- fread(paste0('./con3/',type,'/',dataa$Conditional_3_rsid[i],'_extract_list'),h = F)
  dat3 <- subset(dat3,ref_SNP%in%in_list$V1)
  re <- dat3[,c('ref_SNP','ref_A1','ref_A2','ref_MAF','Effect','StdErr','P')]
  colnames(re) <- c('SNP','A1','A2','freq','b','se','p')
  re$N <- sample_size
  write.table(re,file=paste0('con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i],'/',dataa$Conditional_4_rsid[i],'.ma'),row.names=F,quote=F,sep='\t')
  write.table(dataa$Conditional_4_rsid[i],file=paste0('con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i],'/cond.snplist.',dataa$Conditional_4_rsid[i]),row.names=F,quote=F,col.names=F,sep='\t')
  system(paste0('plink --bfile con1/',people,'/',type,'/',dataa$rsid[i],'/',dataa$rsid[i], ' --extract ./con3/',type,'/',dataa$Conditional_3_rsid[i],'_extract_list --make-bed --out con4/', people,'/',type,'/',dataa$Conditional_4_rsid[i],'/',dataa$Conditional_4_rsid[i]))
  system(paste0('gcta  --bfile con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i],'/',dataa$Conditional_4_rsid[i],' --cojo-wind 10000 --cojo-file con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i],'/',dataa$Conditional_4_rsid[i],'.ma --cojo-cond con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i],'/cond.snplist.',dataa$Conditional_4_rsid[i],' --out con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i],'/',dataa$Conditional_4_rsid[i],'_condition.top'))
  con <- read.table(paste0('con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i],'/',dataa$Conditional_4_rsid[i],'_condition.top.cma.cojo'),h=T)
  dat <- merge(re, con[,c('SNP','bC','bC_se','pC')], by='SNP',all.x=T)
  nrow(re)
  nrow(dat)
  head(dat)
  write.table(dat,file=paste0('con4/',people,'/',type,'/',dataa$Conditional_4_rsid[i],'/',dataa$Conditional_4_rsid[i],'_condition.top-sQTL.txt'),row.names=F,col.names=T,quote=F,sep='\t')

}


# 整理pC<5e-08的结果到EUR_LC_conditional_res4.csv


###______________________________________________________________________________________________________________________________________________________________________________________
# Conditional 5

# 根据*_conditional_res4.csv的结果整理成trans_meta_sig_SNP_con5.csv
R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())

# LC
list <- read.csv('trans_meta_sig_SNP_con4.csv')
list <- subset(list,!Conditional_4_rsid==0)
list1 <- subset(list,Subtype=='LC')
gwas_tag <- '/data1/Trans_meta/meta_dataset_fixed/Trans_meta/EAS_EUR_meta_fixed.txt'
gwas <- fread(gwas_tag)
gwas$chr_bp <- paste0(gwas$CHR,':',gwas$BP)

dat <- data.frame()
for (i in 1:nrow(list1)) {
  dat1 <- fread(paste0('./con4/EAS/LC/',list1$Conditional_4_rsid[i],'/',list1$Conditional_4_rsid[i],'_condition.top.cma.cojo'))
  dat1 <- subset(dat1,pC<5e-08)
  dat2 <- fread(paste0('./con4/EUR/LC/',list1$Conditional_4_rsid[i],'/',list1$Conditional_4_rsid[i],'_condition.top.cma.cojo'))
  dat2 <- subset(dat2,pC<5e-08)
  dat3 <- rbind(dat1,dat2)
  dat3$chr_bp <- paste0(dat3$Chr,':',dat3$bp)
  write.table(dat3$SNP,file = paste0('./con4/LC/',list1$Conditional_4_rsid[i],'_extract_list'),row.names = F,col.names = F,quote = F,sep = '\t')
  gwas1 <- subset(gwas,chr_bp%in%dat3$chr_bp)
  names(gwas1) <- paste0('meta_',names(gwas1))
  names(gwas1)[ncol(gwas1)] <- 'chr_bp'
  gwas1_sig <- subset(gwas1,meta_P<5e-08)
  dat4 <- merge(dat3,gwas1_sig,by = 'chr_bp')
  dat4 <- subset(dat4,!chr_bp%in%list1$Conditional_4_chr_bp[i])
  dat4 <- subset(dat4,!duplicated(SNP))
  write.csv(dat4,file = paste0('./con4/LC/',list1$Conditional_4_rsid[i],'_con4_meta_sig.csv'),row.names = F)
  dat_out <- subset(dat4,meta_P==min(meta_P))
  dat_out$tag <- list1$Conditional_4_rsid[i]
  dat <- rbind(dat,dat_out)
}

write.csv(dat,file = './con4/LC/con5_tag_SNP.csv',row.names = F)

