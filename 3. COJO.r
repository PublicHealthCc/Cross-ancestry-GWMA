# GCTA-COJO analysis
# Conditional 1
# Give an example: EAS-LC
R
library(data.table)
options(stringsAsFactors=F)
setwd('./')	
rm(list=ls())
people <- 'EAS'
type <- 'LC'
sample_size <- 247989
gwas_tag <- './EAS_meta_fixed.txt'
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
  system(paste0('plink --noweb --bfile ./chr_all_1kgv3_2015_EAS --snp ',rsID,' --window 1000 --make-bed --out ',people,'/',type,'/',rsID,'/',rsID))
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
  write.table(dat,file=paste0(people,'/',type,'/',rsID,'/',rsID,'_condition.top.txt'),row.names=F,col.names=T,quote=F,sep='\t')
 
}
