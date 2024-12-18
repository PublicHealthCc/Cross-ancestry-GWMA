###______________________________________________________________________________________________________________________
####1. Enrichment analysis of gene pathways#############################################################################
R
library(data.table)
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

####Step1. level1: susie missense mutation genes + colocalization genes with lung and brain tissues
data<-fread('~/level1.csv')
entrezid_all = mapIds(x = org.Hs.eg.db,  
                      keys = data$symbol, 
                      keytype = "SYMBOL", 
                      column = "ENTREZID")
entrezid_all = data.frame(entrezid_all) 
entrezid_all['AC008537.5','entrezid_all']<-6439
entrezid_all['AL133458.1','entrezid_all']<-18570
entrezid_all  = na.omit(entrezid_all) 
head(entrezid_all)
###GO
GO_enrich = enrichGO(gene = entrezid_all[,1], 
                     OrgDb = org.Hs.eg.db, 
                     keyType = "ENTREZID", 
                     ont = "ALL", 
                     pvalueCutoff = 1,qvalueCutoff = 1, 
                     readable = T) 
GO_enrich  = data.frame(GO_enrich) 
fwrite(GO_enrich,'GO_level1_passway_result.csv',sep=',',col.names=T)

####Step2. level2
data<-fread('~/level2.csv')
entrezid_all = mapIds(x = org.Hs.eg.db,  
                      keys = data$symbol, 
                      keytype = "SYMBOL", 
                      column = "ENTREZID")
entrezid_all['AC008537.5']<-'6439'
entrezid_all['AL133458.1']<-'18570'
entrezid_all['HIST1H2BK']<-'85236'
entrezid_all['HIST1H2BL']<-'8340'
entrezid_all['FGFR1OP']<-'421568'
entrezid_all['OBFC1']<-'79991'
entrezid_all['AMICA1']<-'120425'
entrezid_all['AC027228.2']<-'17708'
entrezid_all  = na.omit(entrezid_all)  
entrezid_all = data.frame(entrezid_all) 
head(entrezid_all)
########GO
GO_enrich = enrichGO(gene = entrezid_all[,1], 
                     OrgDb = org.Hs.eg.db, 
                     keyType = "ENTREZID", 
                     ont = "ALL", 
                     pvalueCutoff = 1,qvalueCutoff = 1, 
                     readable = T) 
GO_enrich  = data.frame(GO_enrich) 
fwrite(GO_enrich,'GO_level2_passway_addBPTF_result.csv',sep=',',col.names=T)


###______________________________________________________________________________________________________________________
####2. Prepare data for gene-set-based PRS#############################################################################
###Step1. extract TCGA genotype data#######################################################
for i in {1..22}; do
plink2 --~/gwas_data_chr$i \
--extract extract_TCGA_list_35snps \
--export A \
--out snps33_for_prs_chr$i
done

# merge data
R
library(data.table)
rm(list = ls())
list <- fread('raw_list',h = F)
res <- fread('snps33_for_prs_chr10.raw')

for (i in 2:nrow(list)) {
  res1 <- fread(list$V1[i])
  res1 <- res1[,-c(1:6)]
  res <- cbind(res,res1)
}
write.table(res,file = 'TCGA_35SNPs_for_prs_dosage.txt',sep = '\t',quote = F,row.names = F,col.names = T)

# covariate data prepare###########################################################
rm(list=ls())
library(data.table)
library(dplyr)

cov_lusc <- fread('TCGA-LUSC.GDC_phenotype.tsv',check.names = F) # UCSC
cov_luad <- fread('TCGA-LUAD.GDC_phenotype.tsv',check.names = F) # UCSC
cov_nsclc <- rbind(cov_luad,cov_lusc)
table(cov_nsclc$tumor_stage.diagnoses)

cov_nsclc$smoking_status1 <- cov_nsclc$tobacco_smoking_history
cov_nsclc$py <- cov_nsclc$pack_years_smoked.exposures
a <- subset(cov_nsclc,smoking_status1==1)
cov_nsclc$py[cov_nsclc$smoking_status1==1 & is.na(cov_nsclc$pack_years_smoked.exposures)] <- 0
a <- subset(cov_nsclc,smoking_status1==2)
cov_nsclc$py[cov_nsclc$smoking_status1==2 & is.na(cov_nsclc$pack_years_smoked.exposures)] <- median(na.omit(a$pack_years_smoked.exposures))
a <- subset(cov_nsclc,smoking_status1==3)
cov_nsclc$py[cov_nsclc$smoking_status1==3 & is.na(cov_nsclc$pack_years_smoked.exposures)] <- median(na.omit(a$pack_years_smoked.exposures))
a <- subset(cov_nsclc,smoking_status1==4)
cov_nsclc$py[cov_nsclc$smoking_status1==4 & is.na(cov_nsclc$pack_years_smoked.exposures)] <- median(na.omit(a$pack_years_smoked.exposures))
a <- subset(cov_nsclc,smoking_status1==5)
cov_nsclc$py[cov_nsclc$smoking_status1==5 & is.na(cov_nsclc$pack_years_smoked.exposures)] <- median(na.omit(a$pack_years_smoked.exposures))


cov_nsclc$smoking_status2 <- ifelse(cov_nsclc$smoking_status1==1,1,
                                    ifelse(cov_nsclc$smoking_status1%in%c(2,3,4,5) & cov_nsclc$py<30,2,
                                           ifelse(cov_nsclc$smoking_status1%in%c(2,3,4,5) & cov_nsclc$py>=30,3,NA)))
cov_nsclc$smoking_status2[is.na(cov_nsclc$smoking_status2)]=999
names(cov_nsclc)[12] <- 'FID'

cov_cell <- fread('cell_TCGA_cov.txt') # Cell 2016
names(cov_cell)[1] <- 'FID'
cov_cell <- subset(cov_cell,type%in%c('LUAD','LUSC'))
cov_cell$type <- ifelse(cov_cell$type=='LUAD',1,2)
cov_nsclc <- merge(cov_nsclc,cov_cell,by = 'FID')

cov_nsclc$ajcc_stage <- ifelse(cov_nsclc$ajcc_pathologic_tumor_stage%in% c('Stage I','Stage IA','Stage IB'),1,
                               ifelse(cov_nsclc$ajcc_pathologic_tumor_stage%in% c('Stage II','Stage IIA','Stage IIB'),2,
                                      ifelse(cov_nsclc$ajcc_pathologic_tumor_stage%in% c('Stage III','Stage IIIA','Stage IIIB'),3,
                                             ifelse(cov_nsclc$ajcc_pathologic_tumor_stage%in% c('Stage IV'),4,999))))

cov_nsclc$gender_use <- ifelse(cov_nsclc$gender=='MALE',1,2)
cov_nsclc$age_use <- -cov_nsclc$birth_days_to/365.25

save(cov_nsclc,file = 'TCGA_NSCLC_cov.rda')

# TMB calculation ###########################################################
library(maftools)
library(dplyr)
rm(list = ls())

mafs = list.files(path = "./", pattern = "*.\\.gz$", full.names = TRUE) 
mafs

all_mut <- merge_mafs(mafs = mafs)

a <- all_mut@data %>%
  .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")] %>%
  as.data.frame() %>%
  mutate(Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode,1,12))

gene <- as.character(unique(a$Hugo_Symbol))
sample <- as.character(unique(a$Tumor_Sample_Barcode))

mat <- as.data.frame(matrix("",length(gene),length(sample),
                            dimnames = list(gene,sample)))
mat_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                dimnames = list(gene,sample)))

for (i in 1:nrow(a)){
  mat[as.character(a[i,1]),as.character(a[i,3])] <- as.character(a[i,2])
}

for (i in 1:nrow(a)){
  mat_0_1[as.character(a[i,1]),as.character(a[i,3])] <- 1
}

write.csv(mat,"NSCLC_all_mut_type.csv")
write.csv(mat_0_1,"NSCLC_all_mut_01.csv")

tmb_table = tmb(maf = all_mut,logScale = F)   # no log change
write.csv(tmb_table,"NSCLC_tmb_results.csv")


# SBS denovo calculation ###########################################################
rm(list=ls())
library(data.table)
library(deconstructSigs)
library(maftools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)
mafs = list.files(path = "./", pattern = "*.\\.gz$", full.names = TRUE)
mafs
all_mut <- merge_mafs(mafs = mafs)
sample.mut.ref <- all_mut@data
table(sample.mut.ref$Variant_Type)
sample.mut.ref <- subset(sample.mut.ref,Variant_Type=='SNP')

sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref,
                                sample.id = "Tumor_Sample_Barcode",
                                chr = "Chromosome",
                                pos = "Start_Position",
                                ref = "Reference_Allele",
                                alt = "Tumor_Seq_Allele2",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)

class(sigs.input)
save(sigs.input,file = 'mut.to.sigs.input.rda')


w=lapply(unique(sample.mut.ref$Tumor_Sample_Barcode) , function(i){
  sample_1 = whichSignatures(tumor.ref = sigs.input[,], 
                             signatures.ref = signatures.cosmic, 
                             sample.id =  i, 
                             contexts.needed = TRUE,
                             tri.counts.method = 'default')
  print(i)
  return(sample_1$weights)
})
w=do.call(rbind,w)
library(pheatmap)
pheatmap(t(w),cluster_rows = F,cluster_cols = T)
pheatmap(w,cluster_rows = T,cluster_cols = F)
mut.wt=w
save(mut.wt,file = 'wgs-mut.wt.Rdata')

###Step2. Calculate gene-set-based PRS
data=read.table("TCGA_35SNPs_for_prs_dosage.txt",h=T,check.names = F)

data$rs71658797 <- 2 - data$`rs71658797:77967507:T:A_T`
data$rs72822431 <- 2 - data$`rs72822431:65521816:T:G_T`
data$rs7621631 <- data$`rs7621631:169512145:C:A_C`
data$rs7631358 <- 2 - data$`rs7631358:189348411:G:A_G`
data$rs2853677 <- data$`rs2853677:1287194:G:A_G`
data$rs3131856 <- 2 - data$`rs3131856:29607101:T:C_T`
data$rs1977358 <- data$`rs1977358:41489854:T:C_T`
data$rs2883418 <- data$`rs2883418:117791751:T:C_T`
data$rs6920364 <- 2 - data$`rs6920364:167376466:G:C_G`
data$rs7820838 <- data$`rs7820838:32405979:T:C_T`
data$rs10811609 <- 2 - data$`rs10811609:21746358:A:G_A`
data$rs4879704 <- data$rs4879704_A
data$rs7902587 <- 2 - data$rs7902587_C
data$rs11196088 <- 2 - data$`rs11196088:114507733:G:A_G`
data$rs1629083 <- data$`rs1629083:118126576:C:T_C`
data$rs4938515 <- 2 - data$`rs4938515:118478404:C:A_C`
data$rs3748522 <- 2 - data$`rs3748522:1058688:A:C_A`
data$rs11612312 <- data$`rs11612312:52349088:T:C_T` 
data$rs11571818 <- 2 - data$`rs11571818:32968810:T:C_T`
data$rs55980174 <- 2 - data$`rs55980174:35216499:T:A_T` 
data$rs77468143 <- data$`rs77468143:49376624:T:G_T`
data$rs8029744 <- 2 - data$`rs8029744:75029407:A:C_A`
data$rs55781567 <- 2 - data$`rs55781567:78857986:C:G_C`
data$rs2967352 <- data$`rs2967352:82196676:T:C_T`
data$rs1976052 <- 2 - data$`rs1976052:65830578:C:T_C`
data$rs56113850 <- 2 - data$`rs56113850:41353107:T:C_T`
data$rs1291141 <- 2 - data$`rs1291141:35518439:T:G_T`
data$rs6011779 <- data$`rs6011779:61984317:C:T_C`
data$rs17879961 <- data$`rs17879961:29121087:A:G_A`
data$rs7705526 <- 2 - data$`rs7705526:1285974:C:A_C`
data$rs401681 <- data$rs401681_C
data$rs13167280 <- 2 - data$`rs13167280:1280477:G:A_G`
data$rs2736109 <- 2 - data$`rs2736109:1296759:C:T_C`
data$rs938682 <- 2 - data$rs938682_G
data$rs76712448 <- 2 - data$`rs76712448:78836719:A:G_A`

data$PRS1 <- 0.1179*data$rs71658797 + 0.0687*data$rs7621631+ 0.1743*data$rs2853677 +
  0.1733*data$rs3131856 + 0.1695*data$rs7902587 + 0.0583*data$rs4938515 +
  0.0538*data$rs3748522 + 0.0892*data$rs2967352 + 0.0587*data$rs1291141 +
  0.17*data$rs7705526 + 0.1823*data$rs13167280 + 0.1325*data$rs2736109

data$PRS2 <-  0.0754*data$rs72822431 + 0.1344*data$rs7631358 + 0.1126*data$rs1977358 +
  0.1023*data$rs2883418 + 0.0771*data$rs7820838 + 0.0913*data$rs4879704 +
  0.1518*data$rs11196088 + 0.0696*data$rs1629083 + 0.0797*data$rs11612312 + 0.0782*data$rs1976052

data$PRS3 <-  0.0577*data$rs6920364 + 0.0845*data$rs10811609 + 0.4627*data$rs11571818 +
  0.0706*data$rs55980174 + 0.0804*data$rs77468143 + 0.0655*data$rs8029744 +
  0.3969*data$rs17879961 + 0.1344*data$rs401681 

data$PRS4 <-  0.2535*data$rs55781567 + 0.0903*data$rs56113850 + 0.078*data$rs6011779 +
  0.1319*data$rs938682 + 0.0633*data$rs76712448

save(data,file = 'PRS_4_flip.rda')

# Combine data################################################################
rm(list=ls())
library(data.table)
library(maftools)
library(dplyr)
options(scipen = 200)

load("PRS_4_flip.rda")
load('TCGA_NSCLC_cov.rda')
cov_nsclc <- subset(cov_nsclc,!duplicated(FID))
dat <- dplyr::left_join(data,cov_nsclc,by = 'FID')
dat$age_use[is.na(dat$age_use)] <- mean(dat$age_use,na.rm = T)

# single gene mutation
list <- read.table('Gene.txt',h=F)
maf <- read.csv('NSCLC_all_mut_01.csv',row.names = 1,check.names = F)
maf <- t(maf)
maf <- cbind(FID=rownames(maf),maf)
maf <- as.data.frame(maf)
names(maf) <- gsub('-','_',names(maf),fixed = T)
maf <- maf[,c('FID',list$V1)]
maf[,2:ncol(maf)] <- apply(maf[,2:ncol(maf)],2,as.integer)
dat=merge(dat,maf,by = 'FID')

# TMB
tmb <- read.csv('NSCLC_tmb_results.csv',row.names = 1)
tmb <- tmb[order(tmb$total_perMB),]
tmb <- subset(tmb,sapply(strsplit(tmb$Tumor_Sample_Barcode,'-'),'[',4)=='01A')
tmb$FID <- paste0(sapply(strsplit(tmb$Tumor_Sample_Barcode,'-'),'[',1),'-',sapply(strsplit(tmb$Tumor_Sample_Barcode,'-'),'[',2),'-',sapply(strsplit(tmb$Tumor_Sample_Barcode,'-'),'[',3))
dat <- dplyr::left_join(dat,tmb,by = 'FID')
dat <- subset(dat,!is.na(dat$total_perMB))

# SBS
load('wgs-mut.wt.Rdata')
mut.wt$FID <- paste0(sapply(strsplit(rownames(mut.wt),'-'),'[',1),'-',sapply(strsplit(rownames(mut.wt),'-'),'[',2),'-',sapply(strsplit(rownames(mut.wt),'-'),'[',3))
dat <- dplyr::left_join(dat,mut.wt,by = 'FID')
dat <- subset(dat,!is.na(dat$Signature.30))

dat$race_use <- ifelse(dat$race%in%c('WHITE'),1,
                       ifelse(dat$race%in%c('AMERICAN INDIAN OR ALASKA NATIVE','BLACK OR AFRICAN AMERICAN','ASIAN'),2,999))

save(dat,file = 'PRS_4_flip_tmb_sbs_age_race_impute.rda')


####Step3. PRS association
####(1) singe gene mutation
rm(list=ls())
library(data.table)
library(dplyr)

load('PRS_4_flip_tmb_sbs_age_race_impute.rda')

prs_list <- c('PRS1','PRS2','PRS3','PRS4')
for (j in prs_list) {
  
  res_out1 <- data.frame()
  
  for (i in c(242:261)){
    dat1 <- dat[,c(names(dat)[i],j,'age_use','gender_use','smoking_status2','type','race_use')]
    dat1[,1] <- as.integer(dat1[,1])
    dat1$gender_use <- as.factor(dat1$gender_use)
    dat1$smoking_status2 <- as.factor(dat1$smoking_status2)
    dat1$type <- as.factor(dat1$type)
    dat1$race_use <- as.factor(dat1$race_use)
    formula <- as.formula(paste0(names(dat)[i], "~ ",j," + age_use + gender_use + smoking_status2 +type +race_use"))
    fit1 <- glm(formula,family=c("binomial"),data=dat1)
    res <- summary(fit1)
    res <- as.data.frame(res$coefficients)
    res <- cbind(group=rownames(res),res)
    res$gene <- names(dat)[i]
    res_out1 <- rbind(res_out1,res)
    print(i)
  }
  write.csv(res_out1,file =  paste0('logit_race_age_impute_NSCLC_',j,'_res.csv'),row.names = F)
}  

####(2) prognosis
rm(list=ls())
library(data.table)
library(survival)
library(dplyr)

load('PRS_4_flip_tmb_sbs_age_race_impute.rda')

prs_list <- c('PRS1','PRS2','PRS3','PRS4')
type <- c('OS','DSS','DFI','PFI')

for (i in type) {
  res_out1 <- data.frame()
  for (j in prs_list) {
    os_name <- paste0(i,'.time')
    dat1 <- dat[,c(j,'age_use','gender_use','smoking_status2','ajcc_stage','race_use','type',i,os_name)]
    dat1$gender_use <- as.factor(dat1$gender_use)
    dat1$smoking_status2 <- as.factor(dat1$smoking_status2)
    dat1$type <- as.factor(dat1$type)
    dat1$ajcc_stage <- as.factor(dat1$ajcc_stage)
    dat1$race_use <- as.factor(dat1$race_use)
    formula <- as.formula(paste0('Surv(',i,'.time,',i,') ~  ',j,' + age_use + gender_use + smoking_status2 + type + ajcc_stage +race_use'))
    fit11 <- coxph(formula ,data=dat1)
    res <- summary(fit11)
    res <- as.data.frame(res$coefficients)
    res <- cbind(group=rownames(res),res)
    res_out1 <- rbind(res_out1,res)
  }
  write.csv(res_out1,file = paste0('race_age_impute_NSCLC_',i,'.csv'),row.names = F)
}


###(3) SBS
rm(list=ls())
library(data.table)
library(dplyr)
load('PRS_4_flip_tmb_sbs_age_race_impute.rda')

prs_list <- c('PRS1','PRS2','PRS3','PRS4')
signature_lsit <- paste0('Signature.',c(1:30))

for (i in prs_list) {
  res_signature <- data.frame()
  for (j in signature_lsit) {
    dat1 <- dat[,c(i,'age_use','gender_use','smoking_status2','race_use','type',j)]
    dat1$gender_use <- as.factor(dat1$gender_use)
    dat1$smoking_status2 <- as.factor(dat1$smoking_status2)
    dat1$type <- as.factor(dat1$type)
    dat1$race_use <- as.factor(dat1$race_use)
    formula <- paste0(j,' ~ ',i,'+ age_use + gender_use + smoking_status2 + type +race_use')
    fit11 <- glm(formula,data=dat1,family=quasipoisson(link = log))
    res <- summary(fit11)
    res <- as.data.frame(res$coefficients)
    res <- cbind(group=rownames(res),res)
    res$SBS <- j
    # res
    res_signature <- rbind(res_signature,res)
  }
  write.csv(res_signature,file = paste0('race_age_impute_NSCLC_',i,'.csv'),row.names = F)
}


###(4) TMB
rm(list=ls())
library(data.table)
library(dplyr)

load('PRS_4_flip_tmb_sbs_age_race_impute.rda')
dat <- subset(dat,!is.na(dat$total_perMB))

prs_list <- c('PRS1','PRS2','PRS3','PRS4')
res_signature <- data.frame()
for (i in prs_list) {
  dat1 <- dat[,c(i,'age_use','gender_use','smoking_status2','type','total_perMB','race_use')]
  dat1$gender_use <- as.factor(dat1$gender_use)
  dat1$smoking_status2 <- as.factor(dat1$smoking_status2)
  dat1$type <- as.factor(dat1$type)
  dat1$race_use <- as.factor(dat1$race_use)
  formula <- paste0('total_perMB ~ ',i,'+ age_use + gender_use + smoking_status2 +type +race_use')
  fit11 <- glm(formula,data=dat1,family=quasipoisson(link = log))
  res <- summary(fit11)
  res <- as.data.frame(res$coefficients)
  res <- cbind(group=rownames(res),res)
  res_signature <- rbind(res_signature,res)
}

write.csv(res_signature,file = paste0('~/TMB/race_age_impute_NSCLC_TMB.csv'),row.names = F)


###__________________________________________________________________________________________________________
###3. EGFR_EAS plot
rm(list = ls())
library(data.table)
library(rms)
library(ggplot2)
load('EGFR.rdata') 
load('PRS_4_fyt.rda') 

data$smoking_status2 <- ifelse(data$smoke==3,1,
                               ifelse(data$smoke%in%c(1,2) & data$packyear<30,2,
                                      ifelse(data$smoke%in%c(1,2) & data$packyear>=30,3,8)))
data_use <- data[,c('age','gender','smoking_status2','histology','EGFR','PRS2')]
data_use$study <- 1
data_f_use <- data_f[,c('age','gender','smoking_status2','histology','EGFR','PRS2')]
data_f_use$study <- 2

dat <- rbind(data_use,data_f_use)
data <- dat[,c('EGFR','PRS2','histology','study')]

data$PRS2.sd=data$PRS2/sd(data$PRS2)
data$histology <- as.factor(data$histology)
data$study <- as.factor(data$study)

pcubic=datadist(data)
options(datadist='pcubic')

fit_3=lrm(EGFR ~ rcs(PRS2.sd,3)+histology+study, data = data,x=T,y=T)

fit_4=lrm(EGFR ~ rcs(PRS2.sd,4)+histology+study, data = data,x=T,y=T)

fit_5=lrm(EGFR ~ rcs(PRS2.sd,5)+histology+study, data = data,x=T,y=T)

AIC(fit_3) #417.2436
AIC(fit_4) #417.8329
AIC(fit_5) #419.5759

#### y-hat =1 adjust ref, OR=1
Pre_HR <-rms::Predict(fit_3,PRS2.sd,fun=exp,type="predictions",ref.zero=T,conf.int = 0.95,digits=2)
ggplot(Pre_HR)
# y-hat=1.00000, PRS2.sd value 
a=subset(Pre_HR,round(Pre_HR$yhat,0)==1,select=c('yhat','PRS2.sd')) 
a$PRS2.sd[which(abs(a$yhat-1)==min(abs(a$yhat-1)))] #4.616512

ddist <- datadist(data)
ddist$limits$PRS2.sd[2] <- 4.616512
options(datadist="ddist")

fit_3=lrm(EGFR ~ rcs(PRS2.sd,3)+histology+study, data = data,x=T,y=T)

anova(fit_3) # P[overall] < 0.0202,  P[non-linear]= 0.3024
pred=Predict(fit_3,PRS2.sd,fun=exp,ref.zero=T)


PP0=ggplot(pred)+
  geom_ribbon(aes(PRS2.sd,ymin = lower, ymax = upper),linetype='blank',linewidth=0.7,fill="#BAE0FC")+theme_bw()+
  geom_line(aes(x=PRS2.sd ,y=yhat),linetype=1,linewidth=1,color="#4A9BCA")+
  geom_hline(yintercept=1, linetype=2,color="#AD002A99",linewidth=1)+
  labs(x="PRS2 (per sd increased)", y="OR of EGFR mutation")+
  theme(panel.background=element_blank(),panel.grid=element_blank(),
        axis.text=element_text(size=10,colour='black'),axis.title.x=element_text(size=12,
                                                                                 vjust=-3),axis.title.y=element_text(size=12,vjust=5, angle = 90),
        plot.caption=element_blank(),plot.margin = unit(c(1, 1, 1, 1),"cm"))

PP0
ggsave(file='EAS_EGFR_PRS2.pdf',width=6,height=4.5)

# OR
fit1 <- glm(EGFR ~ PRS2.sd +histology+study ,family=c("binomial"),data=data)
res <- summary(fit1)
beta <- res$coefficients[2,1]
se <- res$coefficients[2,2]
exp(beta) # 1.411928
exp(beta+1.96*se) # 1.82948
exp(beta-1.96*se) # 1.089676
