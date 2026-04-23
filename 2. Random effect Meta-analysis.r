## ---------------------------------------------------
## Meta-analysis of random effect models for highly heterogeneous loci
## Give an example of LC
## 1. LC
R
library(data.table)
library(tidyr)
library(dplyr)
data<-fread('./EAS_EUR_meta_fixed.txt')
het<-data[data$HetISq>=75,]
BBJ<-fread('./GWAS_database/BBJ/BBJ_SNP_for_meta_preMETA')
NJMU<-fread('./GWAS_database/China/China_NSCLC_SNP_for_meta_nohete_preMETA')
FLCCA<-fread('./GWAS_database/FLCCA/FLCCA_SNP_for_meta_preMETA')
FING<-fread('./GWAS_database/Finngen/Figgen_SNP_for_meta_preMETA')
ILLCO<-fread('./GWAS_database/ILCCO/ILCCO_SNP_for_meta_preMETA')

BBJ_het <- BBJ[BBJ$MarkerName %in% het$MarkerName, ]
NJMU_het <- NJMU[NJMU$MarkerName %in% het$MarkerName, ]
FLCCA_het <- FLCCA[FLCCA$MarkerName %in% het$MarkerName, ]
FING_het <- FING[FING$MarkerName %in% het$MarkerName, ]
ILCCO_het <- ILLCO[ILLCO$MarkerName %in% het$MarkerName, ]

colnames(BBJ_het) <- toupper(colnames(BBJ_het))
colnames(NJMU_het) <- toupper(colnames(NJMU_het))
colnames(FLCCA_het) <- toupper(colnames(FLCCA_het))
colnames(FING_het) <- toupper(colnames(FING_het))
colnames(ILCCO_het) <- toupper(colnames(ILCCO_het))

## prepare input file for GWAMA
fwrite(BBJ_het,'./1_hete_GWAS/BBJ_het.txt',quote=F,row.names=F,col.names=T,sep='\t')
fwrite(NJMU_het,'./1_hete_GWAS/NJMU_het.txt',quote=F,row.names=F,col.names=T,sep='\t')
fwrite(FLCCA_het,'./1_hete_GWAS/FLCCA_het.txt',quote=F,row.names=F,col.names=T,sep='\t')
fwrite(FING_het,'./1_hete_GWAS/FinnGen_het.txt',quote=F,row.names=F,col.names=T,sep='\t')
fwrite(ILCCO_het,'./1_hete_GWAS/ILCCO_het.txt',quote=F,row.names=F,col.names=T,sep='\t')

## GWAMA_meta
./software/GWAMA/GWAMA -i ./1_hete_GWAS/Trans_LC_het.in \
-o ./1_hete_GWAS/LC_results/Trans_LC_het_result \
-qt \
--random \
--name_marker MARKERNAME \
--name_ea EFFECT_ALLELE \
--name_nea OTHER_ALLELE \
--name_eaf EAF \
--name_beta BETA \
--name_se SE 

## step 2. Integrate the results of the highly heterogeneous sites RE with those of the fixed-effect models of other sites for clumping.
