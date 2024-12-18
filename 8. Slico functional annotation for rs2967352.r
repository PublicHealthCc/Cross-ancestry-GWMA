#####Slico functional annotation for rs2967352(16q23.3, MPHOSPH6)#####################################
####Step1. rs2967352 high LD variants
plink --noweb --bfile ~/1KG_EAS_EUR/chr_all_1kgv3_2015_EAS_EUR \
--ld-snp rs2967352 --ld-window 99999 --ld-window-kb 1000 \
--ld-window-r2 0.6 --out rs2967352_1mb_r2_0.6 --r2

R
library(data.table)
library(dplyr)
library(tidyr)
SNP<-fread('~/rs2967352_1mb_r2_0.6.ld')
SNP<-SNP[,c(4,5)]
SNP$CHR_B<-as.character(SNP$CHR_B)
SNP$CHR_B<-'chr16'
SNP$BP<-SNP$BP_B
SNP<-SNP[!SNP$BP =='82196676',]
fwrite(SNP,'~/rs2967352.bed',sep = '\t',col.names = F)

#######Step2. Promoters, enhancers, TFBS and DHS based on A549 cell line
bedtools intersect -a rs2967352.bed -b /data1/Trans_meta/func_annotation/A549/A549_H3k04me3.updated.bed -wa -wb > rs2967352_A549_H3k04me3.bed
bedtools intersect -a rs2967352.bed -b /data1/Trans_meta/func_annotation/A549/A549_H3k09ac.updated.bed -wa -wb > rs2967352_A549_H3k09ac.bed
bedtools intersect -a rs2967352.bed -b /data1/Trans_meta/func_annotation/A549/A549_H3k04me1.updated.bed -wa -wb > rs2967352_A549_H3k04me1.bed
bedtools intersect -a rs2967352.bed -b /data1/Trans_meta/func_annotation/A549/A549_H3k27ac.updated.bed -wa -wb > rs2967352_A549_H3k27ac.bed

bedtools intersect -a rs2967352.bed -b /data1/Trans_meta/func_annotation/A549/ENCODE_A549_35TFBS_ALL.bed -wa -wb > rs2967352_A549_TFBS.bed

bedtools intersect -a rs2967352.bed -b /data1/Trans_meta/func_annotation/A549/wgEncodeAwgDnaseUwHpfUniPk.narrowPeak -wa -wb > rs2967352_DHS1.bed
bedtools intersect -a rs2967352.bed -b /data1/Trans_meta/func_annotation/A549/wgEncodeAwgDnaseUwNhlfUniPk.narrowPeak -wa -wb > rs2967352_DHS2.bed
bedtools intersect -a rs2967352.bed -b /data1/Trans_meta/func_annotation/A549/wgEncodeAwgDnaseUwAg04450UniPk.narrowPeak -wa -wb > rs2967352_DHS3.bed
bedtools intersect -a rs2967352.bed -b /data1/Trans_meta/func_annotation/A549/wgEncodeAwgDnaseUwWi38UniPk.narrowPeak -wa -wb > rs2967352_DHS4.bed
bedtools intersect -a rs2967352.bed -b /data1/Trans_meta/func_annotation/A549/wgEncodeAwgDnaseUwdukeA549UniPk.narrowPeak -wa -wb > rs2967352_DHS5.bed

####Step3. GM12878 and NHLF cell line
R
library(data.table)
library(dplyr)
library(tidyr)
SNP<-fread('~/rs2967352_1mb_r2_0.6.ld')
bim<-fread('~/chr_all_1kgv3_2015_EAS_EUR.bim')
SNP<-SNP[,c(4:6)]
SNP<-SNP[!SNP$BP_B =='82196676',]
bim<-bim[,c(2,5,6)]
colnames(bim)[1]<-'SNP_B'
data<-left_join(SNP,bim,by='SNP_B')
data$BP<-data$BP_B
data<-data[,c(1,2,6,4,5)]
fwrite(data,'~/rs2967352.avinput',sep = '\t',col.names = F)

##GM12878
wget https://hgdownload2.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeBroadHmmGm12878HMM.txt.gz
gzip -d wgEncodeBroadHmmGm12878HMM.txt.gz
perl ~/annovar/annotate_variation.pl -regionanno -buildver hg19 -dbtype wgEncodeBroadHmmGm12878HMM rs2967352.avinput ~/annovar_downdb
###NHLF
wget https://hgdownload2.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeBroadHmmNhlfHMM.txt.gz
gzip -d wgEncodeBroadHmmNhlfHMM.txt.gz
perl ~/annovar/annotate_variation.pl -regionanno -buildver hg19 -dbtype wgEncodeBroadHmmNhlfHMM rs2967352.avinput ~/annovar_downdb

#####Step4. 3DSNP https://omic.tech/3dsnpv2/
#####Step5. CADD: https://cadd.bihealth.org/score
###vcf 
plink  --bfile ~/1KG_EAS_EUR/chr_all_1kgv3_2015_EAS_EUR \
--extract ~/rs2967352_variants.txt  --recode vcf-iid --out rs2967352_anno_vcf

#####Step6. Regulome DB: https://regulome.stanford.edu/regulome-search  
#####Step7. GWAMA-Snp nexus: https://www.snp-nexus.org/v4/
