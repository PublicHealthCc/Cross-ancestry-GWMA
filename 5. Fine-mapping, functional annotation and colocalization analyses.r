##______________________________________________________________________________________________________________________
########1. Fine mapping by Susie
###Give an example: Cross-ancestry LC GWMA
###Step1. Prepare file and PolyFun_susie
R
library(data.table)
list <- fread('../lead_snp_list.txt')
people <- 'EAS_EUR'
type <- 'LC'
list <- subset(list,subtype==type)
meta_file_path <- '/data1/Trans_meta/meta_dataset_fixed/Trans_meta/EAS_EUR_meta_fixed.txt'
meta_sample_size <- 539860
chunk_size <- 500000
soft_dir <- '/home/jinchen/software/polyfun-master'
res_dir <- '/data1/qtl/test/eQTL_for_ZCC/susie/PolyFun'
meta <- fread(meta_file_path)
for (i in 1:nrow(list)) {
  system(paste0('plink --noweb --bfile /data/Public/1000Genome/1kg_new/chr_all_1kgv3_2015_',people,' --snp ',list$rsid[i],' --window 500 --make-bed --out ',people,'/',type,'/',list$rsid[i],'/',list$rsid[i]))
  meta_dat <- subset(meta,CHR == list$chr[i] & BP<list$pos[i]+chunk_size/2 & BP>list$pos[i]-chunk_size/2)
  meta_dat <- meta_dat[,c('MarkerName','Allele1','Allele2','Freq1','Effect','StdErr','P','CHR','BP')]
  names(meta_dat) <- c('MarkerName','ALLELE1','ALLELE0','A1FREQ','BETA','SE','P','CHR','BP')
  write.table(meta_dat,file = paste0(people,'/',type,'/',list$rsid[i],'/',list$rsid[i],'_sumstats'),row.names = F,col.names = T,sep = '\t',quote = F)
  system(paste0('gzip ',people,'/',type,'/',list$rsid[i],'/',list$rsid[i],'_sumstats'))
  system(paste0('python ',soft_dir,'/munge_polyfun_sumstats.py --sumstats ',res_dir,'/',people,'/',type,'/',list$rsid[i],'/',list$rsid[i],'_sumstats.gz --n ',meta_sample_size,' --out ',res_dir,'/',people,'/',type,'/',list$rsid[i],'/',list$rsid[i],'_sumstats_munged.parquet --min-info 0.3 --min-maf 0.01 --keep-hla'))
}

R
library(data.table)
rm(list = ls())
list <- fread('../lead_snp_list.txt')
people <- 'EAS_EUR'
type <- 'LC'
list <- subset(list,subtype==type)
meta_sample_size <- 539860
chunk_size <- 500000
L=5
soft_dir <- '/home/jinchen/software/polyfun-master'
res_dir <- '/data1/qtl/test/eQTL_for_ZCC/susie/PolyFun'
list$start <- list$pos-chunk_size/2
list$end <- list$pos+chunk_size/2
for (i in 1:nrow(list)) {
  system(paste0('python ',soft_dir,'/finemapper.py --geno ',res_dir,'/',people,'/',type,'/',list$rsid[i],'/',list$rsid[i],' --sumstats ',res_dir,'/',people,'/',type,'/',list$rsid[i],'/',list$rsid[i],'_sumstats_munged.parquet --n ',meta_sample_size,' --chr ',list$chr[i],' --start ',list$start[i],' --end ',list$end[i],' --method susie --allow-missing --non-funct --max-num-causal ',L,' --out ',res_dir,'/',people,'/',type,'/',list$rsid[i],'/',list$rsid[i],'_susie.gz'))
}

###Step2. Functional annotation of the potential casual variants of PolyFun_susie
# Following methods：
# 1.Annovar 
# 2.SNP nexus：https://www.snp-nexus.org/v4/
# 3.CADD：https://cadd.gs.washington.edu/score
# 4.PolyPhen :http://genetics.bwh.harvard.edu/pph2/


##______________________________________________________________________________________________________________________
####2. eQTL-eGene overlap##################
###Give an example: NJMU_EAS_338_lung tissue eQTL
plink --bfile ~/1000Genome/chr_all_1kgv3_2015_EAS_EUR \
  --ld-snp-list ~/36loci_500k_0.6.txt \
  --ld-window 9999999 \
  --ld-window-kb 500 \
  --ld-window-r2 0.6 \
  --r2 \
  --out Trans_36loci_500k_0.6 \

###########Extract GWAS high LD points
R
library(dplyr)
library(data.table)
library(tidyr)
gwas<-fread('~/Trans_36loci_500k_0.6.ld')
colnames(gwas)[3]<-'lead_snp'
colnames(gwas)[6]<-'rsid'
bim<-fread('~/1000Genome/chr_all_1kgv3_2015_EAS_EUR.bim')
colnames(bim)[2]<-'rsid'
data<-left_join(gwas,bim,by='rsid')

data<-data %>% mutate(Base_pair=paste(V5,V6,sep=':'))
data$Base_unified <- ifelse((data$Base_pair %in% c("T:G", "G:T", "C:A")),"A:C", 
                                  ifelse((data$Base_pair %in% c("T:C", "C:T", "G:A")), "A:G",
                                         ifelse(data$Base_pair == "T:A","A:T",
                                                ifelse(data$Base_pair == "G:C","C:G",data$Base_pair))))
data<-data %>% mutate(variant_id=paste(CHR_B,BP_B,Base_unified,sep=':'))
data<-data[,c('variant_id','lead_snp','rsid','R2')]

data <- data %>% distinct(rsid, .keep_all = TRUE)
fwrite(data,'~/Trans_36loci_0.6_overlap.txt',sep = '\t',col.names=T)

######Extract overlap eqtl
eqtl<-fread('~/338sample_eqtl_coloc.csv.gz')
gwas<-fread('~/Trans_36loci_0.6_overlap.txt')
eqtl$pfdr <- p.adjust(eqtl$pval_nominal, method = "fdr")
eqtl<-eqtl[eqtl$pfdr < 0.05,]
colnames(eqtl)[2]<-'variant_id'
overlap<-inner_join(gwas,eqtl,by = 'variant_id')
overlap$symbol=sapply(strsplit(overlap$gene_id,"[.]"),"[",1)
fwrite(overlap,'~/GWAS_338eqtl_overlap_all.csv',sep=',',col.names=T)
overlap1 <- overlap %>%
  group_by(lead_snp) %>%
  arrange(pfdr) %>%
  distinct(symbol, .keep_all = TRUE) %>%
  ungroup()
fwrite(overlap1,'~/GWAS_338eqtl_overlap_gene.csv',col.names=T)



#########________________________________________________________________________________________________________________________________________
#####3. Colocalization analyses################################
# Give an example: NJMU_EAS_338_lung tissue eQTL
qtl_coloc <- function(gene_id, coloc_file, sample_size) {
  eqtl_coloc <- list(
    snp = coloc_file$SNP,
    position = coloc_file$BP,
    beta = coloc_file$beta_eqtl,
    varbeta = coloc_file$varbeta_eqtl,
    type = "quant",
    N = as.numeric(sample_size),
    MAF = coloc_file$maf
  )
  gwas_coloc <- list(
    snp = coloc_file$SNP,
    position = coloc_file$BP,
    beta = coloc_file$beta_gwas,
    varbeta = coloc_file$varbeta1_gwas,
    type = "cc"
  )
  result <- coloc.abf(dataset1 = eqtl_coloc, dataset2 = gwas_coloc)
  res <- t(as.data.frame(result$summary))
  rownames(res) <- gene_id
  res <- as.data.frame(res)
  return(res)
}

res_40all <- tibble()
sample_size <- 338
  
eqtl_geno <- fread('~/338sample_eqtl_coloc.csv.gz')
gwas_files <- list.files(path = "~/29loci_GWAS", pattern = "*.txt", full.names = TRUE)
  
  for (file_path in gwas_files) {
    gwas_geno <- fread(file_path)
    coloc_file <- inner_join(eqtl_geno, gwas_geno, by = "SNP")
    
    coloc_file <- coloc_file %>%
      mutate(varbeta1_gwas = varbeta_gwas^2) %>%
      arrange(gene_id) %>%
      drop_na()
    
    egene_list <- coloc_file %>%
      filter(pval_nominal < 0.0001) %>%
      distinct(gene_id) %>%
      pull(gene_id)
    
    for (gene in egene_list) {
      sub <- coloc_file %>%
        filter(gene_id == gene) %>%
        distinct(SNP, .keep_all = TRUE)
      
      if (nrow(sub) > 0) {
        res <- qtl_coloc(gene, sub, sample_size)
        res <- cbind(variant_id = rownames(res), res)
        res$rsid <- gsub("\\.txt$", "", basename(file_path))
        res_40all <- bind_rows(res_40all, res)
      }
    }
  }
  
  fwrite(res_40all, '~/Step3_NJMU_EAS_lung_tissue.csv', sep = ',', row.names = FALSE) 



####_________________________________________________________________________________________________
####4. Colocalization plot
# Give an example: rs2967352(MPHOSPH6), NJMU_lung_eqtl
R
library(data.table)
library(dplyr)
library(tidyr)
library(locuscomparer)
eqtl_geno <- fread('~/338sample_eqtl_coloc.csv.gz')
gwas_geno <- fread('~/rs2967352.txt')
coloc_plot <- inner_join(eqtl_geno, gwas_geno, by = "SNP")
##MPHOSPH6
coloc_plot1<-subset(coloc_plot,coloc_plot$gene_id=='ENSG00000135698.5') 
gwas<-coloc_plot1 %>% select(rsID,P)
eqtl<-coloc_plot1 %>% select(rsID,pval_nominal)
colnames(gwas)[2]<- 'pval'
colnames(eqtl)[2]<- 'pval'
fwrite(gwas,"gwas_eas.tsv",col.names = T,row.names = F,sep="\t",quote = F)
fwrite(eqtl,"eqtl_eas.tsv",col.names = T,row.names = F,sep="\t",quote = F)
plot<- locuscompare(in_fn2 ="gwas_eas.tsv", in_fn1 = "eqtl_eas.tsv",
             marker_col2 = "rsID", pval_col2 = "pval", title2 = "Cross-ancestry GWMA of LUAD",
             marker_col1= "rsID", pval_col1 = "pval", title1 = "NJMU lung tissue",snp = 'rs2967352',combine = FALSE, population = 'EAS')
ggsave("eqtl_LUAD_coloc_eas_compare.pdf", plot =plot$locuscompare, width = 6, height = 5, dpi = 600)
ggsave("eqtl_LUAD_coloc_eas_locus1.pdf", plot =plot$locuszoom1, width = 5, height = 2.7, dpi = 600)
ggsave("eqtl_LUAD_coloc_eas_locus2.pdf", plot =plot$locuszoom2, width = 5, height = 2.7, dpi = 600)


####_________________________________________________________________________________________________
####5. eqtl-gene expression plot
# Give an example: rs2967363-MPHOSPH6
R 
library("ggplot2")
library("ggbeeswarm")
library("ggsci")
library(data.table)
library(tidyverse)
library(ggpubr)
raw_363 <- fread('338_rs2967363.raw')
colnames(raw_363)[7] <- 'C'
expr <- fread('338exp_for_MPHOSPH6.txt')
id <- fread('338_rs2967363.nosex')
expr <- t(expr)
expr <- as.data.frame(expr)
expr <- rownames_to_column(expr)
expr <- expr[-1,]
colnames(expr) <- c('id','MPHOSPH6')
colnames(raw_363)[2] <- 'id'
rs363 <- inner_join(raw_363,expr,by='id')

cov<-fread('338LungCombine.combined_covariates_5pca.txt',header=T)
cov<-as.data.frame(t(cov))
cov <- rownames_to_column(cov)
colnames(cov) <- cov[1,]
cov <- cov[-1,]
data=merge(rs363,cov,by.x=c("id"),by.y=c("ID"))
data$MPHOSPH6 <- as.numeric(data$MPHOSPH6)
model <- glm(MPHOSPH6~C+pc1+pc2+pc3+pc4+pc5+InferredCov1+InferredCov2
             +InferredCov3+InferredCov4+InferredCov5+InferredCov6+InferredCov7+InferredCov8
             +InferredCov9+InferredCov10+InferredCov11+InferredCov12+InferredCov13+InferredCov14
             +InferredCov15+InferredCov16+InferredCov17+InferredCov18+InferredCov19+InferredCov20+InferredCov21
             +InferredCov22+InferredCov23+InferredCov24+InferredCov25+InferredCov26+InferredCov27+InferredCov28
             +InferredCov29+InferredCov30+InferredCov31+InferredCov32+InferredCov33+InferredCov34+InferredCov35
             +InferredCov36+InferredCov37+InferredCov38+InferredCov38+InferredCov40+InferredCov41+InferredCov42
             +InferredCov43+InferredCov44+InferredCov45+sex+age+smoking_status+batch,data=data)

###############plot
data$predict_MPHOSPH6 <- predict(model)
data$group <- ifelse(data$C == 0, 'GG', ifelse(data$C == 1, 'CG', 'CC'))
data$group <- factor(data$group, levels = c('CC', 'CG', 'GG'))

data$MPHOSPH6_plot <- data$MPHOSPH6 - data$predict_MPHOSPH6

P2 <- ggplot(data, aes(x = group, y = MPHOSPH6_plot, fill = group)) + 
  geom_violin(trim = TRUE, color = "white", alpha = 0.5) + 
  stat_compare_means(method = 'anova', label.y = 4, label.x = 0.7, label = 'p') +
  geom_boxplot(width = 0.25, position = position_dodge(0.9), alpha = 1, fill = '#555F66', outlier.shape = NA, color = "white") + 
  scale_fill_manual(values = c("#56B4E9", "#56B4E9", "#56B4E9")) +
  theme_bw() +  
  theme(panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "gray"),  
        axis.text = element_text(color = "gray", margin = margin(0.5, 0.5)),  
        axis.title = element_text(color = "gray", size = 12, face = "bold"), 
        axis.ticks = element_line(color = "gray", size = 1),  
        axis.ticks.length = unit(0.3, "cm"), 
        axis.line.x = element_line(size = 1), 
        axis.line.y = element_line(size = 1),  
        axis.line.x.bottom = element_line(size = 1),  
        axis.line.y.left = element_line(size = 1),  
        axis.line.x.top = element_blank(),  
        axis.line.y.right = element_blank()) +
  guides(fill = "none")
#P<2.00E-16
P2 <- P2 + ylim(min(-3.5E-14), max(3.5E-14))
P2
ggsave("rs2967363_eqtl_cov.pdf", P2, width = 8, height = 5, dpi = 600)



