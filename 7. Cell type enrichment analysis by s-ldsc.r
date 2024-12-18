##Cell type enrichment analysis by s-LDSC
##Step1. Cross-ancestry GWMA LDSC panel
plink --bfile ~/1000Genome/chr_all_1kgv3_2015_EAS_EUR \
--make-bed \
--cm-map genetic_map_b37/genetic_map_chr@_combined_b37.txt \
--out 1000G_EAS_EUR_with_cms/1000G_EAS_EUR_with_cms

R
library(data.table)
dat <- fread('1000G_EAS_EUR_with_cms/1000G_EAS_EUR_with_cms.bim')
dat1 <- subset(dat,nchar(dat$V5)==1 & nchar(dat$V6)==1)
write.table(dat1$V2,file = '1000G_EAS_EUR_Phase3_plink/extract_list',row.names = F,col.names = F,sep = '\t',quote = F)

plink --bfile 1000G_EAS_EUR_with_cms/1000G_EAS_EUR_with_cms \
--extract 1000G_EAS_EUR_Phase3_plink/extract_list \
--make-bed \
--out 1000G_EAS_EUR_Phase3_plink/chr_all_1kgv3_2015_EAS_EUR_single_allele

for i in {1..22} ;
do 
plink --bfile 1000G_EAS_EUR_Phase3_plink/chr_all_1kgv3_2015_EAS_EUR_single_allele \
--chr $i \
--make-bed \
--keep-allele-order \
--maf 0.005 \
--geno 0.05 \
--hwe 1e-06 \
--out 1000G_EAS_EUR_Phase3_plink/1000G_EAS_EUR_QC.chr$i
done 


for i in {1..22} ;
do 
plink --bfile 1000G_EAS_EUR_Phase3_plink/1000G_EAS_EUR_QC.chr$i \
--freq \
--out 1000G_EAS_EUR_Phase3_frq/1000G_EAS_EUR_QC.chr$i
done

# LD panel
plink --bfile /data/winD/Public/1000Genome/1kg_new/chr_all_1kgv3_2015_EAS_EUR \
--make-bed \
--cm-map genetic_map_b37/genetic_map_chr@_combined_b37.txt \
--out 1000G_EAS_EUR_with_cms/1000G_EAS_EUR_with_cms

# extract SNP in Hapmap
plink --bfile 1000G_EAS_EUR_with_cms/1000G_EAS_EUR_with_cms \
--extract /data/home/jinchen/ldsc/EAS/print_snps.txt \
--keep-allele-order \
--make-bed \
--out 1000G_EAS_EUR_with_cms/1000G_EAS_EUR_with_cms_hm3

## disassemble chr
for chr in {1..22} ;
do 
plink --bfile 1000G_EAS_EUR_with_cms/1000G_EAS_EUR_with_cms_hm3 \
--keep-allele-order \
--chr $chr \
--make-bed \
--out 1000G_EAS_EUR_MHC/weights.hm3_noMHC.$chr
done

# exclude MHC
plink --bfile 1000G_EAS_EUR_with_cms/1000G_EAS_EUR_with_cms_hm3 \
--chr 6 \
--from-mb 25 \
--to-mb 34 \
--write-snplist \
--out 1000G_EAS_EUR_MHC/hla

plink --bfile 1000G_EAS_EUR_MHC/weights.hm3_noMHC.6 \
--keep-allele-order \
--chr 6 \
--exclude 1000G_EAS_EUR_MHC/hla.snplist \
--make-bed \
--out 1000G_EAS_EUR_MHC/weights.hm3_noMHC.6

# Calculate LD-scores
for chr in {1..22} ;
do 
python ~/ldsc/ldsc.py \
--bfile 1000G_EAS_EUR_MHC/weights.hm3_noMHC.$chr \
--l2 \
--ld-wind-cm 1 \
--out 1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.$chr 
done


###Step2. Construct LD reference template based on Korean single cell data
R
options(stringsAsFactors=F)
library(data.table)
rm(list=ls())
gen <- read.csv('Korea_Cell_subtype_DEG.csv')
names(gen)[6] <- 'group'
gen$group <- ifelse(gen$group=='Alveolar Mac','Alveolar_Mac',gen$group)
gen$group <- ifelse(gen$group=='B lymphocytes','B_lymphocytes',gen$group)
gen$group <- ifelse(gen$group=='Endothelial cells','Endothelial_cells',gen$group)
gen$group <- ifelse(gen$group=='MAST cells','MAST_cells',gen$group)
gen$group <- ifelse(gen$group=='mo-Mac','mo_Mac',gen$group)
gen$group <- ifelse(gen$group=='NK cells','NK_cells',gen$group)
gen$group <- ifelse(gen$group=='Pleural Mac','Pleural_Mac',gen$group)
gen$group <- ifelse(gen$group=='T lymphocytes','T_lymphocytes',gen$group)
ref <- read.csv('../fetal_lung/ref/gene_change_for_hg19.csv')
res <- data.frame()
for (i in 1:length(names(table(gen$group)))) {
  res1 <- subset(gen,gen$group==names(table(gen$group))[i])
  res1 <- subset(res1,p_val_adj<0.05)
  res1 <- res1[order(res1$avg_log2FC,decreasing = T),]
  names(res1)[7] <- 'gene_name'
  res1 <- dplyr::left_join(res1,ref,by = 'gene_name')
  write.table(res1$gene_id,file = paste0('scRNA_single_',names(table(gen$group))[i],'.txt'),row.names = F,col.names = F,sep = '\t',quote = F)
  res <- rbind(res,res1)
}

write.table(res,file = 'per_scRNA_diff_genes.txt',row.names = F,col.names = T,sep = '\t',quote = F)

# LD panel
file_dir='~/s_ldsc/scRNA/Korea_new'

cat config.scRNA |while read cellType;
do

arr=($cellType)
gene_set=${arr[0]}
cellName=${arr[1]}

for i in {1..22}
do
python ~/ldsc/make_annot.py \
--gene-set-file $gene_set \
--gene-coord-file $file_dir/fetal_lung/ref/gene_coordinates_hg19.txt \
--windowsize 100000 \
--bimfile /data1/qtl/ldsc/plink_file/EAS/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC.${i}.bim \
--annot-file $file_dir/$cellName/tmp_$cellName.${i}.annot

gzip $file_dir/$cellName/tmp_$cellName.${i}.annot
gunzip -c $file_dir/$cellName/tmp_$cellName.${i}.annot.gz |awk 'BEGIN{OFS="\t"}{if(NR==1){print "ALL",$1}else{print 1,$1}}' |gzip -c > $file_dir/$cellName/$cellName.${i}.annot.gz

python ~/ldsc/ldsc.py \
--l2 --ld-wind-cm 1 \
--bfile ~/EAS/1000G_Phase3_EAS_plinkfiles/1000G.EAS.QC.${i} \
--annot $file_dir/$cellName/$cellName.${i}.annot.gz \
--thin-annot \
--out $file_dir/$cellName/$cellName.${i} \
--print-snps ~/ldsc/hapmap3_snps/print_snps.txt 

done
done


##Step2. Cell type enrichment
cat sumstats_list |while read cellType;
do
python ~/ldsc/munge_sumstats.py \
--sumstats ~/${cellType} \
--N-col N \
--chunksize 500000 \
--out ./${cellType}
done

##Give an example: Cross-ancestry LC GWMA
subtype='LC'
gwas_index='Trans_LC_meta_sumstats.txt' 
people='EAS_EUR'
kg_ref_dir='~/ldsc/plink_file/EAS_EUR' 
ref_dir='~/s_ldsc/scRNA/Korea_new_EAS_EUR' 
gwas_dir='~/s-ldsc/2024_07_22/sumstats' 
output_dir='~/s-ldsc/2024_07_22/Korea' 

cat ${ref_dir}/1 |while read cellType;
do
python ~/ldsc/ldsc.py \
--h2 ${gwas_dir}/${gwas_index}.sumstats.gz \
--ref-ld-chr ${ref_dir}/${cellType}/${cellType}. \
--frqfile-chr ${kg_ref_dir}/1000G_EAS_EUR_Phase3_frq/1000G_EAS_EUR_QC.chr \
--print-coefficients \
--overlap-annot \
--w-ld-chr ${kg_ref_dir}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--out ${output_dir}/${people}/${subtype}/s_ldsc_${cellType}
done

