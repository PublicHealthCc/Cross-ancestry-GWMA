# Use an inverse variance-weighted fixed-effects model meta-analysis by METAL
###Give an example cross-ancestry LC GWMA

#bash_name: cross_ancestry_LC_GWMA
#!/bin/bash
# Input columns: 
MARKER MarkerName
ALLELE Effect_allele Other_allele
EFFECT BETA
STDERR SE
FREQ EAF
PVAL P
STRANDLABEL Strand
CUSTOMVARIABLE N
LABEL N AS N
# Metal Options:
SCHEME STDERR
WEIGHT N
USESTRAND ON
AVERAGEFREQ ON
MINMAXFREQ ON
VERBOSE OFF
GENOMICCONTROL ON
# Each GWAS-preMETA database was used for meta-analysis
PROCESS ~/GWAS_database/China/China_NSCLC_SNP_for_meta_nohete_preMETA
PROCESS ~/GWAS_database/BBJ/BBJ_SNP_for_meta_preMETA
PROCESS ~/GWAS_database/FLCCA/FLCCA_SNP_for_meta_preMETA
PROCESS ~/GWAS_database/ILCCO/ILCCO_SNP_for_meta_preMETA
PROCESS ~/GWAS_database/Finngen/Figgen_SNP_for_meta_preMETA

ANALYZE HETEROGENEITY
OUTFILE EAS_EUR_meta_fixed.tbl
QUIT