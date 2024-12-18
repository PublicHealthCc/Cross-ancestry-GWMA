# Use an inverse variance-weighted fixed-effects model meta-analysis by METAL
###1. EAS-specific LC GWMA

#bash_name: EAS_specific_LC_GWMA
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

ANALYZE HETEROGENEITY
OUTFILE EAS_meta_fixed.tbl
QUIT

##Operation metal
~/generic-metal/metal EAS_specific_LC_GWMA.bash

###2. EAS-specific LUAD GWMA
######BBJ has no subtype information
#bash_name: EAS_specific_LUAD_GWMA
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
PROCESS ~/GWAS_database/China/China_AD_SNP_for_meta_nohete_preMETA
PROCESS ~/GWAS_database/FLCCA/FLCCA_AD_SNP_for_meta_preMETA

ANALYZE HETEROGENEITY
OUTFILE EAS_AD_meta_fixed.tbl
QUIT

##Operation metal
~/generic-metal/metal EAS_specific_LUAD_GWMA.bash

###3. EAS-specific LUSC GWMA
######BBJ has no subtype information
#bash_name: EAS_specific_LUSC_GWMA
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
PROCESS ~/GWAS_database/China/China_SC_SNP_for_meta_nohete_preMETA
PROCESS ~/GWAS_database/FLCCA/FLCCA_SC_SNP_for_meta_preMETA

ANALYZE HETEROGENEITY
OUTFILE EAS_SC_meta_fixed.tbl
QUIT

##Operation metal
~/generic-metal/metal EAS_specific_LUSC_GWMA.bash

###_________________________________________________________________________________________________________________________________
###4. EUR-specific LC GWMA

#bash_name: EUR_specific_LC_GWMA
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
PROCESS ~/GWAS_database/ILCCO/ILCCO_SNP_for_meta_preMETA
PROCESS ~/GWAS_database/Finngen/Figgen_SNP_for_meta_preMETA

ANALYZE HETEROGENEITY
OUTFILE EUR_meta_fixed.tbl
QUIT

##Operation metal
~/generic-metal/metal EUR_specific_LC_GWMA.bash


###5. EUR-specific LUAD GWMA

#bash_name: EUR_specific_LUAD_GWMA
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
PROCESS ~/GWAS_database/ILCCO/ILCCO_AD_SNP_for_meta_preMETA
PROCESS ~/GWAS_database/Finngen/Figgen_AD_SNP_for_meta_preMETA

ANALYZE HETEROGENEITY
OUTFILE EUR_AD_meta_fixed.tbl
QUIT

##Operation metal
~/generic-metal/metal EUR_specific_LUAD_GWMA.bash

###6. EUR-specific LUSC GWMA

#bash_name: EUR_specific_LUSC_GWMA
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
PROCESS ~/GWAS_database/ILCCO/ILCCO_SC_SNP_for_meta_preMETA
PROCESS ~/GWAS_database/Finngen/Figgen_SC_SNP_for_meta_preMETA

ANALYZE HETEROGENEITY
OUTFILE EUR_SC_meta_fixed.tbl
QUIT

##Operation metal
~/generic-metal/metal EUR_specific_LUSC_GWMA.bash

###_________________________________________________________________________________________________________________________________
###7. cross-ancestry LC GWMA

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

##Operation metal
~/generic-metal/metal cross_ancestry_LC_GWMA.bash

###8. cross-ancestry LUAD GWMA

#bash_name: cross_ancestry_LUAD_GWMA
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
PROCESS ~/GWAS_database/China/China_AD_SNP_for_meta_nohete_preMETA
PROCESS ~/GWAS_database/FLCCA/FLCCA_AD_SNP_for_meta_preMETA
PROCESS ~/GWAS_database/ILCCO/ILCCO_AD_SNP_for_meta_preMETA
PROCESS ~/GWAS_database/Finngen/Figgen_AD_SNP_for_meta_preMETA

ANALYZE HETEROGENEITY
OUTFILE EAS_EUR_AD_meta_fixed.tbl
QUIT

##Operation metal
~/generic-metal/metal cross_ancestry_LUAD_GWMA.bash

###9. cross-ancestry LUSC GWMA

#bash_name: cross_ancestry_LUSC_GWMA
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
PROCESS ~/GWAS_database/China/China_SC_SNP_for_meta_nohete_preMETA
PROCESS ~/GWAS_database/FLCCA/FLCCA_SC_SNP_for_meta_preMETA
PROCESS ~/GWAS_database/ILCCO/ILCCO_SC_SNP_for_meta_preMETA
PROCESS ~/GWAS_database/Finngen/Figgen_SC_SNP_for_meta_preMETA

ANALYZE HETEROGENEITY
OUTFILE EAS_EUR_SC_meta_fixed.tbl
QUIT

##Operation metal
~/generic-metal/metal cross_ancestry_LUSC_GWMA.bash
