#____________________________________________________________________________________________
#######Figure2A
R
library(data.table)
library(tidyverse)
library("ggsci")
library(scales)
data<-fread('Cross-ancestry_comparison_plot.csv')
data$Effect_EAS_IL=data$Effect_EAS-1.96*data$SE_EAS
data$Effect_EAS_UL=data$Effect_EAS+1.96*data$SE_EAS

data$Effect_EUR_IL=data$Effect_EUR-1.96*data$SE_EUR
data$Effect_EUR_UL=data$Effect_EUR+1.96*data$SE_EUR
data_complete <- na.omit(data)
# glm
model <- glm(Effect_EUR ~ Effect_EAS, data = data_complete)
summary(model)

#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -0.02090    0.01347  -1.552    0.132    
#Effect_EAS   0.64408    0.10590   6.082 1.27e-06 ***
#AIC: -71.593

##r
cor(data_complete$Effect_EAS,data_complete$Effect_EUR)
#[1] 0.748705

#######Group
data$group[data$P_EAS < 5E-8 & data$P_EUR >= 5E-8]<-'EAS'
data$group[data$P_EAS >= 5E-8 & data$P_EUR < 5E-8]<-'EUR'
data$group[data$P_EAS < 5E-8 & data$P_EUR < 5E-8]<-'EAS_EUR'
data$group[data$P_EAS >= 5E-8 & data$P_EUR >= 5E-8]<-'Others'

#######Plot
p <- ggplot(data) + 
    geom_point(aes(x = Effect_EAS, y = Effect_EUR, fill = group, color = group ,fill = group, shape = group), size = 2) + 
    theme(
        legend.position = "right",  
        panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA),   
        axis.line = element_line(colour = "black")
    ) +
    scale_x_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, 0.1), labels = number_format(accuracy = 0.1)) +
    scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, 0.1), labels = number_format(accuracy = 0.1)) +
    geom_segment(aes(y = Effect_EUR_IL, yend = Effect_EUR_UL, x = Effect_EAS, color = group)) + 
    geom_segment(aes(y = Effect_EUR, x = Effect_EAS_IL, xend = Effect_EAS_UL, color = group)) +
    scale_shape_manual(values = c("EAS" = 24, "EUR" = 21,'EAS_EUR'= 22,'Others'=20))+
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) + 
    labs(x = "Effect EAS", y = "Effect EUR", fill = "Group", shape = "Group") +  
    geom_text(aes(x = Effect_EAS + 0.01, y = Effect_EUR + 0.01, label = Gene, color = group)) +  
    geom_abline(intercept = -0.001, slope = 0.589, linetype = "dashed", color = "grey") +  
    geom_abline(intercept = 0, slope = 1, color = "black") +  
    scale_color_manual(values = c("Others" = "gray", "EAS" ="#F8766D", "EUR" = "#619CFF", "EAS_EUR" = "#00BA38")) +
    scale_fill_manual(values = c("Others" = "gray", "EAS" = "#F8766D", "EUR" = "#619CFF", "EAS_EUR" = "#00BA38"))

p
ggsave('Cross-ancestry_comparison.pdf',plot=p, width=20,height=17,unit='cm',dpi=600)



#####————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
###Figure 2B: Cross-trait LDSC
##Step1. For summary data for all phenotypes, we extracted the following columns and wrote as "sumstats.txt"
R
library(data.table)
library(tidyr)
library(dplyr)

data<-fread('~/EAS_for_meta_Age_of_smoking_initiation.txt')
data_final<-data[,c('rsID','CHR','BP','N','Effect_allele','Other_allele','EAF','BETA','SE','P')]  ####A1 is EA, A2 is OA
colnames(data_final)[1:10]<-c('SNP','CHR','BP','N','A1','A2','FRQ','BETA','SE','P')
fwrite(data_final,'~/EAS_age_of_smoking_initiation_sumstats.txt',sep='\t',col.names=T)

#####Prepare documents on age_of_smoking_initiation, smoking_initiation, CPD and Smoking_cessation in two ancestries in turn

##Step2. Prepare ldsc files according to the following code
~/ldsc/munge_sumstats.py \
--sumstats ~/EAS_LC_meta_sumstats.txt \
--N-col N \
--out ~/EAS_LC_meta_sumstats \

##Step3. Cross-trait LDSC analysis
#Give an example: EUR_LC and smoking behaviour
~/ldsc/ldsc.py \
--rg ~/EUR_LC_meta_sumstats.sumstats.gz,EUR_LUAD_meta_sumstats.sumstats.gz,EUR_LUSC_meta_sumstats.sumstats.gz,EUR_age_of_smoking_initiation_sumstats.sumstats.gz,EUR_smoking_initiation_sumstats.sumstats.gz,EUR_CPD_sumstats.sumstats.gz,EUR_smoking_cessation_sumstats.sumstats.gz \
--ref-ld-chr ~/ldsc/ldscore/eur_w_ld_chr/ \
--w-ld-chr ~/ldsc/ldscore/eur_w_ld_chr/ \
--out ~/LC_EUR_smoke_ldsc_result


##Step4. Plot
data<-fread('ldsc_result_plot.csv')
data$signif <- cut(data$P, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", ""))
data$P_ad<-  -log10(data$P)
############Map -log10P to [0,1]
min_val <- min(data$P_ad,na.rm = T)
max_val <- max(data$P_ad,na.rm = T)
data$P_nor <- (data$P_ad - min_val) / (max_val - min_val)
data$P_nor <- ifelse(is.na(data$P_nor), 1, data$P_nor)

#####Set the order of the factors
data$type <- factor(data$type, levels = c("EAS_LUSC", "EAS_LUAD", "EAS_LC","EUR_LUSC", "EUR_LUAD", "EUR_LC"))
data$pheno <- factor(data$pheno, levels = c("LC", "LUAD", "LUSC", "Age of smoking initiation", "Smoking initiation", "Cigarette per day", "Smoking cessation"))
####Convert the data to a long format suitable for ggplot2
data_melt <- melt(data, id.vars = c("type", "pheno", "signif",'P_nor'), measure.vars = "r")
x_labels <- levels(factor(data_melt$pheno))
y_labels <- levels(factor(data_melt$type))

p <- ggplot(data_melt, aes(x = pheno, y = type, color = value)) +
  geom_tile(color = "white", fill = "white") +
  geom_point(aes(size = P_nor), shape = 15) +
  geom_text(aes(label = signif), color = "white", size = 9) +
  scale_color_gradient2(low = "#6D1727", mid = "white", high = "#1E324B", midpoint = 0, limit = c(-1, 1.1), name = "Correlation") +
  scale_size_continuous(range = c(10, 52), guide = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black")
  ) +
  coord_fixed(ratio = 1) +
  labs(x = "Phenotype", y = "Type")

p <- p +
  scale_x_discrete(labels = x_labels) +
  scale_y_discrete(labels = y_labels)
for (x in 1:length(x_labels)) {
  p <- p + geom_vline(xintercept = x + 0.5, color = "gray", linetype = "solid")
}

for (y in 1:length(y_labels)) {
  p <- p + geom_hline(yintercept = y + 0.5, color = "gray", linetype = "solid")
}
p <- p + 
  geom_rect(aes(xmin = 0.5, xmax = length(x_labels) + 0.5, ymin = 0.5, ymax = length(y_labels) + 0.5),
            color = "black", fill = NA, size = 1)
print(p)
ggsave('LC_smoke_ldsc_plot.pdf',plot=p, width=40,height=30,unit='cm',dpi = 600)


#####————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
###Figure 2C: Popcorn
##Give an example: LC
#EAS_meta
EAS_ori<-fread('~/EAS_GWMA.txt')
EAS_ori<-EAS_ori[EAS_ori$Freq1 >=0.01 & EAS_ori$Freq1 <=0.99,]
EAS<-EAS_ori[,c('CHR','BP','Allele2','Allele1','Freq1','N','Effect','StdErr')]
colnames(EAS)[1:8]<-c('chr','pos','a1','a2','af','N','beta','SE')
EAS$a1<-toupper(EAS$a1)
EAS$a2<-toupper(EAS$a2)

#EUR_meta
EUR_ori<-fread('~/EUR_GWMA.txt')
EUR_ori<-EUR_ori[EUR_ori$Freq1 >=0.01 & EUR_ori$Freq1 <=0.99,]
EUR<-EUR_ori[,c('CHR','BP','Allele2','Allele1','Freq1','N','Effect','StdErr')]
colnames(EUR)[1:8]<-c('chr','pos','a1','a2','af','N','beta','SE')
EUR$a1<-toupper(EUR$a1)
EUR$a2<-toupper(EUR$a2)

##########Popcorn file
socre<-fread('~/Popcorn/LD_score/EUR_EAS_all_gen_eff.cscore')
colnames(socre)[1:3]<-c('chr','bp','rsid')
socre<-socre[,c(1:3)]
socre<-socre %>% mutate(chrbp=paste(chr,bp,sep=':'))
EUR<-EUR %>% mutate(chrbp=paste(chr,pos,sep=':'))
EAS<-EAS %>% mutate(chrbp=paste(chr,pos,sep=':'))
EUR<-left_join(socre,EUR,by ='chrbp')
EAS<-left_join(socre,EAS,by ='chrbp')
EUR<-EUR[!is.na(EUR$pos),]
EAS<-EAS[!is.na(EAS$pos),]
EUR<-EUR[,c('rsid','a1','a2','af','N','beta','SE')]
EAS<-EAS[,c('rsid','a1','a2','af','N','beta','SE')]
fwrite(EUR,'~/EUR_meta_popcorn.txt',sep='\t',col.names=T)
fwrite(EAS,'~/EAS_meta_popcorn.txt',sep='\t',col.names=T)

######effect
screen -r popcorn_eff
python -m popcorn fit -v 1 --use_mle --cfile ~/Popcorn/LD_score/EUR_EAS_all_gen_eff.cscore --gen_effect --sfile1 ~/popcorn/EUR_meta_popcorn.txt --sfile2 ~/popcorn/EAS_meta_popcorn.txt LC_eff_cor
######impact
screen -r popcorn_imp
python -m popcorn fit -v 1 --use_mle --cfile ~/Popcorn/LD_score/EUR_EAS_all_gen_imp.cscore --sfile1 ~/popcorn/EUR_meta_popcorn.txt --sfile2 ~/popcorn/EAS_meta_popcorn.txt LC_imp_cor 

######Plot
data<-fread('popcorn_result_plot.csv')
data$lower <- data$Value - 1.96 * data$SE
data$upper <- data$Value + 1.96 * data$SE

p <- ggplot(data, aes(x = Phenotype, y = Value, group = Estimator, color = Estimator, fill = Estimator,shape = Estimator,), size = 3) +
  geom_point(position = position_dodge(width = 0.5), size = 6) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_text(aes(label = ifelse(P < 0.05, "*", "")), position = position_dodge(width = 0.5), vjust = -0.5,size=7) +
  labs(x = "Phenotype", y = "Cross-ancestry genetic correlation", shape = "Estimator") +
  scale_shape_manual(values = c("Genetic effect" = 22, "Genetic impact" = 21)) +
  scale_color_manual(values = c("Genetic effect" = "#6EB8D2", "Genetic impact" = "#4F7ABB")) +
  scale_fill_manual(values = c("Genetic effect" = "#6EB8D2", "Genetic impact" = "#4F7ABB")) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2), limits = c(-0.1, 1.2)) + 
  theme_minimal() +
  theme(axis.text.x = element_text(size = 16, color = "black"),  
        axis.text.y = element_text(size = 16, color = "black"),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(),  
        axis.ticks = element_line(color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(c(0,0,20,20)) )  
p
ggsave('popcorn_plot.pdf',plot=p, width=20,height=22,unit='cm',dpi = 600)

