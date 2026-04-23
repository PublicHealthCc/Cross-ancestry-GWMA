#########Three ways to calculate PRS
####1. GWAS-reported: Get PRS19 and PRS128 information from the PGS catalog
####Give an example for calculating PRS128 in the UK Biobank
cd ~/PRS_UKB/PRS128
for chr in {1..22}  
do
plink2 --pfile ~/UKB_genotype/ukb_imp_chr${chr} \
--extract PRS128.valid.snp \
--make-bed --out ~/PRS_UKB/PRS128/UKB_PRS128_chr${chr}
done

##combine bfile files 
plink --bfile UKB_PRS128_chr1 --merge-list mergelist.txt --make-bed --out UKB_PRS128_chrall
## calculate PRS
plink  --bfile UKB_PRS128_chrall  --score PRS128.txt 1 6 9 header  --extract PRS128.valid.snp  --out PRS128_score


#####——————————————————————————————————————————————————————————————————————————————————————————————————————————
####2. C+T
###Give an example: UKB_5*10E-8
R
library(data.table)
data<-fread('Cross-ancestry_GWMA.csv')
data1<-subset(data,data$P<5E-8)
write.table(data1$V2,file='5_08_SNP.txt',quote=F,row.names=F)
data3<-data1[,c(21,1:19)]
write.table(data3,file='5_08_SNP_SUM.txt',quote=F,row.names=F,sep="\t")

~/plink1.9 --bfile ~/1000Genome/chr_all_1kgv3_2015_EAS_EUR \
--extract 5_08_SNP.txt \
--make-bed \
--out 5_08_1KG

~/plink1.9 --bfile 5_08_1KG \
--clump 5_08_SNP_P.txt \
--clump-snp-field SNP \
--clump-field P \
--clump-p1 0.0001 \
--clump-r2 0.1 \
--clump-kb 10000 \
--out 5_08_indepedent_SNP

##extract id
awk 'NR!=1{print $3}' 5_08_indepedent_SNP.clumped > 5_08.valid.snp
##Extract variants < 5*10-8 in UKB
for chr in {1..22}
do
plink --bfile ~/chr${chr}/ukb_chr${chr} \
--extract ~/5_08.valid.snp \
--make-bed --out ~/UKB_chr${chr}_5_08
done

#combine bfile file
plink --bfile UKB_chr1_5_08 --merge-list mergelist_08.txt --make-bed --out UKB_chrall_5_08
##calculate PRS_5*10-8
cd /data/ZJ/P_value/
plink  --bfile UKB_chrall_5_08  --score 5_08_SNP_SUM.txt 1 21 10 header  --extract 5_08.valid.snp  --out 5_08_snp_score

#####——————————————————————————————————————————————————————————————————————————————————————————————————————————
####3. PRS-CSx
####Give an example in UKB
## Reference panel
wget https://personal.broadinstitute.org/hhuang//public//PRS-CSx/Reference/1KG/ldblk_1kg_eas.tar.gz
wget https://personal.broadinstitute.org/hhuang//public//PRS-CSx/Reference/1KG/ldblk_1kg_eur.tar.gz
wget https://personal.broadinstitute.org/hhuang//public//PRS-CSx/Reference/1KG/snpinfo_mult_1kg_hm3
wget https://personal.broadinstitute.org/hhuang//public//PRS-CSx/Reference/UKBB/ldblk_ukbb_eas.tar.gz
wget https://personal.broadinstitute.org/hhuang//public//PRS-CSx/Reference/UKBB/ldblk_ukbb_eur.tar.gz
wget https://personal.broadinstitute.org/hhuang//public//PRS-CSx/Reference/UKBB/snpinfo_mult_ukbb_hm3
tar -zxvf ldblk_ukbb_eas.tar.gz
tar -zxvf ldblk_ukbb_eur.tar.gz
tar -zxvf ldblk_1kg_eas.tar.gz
tar -zxvf ldblk_1kg_eur.tar.gz

##Analyze
##Give an example: CHR1
## CHR1 （PRS-CSx auto）
python ~/PRScsx-master/PRScsx.py \
--ref_dir=~/ref \
--bim_prefix=~/chr1/ukb_chr1 \
--sst_file=~/EUR_GWAS_sumstats.txt,~/Asain_GWAS_sumstats.txt \
--n_gwas=291871,246963 \
--pop=EUR,EAS \
--chrom=1 \
--meta=META_FLAG  \
--seed=2 \
--out_dir=~/out \
--out_name=UKB_chr1_auto

## Calcualte score
plink  --bfile UKB_chrall_EUR_auto \
--score ~/UKB_chrall_EUR_auto.txt 2 4 6 header \
--extract ~/UKB_chrall_EUR_auto_snp.txt  --out UKB_chrall_EUR_auto_score
## 1076069 variantes

###_____________________________________________________________________________________________________________
###Figure 3A: give an example about PRS_CSx_UKB
data3<-fread('~/UKB_PRS_CSx_model.txt')
data3$SCORE_scale<-data3$SCORE*1076069
data3$SCORE_NORM<-scale(data3$SCORE_scale,center=T,scale=T)
data3[data3$Smoking_status==-3|is.na(data3$Smoking_status),]$Smoking_status <- -9
table(data3$Smoking_status)

fit_PRScsx<-coxph(Surv(surtime,LC_total)~SCORE_NORM+Age+as.factor(Sex)+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=data3)
summary(fit_PRScsx)
summary(fit_PRScsx)$coef

###_____________________________________________________________________________________________________________
###Supplementart Figure 5 PRS-CSx distribution plot: give an example about PRS_CSx_UKB
data<- fread('UKB_PRS_CSx_model.txt')
data$colour=0
data$colour[data$LC_total==0]="#00468B99"
data$colour[data$LC_total==1]="#AD002A99"

data$SCORE_scale<-data$SCORE*1076069
data$SCORE_NORM<-scale(data$SCORE_scale,center=T,scale=T)

p2 <- ggplot(data)+geom_density(aes(x = SCORE_NORM, fill = colour),alpha=0.35)+
  xlab("Polygenic Risk Score")+ylab("Proportion")+
  scale_fill_lancet(labels = c("Non-LC","LC"))+
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        plot.title = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),legend.position=c(0.2,0.85),legend.title=element_blank(),
        legend.text = element_text(size = 10))
p2
ggsave("~/Distribution_UKB_PRS_CSx.pdf", plot = p2, width = 8, height = 6, dpi = 600)
#######rank sum test
res.ftest <- var.test(data$SCORE ~ data$LC_total, data = data)
res.ftest
res <- wilcox.test(data$SCORE ~ data$LC_total, data = data, var.equal = TRUE)
res
p<-res$p.value
p


###_____________________________________________________________________________________________________________
###Figure 3B: Correlation analysis and comparison between PRS-CSx_UKB and PRS_reported
###Give an example: UKB
R
library(tidyverse)
library(survival)
library(pROC)
library(survminer)
library(ggplot2) 
library(data.table)
###PRS-CSX_UKB
data <- fread("~/UKB_PRS_CSx_model.txt")
data$SCORE_17 <- data$SCORE*1076069
data$SCORE_NORM <- scale(data$SCORE_17,center = T,scale = T)
quantile(data$SCORE_NORM,seq(0,1,0.1))
q1=quantile(data$SCORE_NORM,seq(0,1,0.1))[[1]]
q2=quantile(data$SCORE_NORM,seq(0,1,0.1))[[2]]
q3=quantile(data$SCORE_NORM,seq(0,1,0.1))[[3]]
q4=quantile(data$SCORE_NORM,seq(0,1,0.1))[[4]]
q5=quantile(data$SCORE_NORM,seq(0,1,0.1))[[5]]
q6=quantile(data$SCORE_NORM,seq(0,1,0.1))[[6]]
q7=quantile(data$SCORE_NORM,seq(0,1,0.1))[[7]]
q8=quantile(data$SCORE_NORM,seq(0,1,0.1))[[8]]
q9=quantile(data$SCORE_NORM,seq(0,1,0.1))[[9]]
q10=quantile(data$SCORE_NORM,seq(0,1,0.1))[[10]]
q11=quantile(data$SCORE_NORM,seq(0,1,0.1))[[11]]

data$GRS_10<-1
data[data$SCORE_NORM>q2 & data$SCORE_NORM<=q3,]$GRS_10=2
data[data$SCORE_NORM>q3 & data$SCORE_NORM<=q4,]$GRS_10=3
data[data$SCORE_NORM>q4 & data$SCORE_NORM<=q5,]$GRS_10=4
data[data$SCORE_NORM>q5 & data$SCORE_NORM<=q6,]$GRS_10=5
data[data$SCORE_NORM>q6 & data$SCORE_NORM<=q7,]$GRS_10=6
data[data$SCORE_NORM>q7 & data$SCORE_NORM<=q8,]$GRS_10=7
data[data$SCORE_NORM>q8 & data$SCORE_NORM<=q9,]$GRS_10=8
data[data$SCORE_NORM>q9 & data$SCORE_NORM<=q10,]$GRS_10=9
data[data$SCORE_NORM>q10 & data$SCORE_NORM<=q11,]$GRS_10=10

##Cox model 
fit_CSx_UKB_GRS<-coxph(Surv(surtime,LC_total)~as.factor(GRS_10)+as.factor(Sex)+Age+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=data)
summary(fit_CSx_UKB_GRS)
#####result outcome
outcome=function(x){
  HR = x$coef[c(1:9),2]
  HR.confint.lower <- exp(x$coef[c(1:9),1]- x$coef[c(1:9),3]*1.96)
  HR.confint.upper <- exp(x$coef[c(1:9),1]+ x$coef[c(1:9),3]*1.96)
  p.value<-signif(x$coef[c(1:9),"Pr(>|z|)"],3)
  res<-cbind(HR,HR.confint.lower, HR.confint.upper,p.value)
  return(res)}

result1=outcome(summary(fit_CSx_UKB_GRS))
result2 <- c('Reference',"","","") %>% rbind(.,result1)

colnames(result2)=c("hr","hr_l","hr_h","p")
fwrite(result2,file='PRS_CSx_UKB_GRS10_correlation_figure.xlsx',sep = "\t",row.names=T)
######plot#####
resu3<-result2
resu3 <- as.data.frame(resu3)
resu3$hr=as.numeric(resu3$hr)
resu3$hr_l=as.numeric(resu3$hr_l)
resu3$hr_h=as.numeric(resu3$hr_h)
resu3$x=c(1,2,3,4,5,6,7,8,9,10)

reg <- lm(hr~x,data=resu3)
reg #intercept 1.087 ，x= 0.156 
#when x=1, 1.087 +0.156 = 1.243
#when x=10, 1.087 + 0.156*10= 2.647

resu3[1,1] <- 1
p1=ggplot(resu3,aes(x=x,y=hr))+
  geom_linerange(aes(x=x,ymin=hr_l,ymax=hr_h),colour="#ABD9E9",size=0.8)+
  geom_point(shape=15,size=3,colour="#00468BFF")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10),labels=c('Q1','Q2','Q3','Q4','Q5','Q6','Q7',"Q8",'Q9','Q10'))+
  labs(x="Quantile of GRS",y="Hazards Ratio (95%CI)")+
  geom_hline(yintercept=1,color="#AD002A99",linetype=2,size=0.8)+
  geom_segment(aes(x =1, y =1.243, xend=10, yend = 2.647),color="#AD002A99",size=0.8)+
  theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.text=element_text(size=8,colour = "black"),
        axis.title= element_text(size =10,colour = "black",face = "bold" ))

p1
ggsave("PRS_UKB_CSx_GRS10.pdf", plot = p1, width = 8, height = 6, dpi = 600)

###PRS128_UKB
data3 <- fread("~/UKB_PRS128_model.txt")
data3$SCORE_108 <- data3$SCORE*108
data3$SCORE_NORM <- scale(data3$SCORE_108,center = T,scale = T)
quantile(data3$SCORE_NORM,seq(0,1,0.1))
q1=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[1]]
q2=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[2]]
q3=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[3]]
q4=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[4]]
q5=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[5]]
q6=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[6]]
q7=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[7]]
q8=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[8]]
q9=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[9]]
q10=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[10]]
q11=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[11]]

data3$GRS_10<-1
data3[data3$SCORE_NORM>q2 & data3$SCORE_NORM<=q3,]$GRS_10=2
data3[data3$SCORE_NORM>q3 & data3$SCORE_NORM<=q4,]$GRS_10=3
data3[data3$SCORE_NORM>q4 & data3$SCORE_NORM<=q5,]$GRS_10=4
data3[data3$SCORE_NORM>q5 & data3$SCORE_NORM<=q6,]$GRS_10=5
data3[data3$SCORE_NORM>q6 & data3$SCORE_NORM<=q7,]$GRS_10=6
data3[data3$SCORE_NORM>q7 & data3$SCORE_NORM<=q8,]$GRS_10=7
data3[data3$SCORE_NORM>q8 & data3$SCORE_NORM<=q9,]$GRS_10=8
data3[data3$SCORE_NORM>q9 & data3$SCORE_NORM<=q10,]$GRS_10=9
data3[data3$SCORE_NORM>q10 & data3$SCORE_NORM<=q11,]$GRS_10=10

##Cox model 
fit_PRS128_UKB_GRS<-coxph(Surv(surtime,LC_total)~as.factor(GRS_10)+as.factor(Sex)+Age+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=data3)
summary(fit_PRS128_UKB_GRS)
#####result outcome
outcome=function(x){
  HR = x$coef[c(1:9),2]
  HR.confint.lower <- exp(x$coef[c(1:9),1]- x$coef[c(1:9),3]*1.96)
  HR.confint.upper <- exp(x$coef[c(1:9),1]+ x$coef[c(1:9),3]*1.96)
  p.value<-signif(x$coef[c(1:9),"Pr(>|z|)"],3)
  res<-cbind(HR,HR.confint.lower, HR.confint.upper,p.value)
  return(res)}

result5=outcome(summary(fit_PRS128_UKB_GRS))
result6 <- c('Reference',"","","") %>% rbind(.,result5)

colnames(result6)=c("hr","hr_l","hr_h","p")
fwrite(result6,file='PRS128_UKB_GRS10_correlation_figure.xlsx',sep = "\t",row.names=T)
######plot#####
resu2<-result6
resu2 <- as.data.frame(resu2)
resu2$hr=as.numeric(resu2$hr)
resu2$hr_l=as.numeric(resu2$hr_l)
resu2$hr_h=as.numeric(resu2$hr_h)
resu2$x=c(1,2,3,4,5,6,7,8,9,10)

reg <- lm(hr~x,data=resu2)
reg #intercept 1.0746，x= 0.1149 
#when x=1, 1.0746+0.1149 = 1.1895
#when x=10, 1.0746+ 0.1149*10= 2.2236

resu2[1,1] <- 1
p2=ggplot(resu2,aes(x=x,y=hr))+
  geom_linerange(aes(x=x,ymin=hr_l,ymax=hr_h),colour="#ABD9E9",size=0.8)+
  geom_point(shape=15,size=3,colour="#00468BFF")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10),labels=c('Q1','Q2','Q3','Q4','Q5','Q6','Q7',"Q8",'Q9','Q10'))+
  labs(x="Quantile of GRS",y="Hazards Ratio (95%CI)")+
  geom_hline(yintercept=1,color="#AD002A99",linetype=2,size=0.8)+
  geom_segment(aes(x =1, y =1.1895, xend=10, yend = 2.2236),color="#AD002A99",size=0.8)+
  theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),axis.text=element_text(size=8,colour = "black"),
        axis.title= element_text(size =10,colour = "black",face = "bold" ))

p2
ggsave("PRS128_UKB_GRS10.pdf", plot = p2, width = 8, height = 6, dpi = 600)

#################Likelihood-ratio test##############
lrtest <- anova(fit_CSx_UKB_GRS, fit_PRS128_UKB_GRS)
lrtest

##combine plot 
resu3$set<-"A"
resu2$set<-"B"
combined<-rbind(resu3,resu2)

color_mapping <- c(A = "#AD002A99", B = "gray60")
p5 <- ggplot(combined, aes(x = x, y = hr, group = set, color = set)) +
  geom_linerange(aes(x = x, ymin = hr_l, ymax = hr_h, group = set, color = set), size = 0.8, position = position_dodge(width = 0.4)) +
  geom_point(aes(x = x, y = hr, group = set, color = set), shape = 15, size = 3, position = position_dodge(width = 0.4)) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), labels = c('Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7', 'Q8', 'Q9', 'Q10')) +
  labs(x = "Quantile of PRS", y = "Hazards Ratio (95%CI)") +
  geom_hline(yintercept = 1, linetype = 2, size = 0.8, color="#AD002A99") +
  geom_segment(data = combined[combined$set == "A", ], 
               aes(x = 1, y = 1.243, xend = 10, yend = 2.647), size = 0.8, color = "#AD002A99", position = position_dodge(width = 0.4)) +
  geom_segment(data = combined[combined$set == "B", ], 
               aes(x = 1, y = 1.1895, xend = 10, yend = 2.2236), size = 0.8, color = "gray60", position = position_dodge(width = 0.4)) +
  scale_color_manual(values = color_mapping) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 10, colour = "black", face = "bold")
  )
p5
ggsave("CSx_PRS128_UKB_GRS10.pdf", plot = p5, width = 8, height = 6, dpi = 600)


###_____________________________________________________________________________________________________________
##################################Figre 3D: give an example about PRS-CSX_UKB
####Give an example: PRS-CSX_UKB
######(1) Model1, adjust for age, sex and the top ten PC
data3 <- fread("~/UKB_PRS_CSx_model.txt")
data3$SCORE_17 <- data3$SCORE*1076069
data3$SCORE_NORM <- scale(data3$SCORE_17,center = T,scale = T)
quantile(data3$SCORE_NORM,seq(0,1,0.1))
q1=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[1]]
q2=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[2]]
q3=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[3]]
q4=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[4]]
q5=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[5]]
q6=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[6]]
q7=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[7]]
q8=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[8]]
q9=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[9]]
q10=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[10]]
q11=quantile(data3$SCORE_NORM,seq(0,1,0.1))[[11]]

data3$GRS_10<-1
data3[data3$SCORE_NORM>q2 & data3$SCORE_NORM<=q3,]$GRS_10=2
data3[data3$SCORE_NORM>q3 & data3$SCORE_NORM<=q4,]$GRS_10=3
data3[data3$SCORE_NORM>q4 & data3$SCORE_NORM<=q5,]$GRS_10=4
data3[data3$SCORE_NORM>q5 & data3$SCORE_NORM<=q6,]$GRS_10=5
data3[data3$SCORE_NORM>q6 & data3$SCORE_NORM<=q7,]$GRS_10=6
data3[data3$SCORE_NORM>q7 & data3$SCORE_NORM<=q8,]$GRS_10=7
data3[data3$SCORE_NORM>q8 & data3$SCORE_NORM<=q9,]$GRS_10=8
data3[data3$SCORE_NORM>q9 & data3$SCORE_NORM<=q10,]$GRS_10=9
data3[data3$SCORE_NORM>q10 & data3$SCORE_NORM<=q11,]$GRS_10=10

##bottom 1/5, intermediate 2/5-4/5 and top 1/5）
data3$GRS_5<-1
data3[data3$GRS_10>2 & data3$GRS_10<=8,]$GRS_5=2
data3[data3$GRS_10>8 & data3$GRS_10<=10,]$GRS_5=3
table(data3$GRS_5)

data3$surtime <- data3$surtime/365.25
fit_CSx_UKB_GRS_5<-coxph(Surv(surtime,LC_total)~factor(GRS_5)+Age+Sex+PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10,data=data3)
summary(fit_CSx_UKB_GRS_5)
summary(fit_CSx_UKB_GRS_5)$coef

##########logRank test
logrank <- survdiff(Surv(surtime,LC_total) ~ as.factor(GRS_5), data = data3)
logrank
logrank$pvalue

###plot
dat1=data3
dat_df <- with(dat1,
               data.frame(GRS_5 = c(1, 2, 3), 
                          Age=rep(mean(Age,na.rm = TRUE),3),
                          Sex=rep(mean(Sex,na.rm = TRUE),3),
                          PCA1=rep(mean(PCA1,na.rm = TRUE),3),
                          PCA2=rep(mean(PCA2,na.rm = TRUE),3),
                          PCA3=rep(mean(PCA3,na.rm = TRUE),3),
                          PCA4=rep(mean(PCA4,na.rm = TRUE),3),
                          PCA5=rep(mean(PCA5,na.rm = TRUE),3),
                          PCA6=rep(mean(PCA6,na.rm = TRUE),3),
                          PCA7=rep(mean(PCA7,na.rm = TRUE),3),
                          PCA8=rep(mean(PCA8,na.rm = TRUE),3),
                          PCA9=rep(mean(PCA9,na.rm = TRUE),3),
                          PCA10=rep(mean(PCA10,na.rm = TRUE),3)			  
               )
)

dat_df$GRS_5=as.factor(dat_df$GRS_5)
fit1 <- survfit(fit_CSx_UKB_GRS_5, newdata = dat_df)

p2=ggsurvplot(fit1,
              linetype=c(1,1,1),
              xlab=c("Years of follow-up"),
              xlim=c(0,11),
              ylim=c(0,0.010),
              ylab=c("Standardized lung cancer event rate"),
              legend=c(0.30,0.8),
              legend.title="Lung cancer polygenic risk score",
              legend.labs=c("Low (Reference)","Intermediate (HR):1.70 (1.54-1.88)","High (HR):2.09 (1.87-2.34)"),			
              conf.int=TRUE,
              conf.int.fill="gray",
              conf.int.style="ribbon",
              conf.int.alpha=0.2,
              risk.table.col="strata",
              palette=c("#00468B99","#42B540E5","#AD002AFF"),
              ggtheme=theme_bw()+theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(colour = "black", fill = NA)),
              fun="event",
              font.xtickslab=c(12,"bold"),
              font.ytickslab=c(12,"bold"),
              font.x=c(14,"bold"),
              font.y=c(14,"bold"),
              font.legend=c(10,"bold"),
              censor=FALSE,

              dat=dat_df)
p2
ggsave("UKB_Lung_Cancer_risk_predtion_of_PRS_CSx.pdf", plot = p2$plot, width = 8, height = 6, dpi = 600)

#####——————————————————————————————————————————————————————————————————————————————————————————————————————————
####4. PRS-CS
####Give an example in UKB
## download reference panel
wget https://personal.broadinstitute.org/hhuang//public//PRS-CSx/Reference/UKBB/ldblk_ukbb_eas.tar.gz
wget https://personal.broadinstitute.org/hhuang//public//PRS-CSx/Reference/UKBB/ldblk_ukbb_eur.tar.gz
wget https://personal.broadinstitute.org/hhuang//public//PRS-CSx/Reference/UKBB/snpinfo_mult_ukbb_hm3
tar -zxvf ldblk_ukbb_eas.tar.gz
tar -zxvf ldblk_ukbb_eur.tar.gz

### Prepare GWAS summary statistics
R
library(data.table)
library(dplyr)
library(stringr)
## EUR
data<-fread('./EUR_meta_fixed.txt')
data_final<- data[,c('rsID','Allele1','Allele2','Effect','StdErr')]
colnames(data_final)<-c('SNP','A1','A2','BETA','SE')
fwrite(data_final,'./EUR_meta_for_PRScs.txt',sep='\t',col.names=T)

### run PRScs
for chr in {1..22}; do
  python ./software/PRScs/PRScs.py \
    --ref_dir=./software/PRScs/LD_reference_panel/ldblk_ukbb_eur \
    --bim_prefix=./software/PRScs/bim_prefix/ukbb_bim/chr${chr} \
    --sst_file=./EUR_meta_for_PRScs.txt \
    --n_gwas=291871 \
    --chrom=${chr} \
    --out_dir=./results/EUR/EUR
done

## 1082169 variants

## Validate PRS-CSx in UKB
R
library(pROC)
library(survival)
data3 <-fread('./UKB_PRScs_model.txt')
data3$SCORE_scale<-data3$SCORE*1082169
data3$SCORE_NORM<-scale(data3$SCORE_scale,center=T,scale=T)
fit_PRScs<-coxph(Surv(surtime,LC_total)~SCORE_NORM,data=data3)
summary(fit_PRScs)
summary(fit_PRScs)$coef
logit_3 <- glm(LC_total ~ SCORE,family = binomial,data=data3)
summary(logit_3)


