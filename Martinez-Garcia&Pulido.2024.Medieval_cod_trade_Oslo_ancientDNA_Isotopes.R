#=============================================================================================================#
#
# Script created by Angélica Pulido and Lourders Martinez-Garcia
# 
# Study published in :
# Cite as :  Martinez-Garcia and Pulido et al., (2024) "Tracing 600 years of long-distance Atlantic cod trade in medieval and post-medieval Oslo using stable isotopes and ancient DNA"
#
#
# This script: run all analysis and plots in Martinez-Garcia and Pulido et al., (2024): 
#       Produce statistical analysis, figures and supplementary material
#	
# Usage notes: run line by line
#=============================================================================================================#

#===========================================================================================================================================================================####

rm(list = ls()) #clears environment

working.directory<-setwd("/path/to/your/working/directory/files.and.scripts")
working.directory
source('summarySE.R', chdir = F)

#------ Loading required libraries -----#####

# install.packages("dplyr")
# install.packages("factoextra")
# install.packages("ggpubr")

library(dplyr)
library(ggplot2)
library(reshape2)
library(forcats)
library(tidyverse)
library(factoextra)
library(ggpubr)


#--------------------------------------#####
#-------------- Paleomix --------------#####
#--------------------------------------#####

paleo<-read.table("Paleomix_summary.txt", header = T)

# Clonality vs # reads
ggplot(paleo,aes(Pair,Clonality, color=Locality)) + geom_point(size=3)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

# Endogenous DNA fraction
box_frac<-ggplot(paleo, aes(x=Locality, y=Frac,color=Locality)) + 
  # scale_fill_manual(values=c("gray", "paleturquoise3"))+
  # scale_color_manual(values=c("gray", "paleturquoise3"))+
  geom_boxplot()  + 
  ggtitle("Endogenous DNA fraction")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

box_frac + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

paleo %>% group_by(Locality) %>%
  summarise(Count = n())

detach("package:dplyr", unload = TRUE)

res<-summarySE(paleo, measurevar="Frac", groupvars="Locality")
res

#--------------------------------------#####
#---- BAM Scorer assigments -----#####
#--------------------------------------#####
library(dplyr)

LG01<-read.table("InversionCaller_LG01.txt", header = T)
LG02<-read.table("InversionCaller_LG02.txt", header = T)
LG07<-read.table("InversionCaller_LG07.txt", header = T)
LG12<-read.table("InversionCaller_LG12.txt", header = T)

head(LG01)
dim(LG01)
dim(LG02)
dim(LG07)
dim(LG12)
dim(WGS)

# removing individuals with low SNPs supporting the haplotype asssigment ####
LG01_fil<-filter(LG01,SNPs>=5)
dim(LG01_fil) # # 53 ID
LG02_fil<-filter(LG02,SNPs>=5)
dim(LG02_fil) # 50 ID
LG07_fil<-filter(LG07,SNPs>=5)
dim(LG07_fil) # 51 ID
LG12_fil<-filter(LG12,SNPs>=5)
dim(LG12_fil) # 49 ID

# View(LG01)

# LG01 ####
LG01_AA<-LG01_fil %>% 
  dplyr::group_by(Locality,AA) %>%
  dplyr::summarise(length(AA))

LG01_BB<-LG01_fil %>% 
  dplyr::group_by(Locality,BB) %>%
  dplyr::summarise(length(BB))

LG01_AB<-LG01_fil %>% 
  dplyr::group_by(Locality,AB) %>%
  dplyr::summarise(length(AB))

LG01_counts<-rbind(LG01_AA,LG01_BB,LG01_AB)

# LG02 ####
LG02_AA<-LG02_fil %>% 
  dplyr::group_by(Locality,AA) %>%
  dplyr::summarise(length(AA))

LG02_BB<-LG02_fil %>% 
  dplyr::group_by(Locality,BB) %>%
  dplyr::summarise(length(BB))

LG02_AB<-LG02_fil %>% 
  dplyr::group_by(Locality,AB) %>%
  dplyr::summarise(length(AB))

LG02_counts<-rbind(LG02_AA,LG02_BB,LG02_AB)

# LG07 ####
LG07_AA<-LG07_fil %>% 
  dplyr::group_by(Locality,AA) %>%
  dplyr::summarise(length(AA))

LG07_BB<-LG07_fil %>% 
  dplyr::group_by(Locality,BB) %>%
  dplyr::summarise(length(BB))

LG07_AB<-LG07_fil %>% 
  dplyr::group_by(Locality,AB) %>%
  dplyr::summarise(length(AB))

LG07_counts<-rbind(LG07_AA,LG07_BB,LG07_AB)


# LG12 ####
LG12_AA<-LG12_fil %>% 
  dplyr::group_by(Locality,AA) %>%
  dplyr::summarise(length(AA))

LG12_BB<-LG12_fil %>% 
  dplyr::group_by(Locality,BB) %>%
  dplyr::summarise(length(BB))

LG12_AB<-LG12_fil %>% 
  dplyr::group_by(Locality,AB) %>%
  dplyr::summarise(length(AB))

LG12_counts<-rbind(LG12_AA,LG12_BB,LG12_AB)
head(LG12_counts)

# write.table(LG01_counts,file = "LG01_counts.txt")
# write.table(LG02_counts,file = "LG02_counts.txt")
# write.table(LG07_counts,file = "LG07_counts.txt")
# write.table(LG12_counts,file = "LG12_counts.txt")

## read modified counts ####
# treshold 0.8 of assigment probability!
LG01_c<-read.table( "LG01_counts_m.txt", header = T)
LG02_c<-read.table( "LG02_counts_m.txt", header = T)
LG07_c<-read.table( "LG07_counts_m.txt", header = T)
LG12_c<-read.table( "LG12_counts_m.txt", header = T)

# head(LG01_c)
# head(LG02_c)
# head(LG07_c)
# head(LG12_c)

#####################################
#Probability of a genotypic affinity#
#####################################
library(tidyverse)
prob<-read.table("Probability_ancient.genotypic.affinity_FoodImpact2.txt", header = T)

head(prob)
levels(prob$Locality)

oslo <- filter(prob, Locality =='Oslo_Brann_1'|Locality == 'Oslo_Brann_2'|Locality =='Oslo_Brann_3'|
                 Locality =='Oslo_Brann_4' |Locality == 'Oslo_Brann_5' |
                 Locality =='Oslo_Brann_6' |Locality =='Oslo_Brann_7' |Locality == 'Oslo_Brann_8'|
                 Locality =='Oslo_Mindets_tomt')
oslo<-na.omit(oslo)

a<-ggplot(data=oslo, aes(x=probability, y=Sample, fill=type)) +
  geom_bar(stat="identity") + ggtitle("Oslo") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
a

######################################################################################################
#Independently testing Size, Layer, bone structure in determining the NEAC or Iceland / LG01 Genotype#
######################################################################################################
#### for genotype these individuals were removed because of low endogenous DNA and/or SNPs supporting population assigments: 
### aDNA_ID (Old_ID)
### COD348	(2535)
### COD379	(2570)
### COD396	(2588)
### COD406	(2598)
### COD407	(2599)
### COD410	(2602)
### COD412	(2604)
### COD415	(2607)
### COD416	(2608)
### COD417	(2609)
### COD418	(2610)
### COD420	(2612)
### COD422	(2614)
### COD423	(2615)
### COD427	(2619)

library(dplyr)
library(ggplot2)
library(reshape2)
library(forcats)
library(tidyverse)
library(factoextra)
library(dplyr)
library(ggpubr)
library(car)

genotype<- read.table("/Users/loma0114/Documents/Cod/AngelicasPaper/Analysis/NewDataSets/Files_Oslo_InUse/Data_Analysis_Angelica/All_genotypes_LMG.txt", header = T) 
genotype<- read.table("All_genotypes_LMG.txt", header = T) 
dim(genotype)# ALL 50 individuals sequenced

genotype_table <- subset(genotype,aDNA_ID!= "COD348" & aDNA_ID!= "COD379" & aDNA_ID!= "COD396" &
                           aDNA_ID!= "COD406" & aDNA_ID!= "COD407" & aDNA_ID!= "COD410" & aDNA_ID!= "COD412" &
                           aDNA_ID!= "COD415" & aDNA_ID!= "COD416" & aDNA_ID!= "COD417" &
                           aDNA_ID!= "COD418" & aDNA_ID!= "COD420" & aDNA_ID!= "COD422" &
                           aDNA_ID!= "COD423" & aDNA_ID!= "COD427", 
                         select = c(aDNA_ID, Old_ID, Element, Layer, Layer2, Age, Size, 
                                    LG01, LG02, LG07, LG12, 
                                    LofotenSkrei, LofotenCoastal, NorthSea, Oresund,
                                    IrishSea, Iceland, IcFrontal, IcCoastal, WestCoast,
                                    northcentral, northernmost, population.assigment, ncno))
head(genotype_table)
dim(genotype_table) # 35 individuals
genotype_composite<-na.omit(genotype_table) #COD408 and COD21 are removed because they have a undistinguisable northcentral or northermost origin.
dim(genotype_composite) # 33 individuals
# View(genotype_table)

#Transform the genotypic probabilities to  logit() -> logit(p) = log[p/(1 - p)] for the proportion p.
library(car)
#Logit fuction
genotype_table$LofotenSkrei_log <- logit(genotype_table$LofotenSkrei) #In logit(probNEA) : proportions remapped to (0.025, 0.975)
genotype_table$Iceland_log <- logit(genotype_table$Iceland)
genotype_table$IcFrontal_log <- logit(genotype_table$IcFrontal)
genotype_composite$LofotenSkrei_log <- logit(genotype_composite$LofotenSkrei) #In logit(probNEA) : proportions remapped to (0.025, 0.975)
genotype_composite$Iceland_log <- logit(genotype_composite$Iceland)
genotype_composite$IcFrontal_log <- logit(genotype_composite$IcFrontal)
#Checking new distribution
hist(genotype_table$LofotenSkrei_log)
shapiro.test(genotype_table$LofotenSkrei_log) #W = 0.73974, p-value = 1.638e-06* significantly deviate from a normal distribution
qqnorm(genotype_table$LofotenSkrei_log)
qqline(genotype_table$LofotenSkrei_log)
hist(genotype_composite$LofotenSkrei_log)
shapiro.test(genotype_composite$LofotenSkrei_log) #W = 0.74765, p-value = 3.747e-06* significantly deviate from a normal distribution
qqnorm(genotype_composite$LofotenSkrei_log)
qqline(genotype_composite$LofotenSkrei_log)

hist(genotype_table$Iceland_log)
shapiro.test(genotype_table$Iceland_log) #W = 0.91097, p-value = 0.007901* significantly deviate from a normal distribution
qqnorm(genotype_table$Iceland_log)
qqline(genotype_table$Iceland_log)
hist(genotype_composite$Iceland_log)
shapiro.test(genotype_composite$Iceland_log) #W = 0.91521, p-value = 0.01353* significantly deviate from a normal distribution
qqnorm(genotype_composite$Iceland_log)
qqline(genotype_composite$Iceland_log)

hist(genotype_table$IcFrontal_log)
shapiro.test(genotype_table$IcFrontal_log) #W = 0.86108, p-value = 0.0004157* significantly deviate from a normal distribution
qqnorm(genotype_table$IcFrontal_log)
qqline(genotype_table$IcFrontal_log)
hist(genotype_composite$IcFrontal_log)
shapiro.test(genotype_composite$IcFrontal_log) #W = 0.85665, p-value = 0.0004785* significantly deviate from a normal distribution
qqnorm(genotype_composite$IcFrontal_log)
qqline(genotype_composite$IcFrontal_log)

###############################
#Size vs Genetic probabilities#
###############################
#Size vs NEA probability
dim(genotype_composite) #33 individuals
#Testing for homoscedasticity
leveneTest(LofotenSkrei_log~Size, data = genotype_composite) #F-value = 1.4367; p-value = 0.2536 = equal variance
#Kruskal Wallis test
kruskal.test(LofotenSkrei_log ~ Size, data = genotype_composite)
#Kruskal-Wallis chi-squared = 2.005, df = 2, p-value = 0.367

level_order <- c("50-80cm", "80-100cm", ">100cm")
p<-genotype_composite %>%
  ggplot(aes(y=LofotenSkrei_log , x= factor(Size,level = level_order), colour=Size, fill =Size)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = 0.1), size = 6.5, colour="black") + 
  #geom_text(aes(label=aDNA_ID),size=3,hjust=1,vjust=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

#Size vs Iceland probability 
#Testing for homoscedasticity
leveneTest(Iceland_log~Size, data = genotype_composite) #F-value = 2.1839; p-value = 0.1302 = equal variance
#Kruskal Wallis test
kruskal.test(Iceland_log ~ Size, data = genotype_composite)
#Kruskal-Wallis chi-squared = 0.25356, df = 2, p-value = 0.8809
p<-genotype_composite %>%
  ggplot(aes(y=Iceland_log , x= factor(Size,level = level_order), colour=Size, fill = Size)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = 0.1), size = 6.5, colour="black") +
  #geom_text(aes(label=aDNA_ID),size=3,hjust=1,vjust=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

#Size vs Iceland Frontal probability 
#Testing for homoscedasticity
leveneTest(IcFrontal_log~Size, data = genotype_composite) #F-value = 0.874; p-value = 0.4276 = equal variance
#Kruskal Wallis test
kruskal.test(IcFrontal_log ~ Size, data = genotype_composite)
#Kruskal-Wallis chi-squared = 0.035138, df = 2, p-value = 0.9826
p<-genotype_composite %>%
  ggplot(aes(y=IcFrontal_log , x= factor(Size,level = level_order), colour=Size, fill=Size)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") +
  #geom_text(aes(label=aDNA_ID),size=3,hjust=1,vjust=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

########################################
#Bone element vs. Genetic probabilities#
########################################
head(genotype_composite) #Still using the data set without COD408 and COD421
# some elements can be combined and categorized as belonging to the a cranial structure
# adding a column summarizing the cranial bones vs cleithrum and vertebrate
genotype_composite$bone_element<-NA
levels(genotype_composite$Element)
genotype_mod<-genotype_composite %>%
  mutate(bone_element = case_when(Element=="Articular" ~ "Cranial",
                                  Element=="Dentary" ~ "Cranial",
                                  Element=="Maxilla" ~ "Cranial",
                                  Element=="Premaxilla" ~ "Cranial",
                                  Element=="Cleithrum" ~ "Postcranial",
                                  Element=="Vertebrae" ~ "Postcranial",
                                  Element=="Ceratohyal" ~ "Cranial",
                                  Element=="Parasphenoid" ~ "Cranial")) 
# Bone Element vs. NEA probability
head(genotype_mod)
dim(genotype_mod) #33 individuals 

#Testing for homoscedasticity
leveneTest(LofotenSkrei_log~bone_element, data = genotype_mod) #F-value = 17.369; p-value = 0.0002288 *** = unequal variance

#Kruskal Wallis test
kruskal.test(LofotenSkrei_log ~ bone_element, data = genotype_mod)
#Kruskal-Wallis chi-squared = 2.288, df = 1, p-value = 0.1304

#Plot
p<-genotype_mod %>%
  ggplot(aes(y=LofotenSkrei_log , x=bone_element, colour=bone_element, fill=bone_element)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour="black") +
  #geom_text(aes(label=aDNA_ID),size=3,hjust=1,vjust=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

# Bone Element vs. Iceland probability
dim(genotype_mod) #33 individuals

#Testing for homoscedasticity
leveneTest(Iceland_log~bone_element, data = genotype_mod) #F-value = 0.1312; p-value = 0.7197 = equal variance

#Kruskal Wallis test
kruskal.test(Iceland_log ~ bone_element, data = genotype_mod)
#Kruskal-Wallis chi-squared = 0.12667, df = 1, p-value = 0.7219

#Bone element vs. Iceland Frontal probability
#Testing for homoscedasticity
leveneTest(IcFrontal_log~bone_element, data = genotype_mod) #F-value = 0.064; p-value = 0.802 = equal variance

#Kruskal Wallis test
kruskal.test(IcFrontal_log ~ bone_element, data = genotype_mod)
#Kruskal-Wallis chi-squared = 0.7917, df = 1, p-value = 0.3736

#Plot
p<-genotype_mod %>%
  ggplot(aes(y=Iceland_log, x=bone_element, colour=bone_element, fill=bone_element)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") +
  #geom_text(aes(label=aDNA_ID),size=3,hjust=1,vjust=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

################################
#Time vs. Genetic probabilities#
################################
# Fire (brann) Layer vs NEA probability
dim(genotype_table) #35 individuals
genotype_composite3 <-subset(genotype_table,aDNA_ID!="COD408" & aDNA_ID!="COD421" &
                               aDNA_ID!="COD409" & aDNA_ID!="COD411" & aDNA_ID!="COD413" & 
                               aDNA_ID!="COD414")
#These individuals are taking away as they have an uncertain northernmost-northcentral affinity, and because their fire layers 5, 6, 7, 8 are ignored because of small sample size for this test (<2 or 1)
dim(genotype_composite3) #29 individuals
#Testing for homoscedasticity
leveneTest(LofotenSkrei_log~Layer2, data = genotype_composite3) #F-value = 1.7199; p-value = 0.1885 = equal variance (with 29 individuals)
#F-value = 2.2036; p-value = 0.1116 = equal variance (with 30 individuals)
#Kruskal Wallis test
kruskal.test(LofotenSkrei_log ~ Layer2, data = genotype_composite3)
#Kruskal-Wallis chi-squared = 0.74273, df = 3, p-value = 0.8631 (with 29 individuals)
#Kruskal-Wallis chi-squared = 0.73402, df = 3, p-value = 0.8652 (with 30 individuals)
#Plot
p<-ggplot(genotype_composite3,aes(y=LofotenSkrei_log , x=Layer2, colour=Layer2, fill=Layer2)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour="black") + 
  #scale_color_manual(values=c("blue","paleturquoise3", "yellow3"))+
  #geom_text(aes(label=aDNA_ID),size=3,hjust=1,vjust=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

# Fire (brann) Layer vs Iceland probability
dim(genotype_composite3) #29 individuals
#Testing for homoscedasticity
leveneTest(Iceland_log~Layer2, data = genotype_composite3) #F-value = 0.5911; p-value = 0.6266 = equal variance
#Kruskal Wallis test
kruskal.test(Iceland_log ~ Layer2, data = genotype_composite3)
#Kruskal-Wallis chi-squared = 0.035607, df = 3, p-value = 0.9982
#Plot 
p<-ggplot(genotype_composite3,aes(y=Iceland_log , x=Layer2, colour=Layer2, fill=Layer2)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") + 
  #scale_color_manual(values=c("blue","paleturquoise3", "yellow3"))+
  #geom_text(aes(label=aDNA_ID),size=3,hjust=1,vjust=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

# Fire (brann) Layer vs Iceland Frontal probability
dim(genotype_composite3) #29 individuals
#Testing for homoscedasticity
leveneTest(IcFrontal_log~Layer2, data = genotype_composite3) #F-value = 0.5045; p-value = 0.6827 = equal variance
#Kruskal Wallis test
kruskal.test(IcFrontal_log ~ Layer2, data = genotype_composite3)
#Kruskal-Wallis chi-squared = 0.96406, df = 3, p-value = 0.8099
#Plot 
p<-ggplot(genotype_composite3,aes(y=IcFrontal_log , x=Layer2, colour=Layer2)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") + 
  #scale_color_manual(values=c("blue","paleturquoise3", "yellow3"))+
  #geom_text(aes(label=aDNA_ID),size=3,hjust=1,vjust=1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

# Fire (brann) Layer vs genetic probability (as a binary variable = population.assigment)
dim(genotype_composite3) #29 individuals
#Plot & fisher test (small sample size)
p <- data.frame(
  "A" = c(1,1,2,4,3,1),
  "B" = c(0,1,2,2,1,0),
  "C" = c(0,0,3,2,2,0),
  "D" = c(0,3,0,1,0,0),
  row.names = c("NC-swNC", "NC-NS", "NS", "IC", "NEA", "IS"),
  stringsAsFactors = FALSE)
p

mosaicplot(p,
           main = "Mosaic plot",
           color = TRUE)

DATA <- c()
for (row in rownames(p)) {
  for (col in colnames(p)) {
    DATA <- rbind(DATA, matrix(rep(c(row, col), p[row, col]), ncol = 2, byrow = TRUE))
  }
}
DATAGEN <- as.data.frame(DATA)
colnames(DATAGEN) <- c("Population", "Time")
df
SelectedTest <- fisher.test(table(DATAGEN))
SelectedTest

install.packages("ggstatsplot")
library(ggstatsplot)
ggbarstats(
  DATAGEN, Population, Time,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(SelectedTest$p.value < 0.001, "< 0.001", round(SelectedTest$p.value, 3))
  )
)

# Fire (brann) Layer vs genetic probability (as a binary variable = ncno)
dim(genotype_composite3) #29 individuals
#Plot & fisher test (small sample size)
p <- data.frame(
  "A" = c(5,7),
  "B" = c(3,3),
  "C" = c(3,4),
  "D" = c(3,1),
  row.names = c("nc", "no"),
  stringsAsFactors = FALSE)
p

mosaicplot(p,
           main = "Mosaic plot",
           color = TRUE)

DATA <- c()
for (row in rownames(p)) {
  for (col in colnames(p)) {
    DATA <- rbind(DATA, matrix(rep(c(row, col), p[row, col]), ncol = 2, byrow = TRUE))
  }
}
DATAGEN <- as.data.frame(DATA)
colnames(DATAGEN) <- c("Population", "Time")
df
SelectedTest <- fisher.test(table(DATAGEN))
SelectedTest

library(ggstatsplot)
ggbarstats(
  DATAGEN, Population, Time,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(SelectedTest$p.value < 0.001, "< 0.001", round(SelectedTest$p.value, 3))
  )
)

# Fire (brann) Layer vs genetic probability (as a binary variable = Norway or Iceland). Only 70% or higher probability and only first 4 fire layers.
#Plot & fisher test (small sample size)
p <- data.frame(
  "A" = c(3,2),
  "B" = c(1,1),
  "C" = c(2,2),
  "D" = c(0,1),
  row.names = c("NEA", "Iceland"),
  stringsAsFactors = FALSE)
p

mosaicplot(p,
           main = "Mosaic plot",
           color = TRUE)

DATA <- c()
for (row in rownames(p)) {
  for (col in colnames(p)) {
    DATA <- rbind(DATA, matrix(rep(c(row, col), p[row, col]), ncol = 2, byrow = TRUE))
  }
}
DATAGEN <- as.data.frame(DATA)
colnames(DATAGEN) <- c("Population", "Time")
df
SelectedTest <- fisher.test(table(DATAGEN))
SelectedTest

library(ggstatsplot)
ggbarstats(
  DATAGEN, Population, Time,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(SelectedTest$p.value < 0.001, "< 0.001", round(SelectedTest$p.value, 3))
  )
)

#################################################
#### Integrating isotopes, size and genotype ####
#################################################
master_table<-read.table("MASTER_isotopes.genotype_LMG.txt", header = T) 
dim(master_table) # 93 individuals 

#Checking for Normal distributed data for Carbon
hist(master_table$d13C) #Skewed to the right side
shapiro.test(master_table$d13C) #W = 0.67512, p-value = 4.865e-13 :  significantly deviate from a normal distribution
qqnorm(master_table$d13C)
qqline(master_table$d13C) #Not fitted line
#Checking for Normal distributed data for Nitrogen
hist(master_table$d15N)
shapiro.test(master_table$d15N) #W = 0.97731, p-value = 0.1051 :  normal distribution
qqnorm(master_table$d15N)
qqline(master_table$d15N) #Normal distributed
#Checking for normal distributed data for Hydrogen
hist(master_table$Meand2HV.SMOW.corrected)
shapiro.test(master_table$Meand2HV.SMOW.corrected) #W = 0.9904, p-value = 0.8021 : normal distribution
qqnorm(master_table$Meand2HV.SMOW.corrected)
qqline(master_table$Meand2HV.SMOW.corrected) #Normal distributed
#Checking for normal distributed data for Sulphur
hist(master_table$Meand34SV.corrected)
shapiro.test(master_table$Meand34SV.corrected) #W = 0.92747, p-value = 0.0001649 :  significantly deviate from a normal distribution
qqnorm(master_table$Meand34SV.corrected)
qqline(master_table$Meand34SV.corrected) #Not fitted line

#Density plots and correlations of isotopic variables 
library(GGally)
ggpairs(master_table[,c(10:11,13,15:18)])

#### for genotype these individuals need to be removed: COD348, COD379, COD420, COD415
#### and for isotopes the individuals removed follow these thresholds:
# C.N_ratio: 2.9 - 3.6
# C.S_ratio: 175±50 = 125 - 225
# N.S_ratio: 60±20 = 40 - 80

##############
#PCA Isotopes#
##############
head(master_table)
dim(master_table) #93 individuals
pca_master_table <- subset(master_table, C.N_ratio>2.9& C.N_ratio<3.6 & C.S_ratio>125 & 
                             C.S_ratio<225 & N.S_ratio>40 & N.S_ratio<80,
                           select = c(aDNA_ID,Old_ID, Element,Layer,Layer2,Age,Size,d13C,d15N,Meand2HV.SMOW.corrected,Meand34SV.corrected,population.assigment, northcentral, northernmost,ncno, C.N_ratio, Size2)) 
dim(pca_master_table) # 64 individuals, passing all isotope quality checks
head(pca_master_table)
C.N_table <- subset(master_table, C.N_ratio>2.9& C.N_ratio<3.6,
                    select = c(aDNA_ID, Old_ID, Element,Layer,Layer2,Age,Size, Size2,d13C,d15N,Meand2HV.SMOW.corrected,Meand34SV.corrected,
                               LG01, LG02, LG07, LG12, 
                               LofotenSkrei, LofotenCoastal, NorthSea,  Oresund, IrishSea, Iceland, IcFrontal, IcCoastal, WestCoast, northcentral, northernmost,
                               population.assigment,ncno, C.N_ratio)) 
dim(C.N_table) #86 individuals - For only Carbon, Nitrogen and Hydrogen data (passed the C.N_ratio threshold)

all.thre_isotope_table <- subset(master_table, C.N_ratio>2.9& C.N_ratio<3.6 & C.S_ratio>125 & 
                                   C.S_ratio<225 & N.S_ratio>40 & N.S_ratio<80,
                                 select = c(aDNA_ID,Old_ID, Element,Layer,Layer2,Age,Size, Size2,d13C,d15N,Meand2HV.SMOW.corrected,Meand34SV.corrected,
                                            LG01, LG02, LG07, LG12, 
                                            LofotenSkrei, LofotenCoastal, NorthSea,  Oresund, IrishSea, Iceland, IcFrontal, IcCoastal, WestCoast, northcentral, northernmost,
                                            population.assigment, ncno, C.N_ratio)) 
dim(all.thre_isotope_table) #64 individuals but with all information

#Transforming Carbon and Sulphur data
#First, change carbon values to absolute and then transform to logarithm to achieve a normal distribution
master_table$absd13C <- abs(master_table$d13C)
master_table$logd13C <- log(master_table$absd13C)
pca_master_table$absd13C <- abs(pca_master_table$d13C)
pca_master_table$logd13C <- log(pca_master_table$absd13C)
C.N_table$absd13C <- abs(C.N_table$d13C)
C.N_table$logd13C <- log(C.N_table$absd13C)
#Normal distributed data for log Carbon
hist(master_table$logd13C) #Skewed to the right side
shapiro.test(master_table$logd13C) #W = 0.7859, p-value = 2.601e-10 :  significantly deviate from a normal distribution
qqnorm(master_table$logd13C)
qqline(master_table$logd13C)
hist(pca_master_table$logd13C) #Skewed to the right side
shapiro.test(pca_master_table$logd13C) #W = 0.95501, p-value = 0.02033 :  significantly deviate from a normal distribution
qqnorm(pca_master_table$logd13C)
qqline(pca_master_table$logd13C) 
hist(C.N_table$logd13C) #Skewed to the right side
shapiro.test(C.N_table$logd13C) #W = 0.97875, p-value = 0.1681 :  normal distribution
qqnorm(C.N_table$logd13C)
qqline(C.N_table$logd13C)

#Sulphur
master_table$logd34S <- log(master_table$Meand34SV.corrected)
pca_master_table$logd34S <- log(pca_master_table$Meand34SV.corrected)
all.thre_isotope_table$logd34S <- log(all.thre_isotope_table$Meand34SV.corrected)
#Normal distributed data for log Sulphur
hist(master_table$logd34S) #Skewed to the right side
shapiro.test(master_table$logd34S) #W = 0.85226, p-value = 1.284e-07 :  significantly deviate from a normal distribution
qqnorm(master_table$logd34S)
qqline(master_table$logd34S)
hist(pca_master_table$logd34S) #Skewed to the right side
shapiro.test(pca_master_table$logd34S) #W = 0.81449, p-value = 1.595e-07 :  significantly deviate from a normal distribution
qqnorm(pca_master_table$logd34S)
qqline(pca_master_table$logd34S) 
hist(all.thre_isotope_table$logd34S) #Skewed to the right side
shapiro.test(all.thre_isotope_table$logd34S) #W = 0.81449, p-value = 1.595e-07 :  significantly deviate from a normal distribution
qqnorm(all.thre_isotope_table$logd34S)
qqline(all.thre_isotope_table$logd34S) 

# PCA - Size - Population 
iso_pca<-princomp(pca_master_table[,c(8:11)], cor = FALSE, scores = TRUE) #Original values for all isotopes
iso_pca$loadings
#iso_pca<-princomp(pca_master_table[,c(9,10,16,17)], cor = FALSE, scores = TRUE) #Note I have values for log C and S
#iso_pca$loadings
#Loadings:
#                         Comp.1 Comp.2 Comp.3 Comp.4
#d13C                                   0.643  0.764
#d15N                            0.129  0.754 -0.643
#Meand2HV.SMOW.corrected  0.998                     
#Meand34SV.corrected            -0.990  0.124       
#               Comp.1 Comp.2 Comp.3 Comp.4
#SS loadings      1.00   1.00   1.00   1.00
#Proportion Var   0.25   0.25   0.25   0.25
#Cumulative Var   0.25   0.50   0.75   1.00

# eigen values
a<-fviz_eig(iso_pca) # from the package "factoextra"
a$data
#  dim        eig
#1   1 95.2830033
#2   2  3.9026781
#3   3  0.6054318
#4   4  0.2088868

# merging scores with descriptive info
iso<-cbind(pca_master_table[,c(1:7,12,15)],iso_pca$scores)
head(iso)

# plot - PCs vs Body Size vs Population assignment
#PC1 vs PC2
p.iso1<-iso%>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  ggplot(aes(Comp.1,Comp.2 )) + geom_point(aes(size=Size))+
  geom_text(aes(label=population.assigment),size=2.5,hjust=0.5,vjust=1.8)+ 
  #scale_color_manual(values=c("#F4C67D","#FFA500","#FA9837", "#C4C4C4", "666666", "949494")) +
  #stat_ellipse(geom="polygon", aes(fill = Size), alpha = 0.2, show.legend = FALSE, level = 0.95) +
  xlab(paste0("PC1 (95.28%)")) +
  ylab(paste0("PC2 (3.90%)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p.iso1

# plot of PC1 vs PC2 to have mean points
new<-iso%>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm")))
groups <- as.factor(new$Size)
p <- fviz_pca_ind(iso_pca,
             geom.ind = c("point"),
             pointsize = 5,
             pointshape = 16,
             col.ind = groups, #color by groups
             #col = "Set1",
             #legend.title = "Groups",
             mean.point.size = 9,
             mean.point.shape = 17)
             #ggtheme = theme(legend.text = element_text(size = 7))) 
             #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             #scale_shape_manual(values=c(19,20,21)),
             #addEllipses = TRUE, # Concentration ellipses
             #ellipse.type = "confidence")
p
p + theme_classic() + geom_point(aes(#shape = factor(new$Size),
                   color = factor(new$population.assigment),
                   size = new$Size))

#PC2 vs PC3
p.iso23<-iso%>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  ggplot(aes(Comp.3,Comp.2)) + geom_point(aes(size=Size, colour = ncno))+
  #geom_text(aes(label=Old_ID),size=5,hjust=0.5,vjust=1.8)+ 
  #scale_color_manual(values=c("#F4C67D","#FFA500","#FA9837", "#C4C4C4", "666666", "949494")) +
  #stat_ellipse(geom="polygon", aes(fill = Size), alpha = 0.2, show.legend = FALSE, level = 0.95) +
  xlab(paste0("PC3 (0.60%)")) +
  ylab(paste0("PC2 (3.90%)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p.iso23

#PC2 vs PC4
p.iso24<-iso%>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  ggplot(aes(Comp.2,Comp.4)) + geom_point(aes(size=Size, colour= ncno))+
  #geom_text(aes(label=Old_ID),size=5,hjust=0.5,vjust=1.8)+ 
  #scale_color_manual(values=c("#F4C67D","#FFA500","#FA9837", "#C4C4C4", "666666", "949494")) +
  #stat_ellipse(geom="polygon", aes(fill = Size), alpha = 0.2, show.legend = FALSE, level = 0.95) +
  xlab(paste0("PC2 (3.90%)")) +
  ylab(paste0("PC4 (0.21%)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p.iso24

#PC3 vs PC4
p.iso34<-iso%>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  ggplot(aes(Comp.3,Comp.4 )) + geom_point(aes(size=Size, colour = ncno))+
  geom_text(aes(label=Old_ID),size=5,hjust=0.5,vjust=1.8)+ 
  #scale_color_manual(values=c("#F4C67D","#FFA500","#FA9837", "#C4C4C4", "666666", "949494")) +
  #stat_ellipse(geom="polygon", aes(fill = Size), alpha = 0.2, show.legend = FALSE, level = 0.95) +
  xlab(paste0("PC3 (0.6054%)")) +
  ylab(paste0("PC4 (0.2089%)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p.iso34

# plot of PC3 vs PC4 to have mean points
new34<-iso%>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm")))
groups <- as.factor(new34$Size)
p <- fviz_pca_ind(iso_pca, axes = c(3, 4),
                  geom.ind = c("point"),
                  pointsize = 5,
                  pointshape = 16,
                  col.ind = groups, #color by groups
                  #col = "Set1",
                  #legend.title = "Groups",
                  mean.point.size = 9,
                  mean.point.shape = 17)
#ggtheme = theme(legend.text = element_text(size = 7))) 
#palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#scale_shape_manual(values=c(19,20,21)),
#addEllipses = TRUE, # Concentration ellipses
#ellipse.type = "confidence")
p
p + theme_classic() + geom_point(aes(#shape = factor(new$Size),
  color = factor(new34$ncno),
  size = new34$Size))

#Contribution of each variable to the PCs
variants<-fviz_pca_var(iso_pca)
variants$data 
#                                          name          x           y       coord        cos2    contrib
#d13C                                       d13C  0.4404281  0.08045176   0.2004494   0.2004494  0.1500934
#d15N                                       d15N  0.4710297  0.29664195   0.3098655   0.3098655  0.2320225
#Meand2HV.SMOW.corrected Meand2HV.SMOW.corrected 11.3000814 -0.10239578 127.7023249 127.7023249 95.6215633
#Meand34SV.corrected         Meand34SV.corrected -0.4326919 -2.26932865   5.3370748   5.3370748  3.9963207

#PC1
variants2<-fviz_contrib(iso_pca, choice = "var", axes = 1)
variants2$data
#                                           name    contrib
#d13C                                       d13C  0.1511961
#d15N                                       d15N  0.1729367
#Meand2HV.SMOW.corrected Meand2HV.SMOW.corrected 99.5299360
#Meand34SV.corrected         Meand34SV.corrected  0.1459312
#PC2
variants3<-fviz_contrib(iso_pca, choice = "var", axes = 2)
variants3$data
#name    contrib
#d13C                                       d13C  0.1231727
#d15N                                       d15N  1.6745897
#Meand2HV.SMOW.corrected Meand2HV.SMOW.corrected  0.1995296
#Meand34SV.corrected         Meand34SV.corrected 98.0027080
#PC3
variants4<-fviz_contrib(iso_pca, choice = "var", axes = 3)
variants4$data
#                                           name    contrib
#d13C                                       d13C 41.3631508
#d15N                                       d15N 56.8391390
#Meand2HV.SMOW.corrected Meand2HV.SMOW.corrected  0.2678776
#Meand34SV.corrected         Meand34SV.corrected  1.5298326
#PC4
variants5<-fviz_contrib(iso_pca, choice = "var", axes = 4)
variants5$data
#                                           name      contrib
#d13C                                       d13C 58.362480452
#d15N                                       d15N 41.313334658
#Meand2HV.SMOW.corrected Meand2HV.SMOW.corrected  0.002656787
#Meand34SV.corrected         Meand34SV.corrected  0.321528103

#Anova in each PC
#PC vs. Size
a<-aov(data = iso, formula = Comp.1~Size)
summary(a) # F-value = 8.85 ; Pr(>F) = 0.000422 ***
TukeyHSD(a)
#                       diff        lwr       upr     p adj
#50-80cm->100cm   -14.306380 -22.737046 -5.875715 0.0003908 * (small vs. large) 
#80-100cm->100cm   -4.911067 -12.506927  2.684793 0.2737513
#80-100cm-50-80cm   9.395314   2.079841 16.710786 0.0084684 * (small vs. medium)

b<-aov(data = iso, formula = Comp.2~Size)
summary(b) # F-value = 0.122 ; Pr(>F) = 0.886
TukeyHSD(b)
#                        diff       lwr      upr     p adj
#50-80cm->100cm   -0.31120390 -2.245374 1.622966 0.9210635
#80-100cm->100cm  -0.34479318 -2.087441 1.397855 0.8831740
#80-100cm-50-80cm -0.03358928 -1.711910 1.644732 0.9987265

c<-aov(data = iso, formula = Comp.3~Size)
summary(c) # F-value = 2.252 ; Pr(>F) = 0.114
TukeyHSD(c)
#                         diff        lwr       upr     p adj
#50-80cm->100cm   -0.551197667 -1.2878134 0.1854180 0.1789607
#80-100cm->100cm  -0.544043991 -1.2077199 0.1196319 0.1285699
#80-100cm-50-80cm  0.007153675 -0.6320238 0.6463312 0.9996016

d<-aov(data = iso, formula = Comp.4~Size)
summary(d) # F-value = 0.493 ; Pr(>F) = 0.613
TukeyHSD(d)
#                        diff        lwr       upr     p adj
#50-80cm->100cm   -0.07539305 -0.5201812 0.3693951 0.9127984
#80-100cm->100cm   0.08230427 -0.3184409 0.4830495 0.8747363
#80-100cm-50-80cm  0.15769733 -0.2282551 0.5436497 0.5912519

#PC vs. Time (note that all layer for time are considered in this test)
a<-aov(data = iso, formula = Comp.1~Layer2)
summary(a) # F-value = 1.53; Pr(>F) = 0.185
TukeyHSD(a)

b<-aov(data = iso, formula = Comp.2~Layer2)
summary(b) # F-value = 3.662; Pr(>F) = 0.00382 ***
TukeyHSD(b)
#diff        lwr       upr     p adj
#B-A  0.6140418 -1.9605567  3.188640 0.9901889
#C-A  1.2836699 -1.3960572  3.963397 0.7644438
#D-A  2.4221850 -0.5692282  5.413598 0.1881426
#E-A  1.3982213 -3.5874673  6.383910 0.9773563
#F-A  3.0378870 -3.6511161  9.726890 0.8058938
#G-A  8.6362291  1.9472259 15.325232 0.0039307*
#C-B  0.6696281 -1.2967540  2.636010 0.9421271
#D-B  1.8081431 -0.5655194  4.181806 0.2490055
#E-B  0.7841794 -3.8572441  5.425603 0.9985144
#F-B  2.4238452 -4.0126511  8.860341 0.9091312
#G-B  8.0221872  1.5856909 14.458683 0.0060019*
#D-C  1.1385151 -1.3487834  3.625814 0.8001620
#E-C  0.1145513 -4.5860009  4.815104 1.0000000
#F-C  1.7542171 -4.7250471  8.233481 0.9810520
#G-C  7.3525591  0.8732949 13.831823 0.0164442*
#E-D -1.0239637 -5.9089209  3.860994 0.9950916
#F-D  0.6157021 -5.9985613  7.229965 0.9999527
#G-D  6.2140441 -0.4002193 12.828307 0.0788427
#F-E  1.6396658 -6.0841298  9.363461 0.9947419
#G-E  7.2380078 -0.4857878 14.961803 0.0802226
#G-F  5.5983420 -3.3203289 14.517013 0.4770609

c<-aov(data = iso, formula = Comp.3~Layer2)
summary(c) # F-value = 0.731; Pr(>F) = 0.627
TukeyHSD(c)

d<-aov(data = iso, formula = Comp.4~Layer2)
summary(d) # F-value = 1.892; Pr(>F) = 0.0978
TukeyHSD(d)

# taking away fire layers E,F,G: 
isosubset<-subset(iso, Layer2!= "E" & Layer2!= "F" &Layer2!= "G")
a<-aov(data = isosubset, formula = Comp.1~Layer2)
summary(a) # F-value = 2.517; Pr(>F) = 0.0674
TukeyHSD(a)

b<-aov(data = isosubset, formula = Comp.2~Layer2)
summary(b) # F-value = 2.552; Pr(>F) = 0.0647
TukeyHSD(b)

c<-aov(data = isosubset, formula = Comp.3~Layer2)
summary(c) # F-value = 0.811; Pr(>F) = 0.493
TukeyHSD(c)

d<-aov(data = isosubset, formula = Comp.4~Layer2)
summary(d) # F-value = 2.178; Pr(>F) = 0.101
TukeyHSD(d)

#PCA with only NEA or Iceland individuals with >70% probability
#First removing individuals with a source populations with less than 70% probability.
#dim(master_table) #93 individuals
#dim(pca_master_table) #64 individuals
#dim(all.thre_isotope_table) #64 individuals but with all information

#all.thre.isotope_genotype_table2<-na.omit(all.thre_isotope_table) #20 individuals
#dim(all.thre.isotope_genotype_table2) #20 individuals
#Leaving just NEA and Iceland
#all.thre.isotope_genotype_table4 <-subset(all.thre.isotope_genotype_table2, 
#                                          aDNA_ID!="COD332" & aDNA_ID!="COD334" & aDNA_ID!="COD338" &
#                                            aDNA_ID!="COD339" & aDNA_ID!="COD341" & aDNA_ID!="COD082" &
#                                            aDNA_ID!="COD086" & aDNA_ID!="COD343" & aDNA_ID!="COD356" & 
#                                            aDNA_ID!="COD345" & aDNA_ID!="COD365" & aDNA_ID!="COD367" &
#                                            aDNA_ID!="COD377" & aDNA_ID!="COD383" & aDNA_ID!="COD389" &
#                                            aDNA_ID!="COD084" & aDNA_ID!="COD398" & aDNA_ID!="COD405" &
#                                            aDNA_ID!="COD421" & aDNA_ID!="COD408" & aDNA_ID!="COD411" &
#                                            aDNA_ID!="COD414" & aDNA_ID!="COD409") 
#dim(all.thre.isotope_genotype_table4) #9 individuals
#Only leaving NEA and Iceland Frontal
#all.thre.isotope_genotype_table5 <-subset(all.thre.isotope_genotype_table4, 
#                                          aDNA_ID!="COD372" & aDNA_ID!="COD394" & aDNA_ID!="COD089")
#dim(all.thre.isotope_genotype_table5) # 7 individuals with genotype - those that passed all 3 ratio's thresholds (C.N_ratio, C.S_ratio, N.S_ratio)

# PCA - Size - NEA/Iceland 
#iso_pca<-princomp(all.thre.isotope_genotype_table4[,c(8:11)], cor = FALSE, scores = TRUE)
#iso_pca$loadings
#Loadings:
#                         Comp.1 Comp.2 Comp.3 Comp.4
#d13C                            0.214  0.246  0.945
#d15N                            0.102  0.955 -0.274
#Meand2HV.SMOW.corrected  0.996                     
#Meand34SV.corrected            -0.968  0.159  0.181
#                 Comp.1 Comp.2 Comp.3 Comp.4
#SS loadings      1.00   1.00   1.00   1.00
#Proportion Var   0.25   0.25   0.25   0.25
#Cumulative Var   0.25   0.50   0.75   1.00

# eigen values
#a<-fviz_eig(iso_pca) # from the package "factoextra"
#a$data
#   dim         eig
#1   1 98.46744476
#2   2  1.18488441
#3   3  0.29341729
#4   4  0.05425354

# merging scores with descriptive info
#iso<-cbind(all.thre.isotope_genotype_table4[,c(1:7,27)],iso_pca$scores)
#head(iso)
# plot - PCs vs NEA/Iceland
#p.iso1<-iso%>%
#  ggplot(aes(Comp.1,Comp.2)) + geom_point(size = 4, aes(fill = factor(population.assigment), colour = factor(population.assigment)))+
#  geom_text(aes(label=Old_ID),size=2,hjust=0.5,vjust=1.8)+ 
  #stat_ellipse(geom="polygon", aes(fill = population.assigment), alpha = 0.2, show.legend = FALSE, level = 0.95) +
  #scale_color_manual(values=c("#F4C67D","#FFA500","#FA9837", "#C4C4C4", "666666", "949494")) +
#  xlab(paste0("PC1 (98.47%)")) +
#  ylab(paste0("PC2 (1.18%)")) +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#       panel.background = element_blank(), axis.line = element_line(colour = "black"))
#p.iso1
# contribution of each variable to the PCs
#variants<-fviz_pca_var(iso_pca)
#variants$data 
#                                          name          x          y       coord        cos2    contrib
#d13C                                       d13C  0.3731757  0.2637264   0.2088118   0.2088118  0.1640894
#d15N                                       d15N  0.5236123  0.1254575   0.2899095   0.2899095  0.2278179
#Meand2HV.SMOW.corrected Meand2HV.SMOW.corrected 11.1663198 -0.1001148 124.6967198 124.6967198 97.9897311
#Meand34SV.corrected         Meand34SV.corrected -0.8010459 -1.1907013   2.0594442   2.0594442  1.6183616

###############
#Time vs. Size #
################
dim(master_table) # 93individuals - all isotope data
head(all.thre_isotope_table)
dim(all.thre_isotope_table)# 64individuals - those that passed all 3 ratio's thresholds (C.N_ratio, C.S_ratio, N.S_ratio)
#p<-master_table %>%
#  mutate(Size = factor(Size, levels=c("50-80cm", "80-100cm", ">100cm"))) %>%
#  ggplot(aes(Size, fill=Layer)) +
#  geom_bar() + ggtitle("Oslo") + 
#  scale_color_manual(values=c("blue","paleturquoise3", "yellow3"))+
  #geom_text(aes(label=aDNA_ID),size=3,hjust=1,vjust=1) +
#  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#p

#Fisher test for Body Size across time
#NOT including: COD339, COD356, COD398, COD405
#TWO <- data.frame(
#  "1" = c(0, 8, 3),
#  "2" = c(12, 13, 4),
#  "3" = c(13, 9, 5),
#  "4" = c(3, 7, 9),
#  "5" = c(0, 2, 0),
#  "6" = c(0, 0, 1),
#  "7" = c(1, 1, 0),
#  row.names = c("Small", "Medium", "Big"),
#  stringsAsFactors = FALSE)
#TWO
#mosaicplot(TWO,
#           main = "Mosaic plot",
#           color = TRUE)
#testTWO <- fisher.test(TWO, simulate.p.value=TRUE)
#testTWO #p-value = 0.01399 with not simulated p value
#p-value = 0.01099 with simulated p value

#Fire (brann) Layer vs Size ####
master_table_2 <-subset(master_table, Layer2!="E" & Layer2!="F" & Layer2!="G" & Layer2!="H") #Fire layers 5, 6, 7, 8 were ignored because of small sample size for this test
dim(master_table_2) #86 individuals
#master_table_2 <-subset(C.N_table, Layer2!="E" & Layer2!="F" & Layer2!="G" & Layer2!="H")
#dim(master_table_2) #81 individuals ###### Kruskal-Wallis chi-squared = 12.796, df = 3, p-value = 0.0051, differences are still between B-D, and C-D (0.02) 
master_table_2 <-subset(all.thre_isotope_table, Layer2!="E" & Layer2!="F" & Layer2!="G" & Layer2!="H")
dim(master_table_2) #60 individuals ###### Kruskal-Wallis chi-squared = 9.3886, df = 3, p-value = 0.02455, differences are only between B-D (0.04).

#Checking normality
hist(master_table_2$Size2)
shapiro.test(master_table_2$Size2) #W = 0.80635, p-value = 2.896e-09 : significantly deviate from a normal distribution
qqnorm(master_table_2$Size2)
qqline(master_table_2$Size2) #Not normal distributed

#Kruskall-Wallis test
leveneTest(Size2~Layer2, data = master_table_2) #F-value = 1.8053; p-value = 0.1527 = equal variance
kruskal.test(Size2~Layer2, data = master_table_2) #Kruskal-Wallis chi-squared = 11.931, df = 3, p-value = 0.007623

#Dunn test
install.packages("FSA")
library(FSA)
i <- dunnTest(Size2~Layer2, data = master_table_2, method="holm")
i
#    Comparison       Z     P.unadj      P.adj
#1      A - B  2.10290634 0.035473958 0.10642187
#2      A - C  2.18129850 0.029161345 0.11664538
#3      B - C  0.13306842 0.894139269 1.00000000
#4      A - D -0.07883153 0.937166623 0.93716662
#5      B - D -2.62414768 0.008686611 0.04343306* 1125-1275 vs. 1275-1375
#6      C - D -2.70534368 0.006823375 0.04094025* 1250-1325 vs. 1275-1375

#Correlations between isotopes and Size
# Carbon vs. Size
dim(C.N_table) #N=86
hist(C.N_table$d13C)
shapiro.test(C.N_table$d13C) #W = 0.9817, p-value = 0.2636 = normal distribution
qqnorm(C.N_table$d13C)
qqline(C.N_table$d13C) #Normal distributed
#leveneTest(logd13C~Size, data = C.N_table) #F-value = 1.3407; p-value = 0.2673 = equal variance
leveneTest(d13C~Size, data = C.N_table) #F-value = 1.4968; p-value = 0.2298 = equal variance
#kruskal.test(d13C~Size, data = C.N_table) #Kruskal-Wallis chi-squared = 11.804, df = 2, p-value = 0.002734*
#Dunn test
#library(FSA)
#i <- dunnTest(logd13C~Size, data = C.N_table, method="bonferroni")
#i
#          Comparison         Z      P.unadj       P.adj
#1   >100cm - 50-80cm -3.435348 0.0005917927 0.001775378*
#2  >100cm - 80-100cm -1.913916 0.0556309308 0.166892793
#3 50-80cm - 80-100cm  1.890864 0.0586424619 0.175927386
size.model<-aov(data = C.N_table, formula = d13C~Size)
summary(size.model) # F-value=6.698 p.value= 0.00201 **
TukeyHSD(size.model)
#                        diff          lwr         upr     p adj
#50-80cm->100cm   -0.8732174 -1.44324659 -0.30318819 0.0012917*
#80-100cm->100cm  -0.4231121 -0.94433037  0.09810611 0.1346997
#80-100cm-50-80cm  0.4501053 -0.05795927  0.95816980 0.0931114

# Hydrogen vs. Size
C.N_table_2 <-subset(C.N_table, Old_ID!="2611" & Old_ID!="2620" & Old_ID!="2599") #Samples with NAs
dim(C.N_table_2) #83 individuals
hist(C.N_table_2$Meand2HV.SMOW.corrected)
shapiro.test(C.N_table_2$Meand2HV.SMOW.corrected) #W = 0.9904, p-value = 0.8021 = normal distribution
qqnorm(C.N_table_2$Meand2HV.SMOW.corrected)
qqline(C.N_table_2$Meand2HV.SMOW.corrected) #Normal distributed

leveneTest(Meand2HV.SMOW.corrected~Size, data = C.N_table_2) #F-value = 0.8861; p-value = 0.4163 = equal variance
size.model<-aov(data = C.N_table_2, formula = Meand2HV.SMOW.corrected~Size)
summary(size.model) # F-value=7.385 p.value= 0.00114 **
TukeyHSD(size.model)
#                      diff        lwr       upr     p adj
#50-80cm->100cm   -11.406094 -18.5606622 -4.251526 0.0007936*
#80-100cm->100cm   -5.025710 -11.6294146  1.577995 0.1704269
#80-100cm-50-80cm   6.380384   0.1227649 12.638003 0.0446526

#Nitrogen vs. Size
hist(C.N_table$d15N)
shapiro.test(C.N_table$d15N) #W = 0.98227, p-value = 0.2867 = normal distribution
qqnorm(C.N_table$d15N)
qqline(C.N_table$d15N) #Normal distributed

leveneTest(d15N~Size, data = C.N_table) #F-value = 2.4251; p-value = 0.09471 = equal variance
size.model<-aov(d15N~Size, data = C.N_table)
summary(size.model) # F-value=5.473 p.value= 0.00585 **
TukeyHSD(size.model)
#                      diff        lwr       upr     p adj
#50-80cm->100cm   -0.8172307 -1.4175098 -0.216951680 0.0047176*
#80-100cm->100cm  -0.5436034 -1.0924812  0.005274471 0.0528171
#80-100cm-50-80cm  0.2736274 -0.2613987  0.808653456 0.4445175

# Sulphur vs. Size
dim(all.thre_isotope_table) #N=64
hist(all.thre_isotope_table$Meand34SV.corrected)
shapiro.test(all.thre_isotope_table$Meand34SV.corrected) #W = 0.91774, p-value = 0.0004005 = not equal variance
qqnorm(all.thre_isotope_table$Meand34SV.corrected)
qqline(all.thre_isotope_table$Meand34SV.corrected) #Not normal distributed

#leveneTest(logd34S~Size, data = all.thre_isotope_table) #F-value = 1.0887; p-value = 0.3431 = equal variance
#kruskal.test(logd34S~Size, data = all.thre_isotope_table) #Kruskal-Wallis chi-squared = 0.1058, df = 2, p-value = 0.9485
leveneTest(Meand34SV.corrected~Size, data = all.thre_isotope_table) #F-value = 1.0887; p-value = 0.3431 = equal variance
kruskal.test(Meand34SV.corrected~Size, data = all.thre_isotope_table) #Kruskal-Wallis chi-squared = 0.1058, df = 2, p-value = 0.9485

#Plots isotope vs. size
install.packages("Rmisc")
library("Rmisc")
#Carbon
level_order <- c("50-80cm", "80-100cm", ">100cm")
#tgc <- summarySE(C.N_table, measurevar="d13C", groupvars="Size") #To get the ci
ggplot(data = C.N_table, aes(x=factor(Size, level = level_order), y=d13C, color = Size, fill = Size)) +
  #geom_point(size = 4, position=position_dodge(width=0.5)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") +
  #geom_errorbar(data = tgc,
  #              aes(ymin=d13C-ci, ymax=d13C+ci),
  #              width = 0.1, color = "black",
  #              #linetype = "dotted",
  #              position=position_dodge(width=0.5)) +
  #geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Hydrogen
tgc <- summarySE(C.N_table_2, measurevar="Meand2HV.SMOW.corrected", groupvars="Size") #To get the ci
ggplot(data = C.N_table_2, aes(x=factor(Size, level = level_order), y=Meand2HV.SMOW.corrected, color = Size, fill = Size)) +
  #geom_point(size = 4, position=position_dodge(width=0.5)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") +
  #geom_errorbar(data = tgc,
  #              aes(ymin=Meand2HV.SMOW.corrected-ci, ymax=Meand2HV.SMOW.corrected+ci),
  #              width = 0.1, color = "black",
  #              #linetype = "dotted",
  #              position=position_dodge(width=0.5)) +
  #geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Nitrogen
tgc <- summarySE(C.N_table, measurevar="d15N", groupvars="Size") #To get the ci
ggplot(data = C.N_table_2, aes(x=factor(Size, level = level_order), y=d15N, color = Size, fill = Size)) +
  #geom_point(size = 4, position=position_dodge(width=0.5)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") +
  #geom_errorbar(data = tgc,
  #              aes(ymin=d15N-ci, ymax=d15N+ci),
  #              width = 0.1, color = "black",
  #              #linetype = "dotted",
  #              position=position_dodge(width=0.5)) +
  #geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

!!!!!!#Sulphur ERROR!!!!!!!
tgc <- summarySE(all.thre_isotope_table, measurevar="Meand34SV.corrected", groupvars="Size") #To get the ci
ggplot(data = C.N_table_2, aes(x=factor(Size, level = level_order), y=Meand34SV.corrected, color = Size, fill = Size)) +
  #geom_point(size = 4, position=position_dodge(width=0.5)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") +
  #geom_errorbar(data = tgc,
  #              aes(ymin=Meand34SV.corrected-ci, ymax=Meand34SV.corrected+ci),
  #              width = 0.1, color = "black",
  #              #linetype = "dotted",
  #              position=position_dodge(width=0.5)) +
  #geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Correlations between isotopes and Time
#Removing Brann 5,6,7,8
C.N_table_timeupdated <-subset(C.N_table, Layer2!="E" & Layer2!="F" & Layer2!="G" & Layer2!="H")
C.N_table_timeupdated_2 <-subset(C.N_table_timeupdated, Old_ID!="2611" & Old_ID!="2620") #Samples with NAs
dim(C.N_table_timeupdated) #81 individuals
dim(C.N_table_timeupdated_2) #79 individuals
all.thre_isotope_table_timeupdated <-subset(all.thre_isotope_table, Layer2!="E" & Layer2!="F" & Layer2!="G" & Layer2!="H")
dim(all.thre_isotope_table_timeupdated) #60 individuals

# Carbon vs. Time
hist(C.N_table_timeupdated$d13C)
shapiro.test(C.N_table_timeupdated$d13C) #W = 0.97701, p-value = 0.1532 = normal distribution
qqnorm(C.N_table_timeupdated$d13C)
qqline(C.N_table_timeupdated$d13C) #Normal distributed

leveneTest(d13C~Layer2, data = C.N_table_timeupdated) #F-value = 0.8259; p-value = 0.4836 = equal variance
#kruskal.test(d13C~Layer2, data = C.N_table_timeupdated) #Kruskal-Wallis chi-squared = 13.902, df = 3, p-value = 0.003042*
#Dunn test
#library(FSA)
#i <- dunnTest(logd13C~Layer2, data = C.N_table_timeupdated, method="holm")
#i
#       Comparison    Z      P.unadj       P.adj
#1      A - B -1.5713890 0.1160923274 0.464369310
#2      A - C -1.1388613 0.2547609905 0.509521981
#3      B - C  0.5490634 0.5829619139 0.582961914
#4      A - D  1.3757416 0.1689016454 0.506704936
#5      B - D  3.5036329 0.0004589576 0.002753745*
#6      C - D  2.9850099 0.0028356917 0.014178459*
size.model<-aov(data = C.N_table_timeupdated, formula = d13C~Layer2)
summary(size.model) # F-value=5.911, p.value= 0.0011 **
TukeyHSD(size.model)
#           diff        lwr       upr     p adj
#B-A -0.46259740 -1.2333224 0.3081276 0.3981527
#C-A -0.40891608 -1.1879571 0.3701250 0.5166234
#D-A  0.53204545 -0.3162908 1.3803817 0.3588437
#C-B  0.05368132 -0.5362123 0.6435749 0.9951593
#D-B  0.99464286  0.3158616 1.6734241 0.0013697* 1125-1275 vs. 1275-1375
#D-C  0.94096154  0.2527523 1.6291707 0.0031839* 1250-1325 vs. 1275-1375

# Hydrogen vs. Time
hist(C.N_table_timeupdated_2$Meand2HV.SMOW.corrected)
shapiro.test(C.N_table_timeupdated_2$Meand2HV.SMOW.corrected) #W = 0.99087, p-value = 0.8522 = normal distribution
qqnorm(C.N_table_timeupdated_2$Meand2HV.SMOW.corrected)
qqline(C.N_table_timeupdated_2$Meand2HV.SMOW.corrected) #Normal distributed

leveneTest(Meand2HV.SMOW.corrected~Layer2, data = C.N_table_timeupdated_2) #F-value = 0.5106; p-value = 0.6762 = equal variance
size.model<-aov(data = C.N_table_timeupdated_2, formula = Meand2HV.SMOW.corrected~Layer2)
summary(size.model) # F-value=3.675, p.value= 0.0158 *
TukeyHSD(size.model)
#         diff        lwr       upr     p adj
#B-A -4.889523 -14.647996  4.868949 0.5554066
#C-A -1.446189 -11.309955  8.417576 0.9804176
#D-A  6.375464  -4.673830 17.424758 0.4330903
#C-B  3.443334  -4.025557 10.912224 0.6216293
#D-B 11.264987   2.288498 20.241476 0.0079705* 1125-1275 vs. 1275-1375
#D-C  7.821654  -1.269190 16.912497 0.1166671

#Nitrogen vs. Time
hist(C.N_table_timeupdated$d15N)
shapiro.test(C.N_table_timeupdated$d15N) #W = 0.9777, p-value = 0.1695 = normal distribution
qqnorm(C.N_table_timeupdated$d15N)
qqline(C.N_table_timeupdated$d15N) #Normal distributed

leveneTest(d15N~Layer2, data = C.N_table_timeupdated) #F-value = 2.3967; p-value = 0.07456 = equal variance
size.model<-aov(d15N~Layer2, data = C.N_table_timeupdated)
summary(size.model) # F-value=3.943 p.value= 0.0113 *
TukeyHSD(size.model)
#       diff         lwr      upr     p adj
#B-A 0.28198377 -0.54898528 1.112953 0.8094646
#C-A 0.37409091 -0.46584422 1.214026 0.6477246
#D-A 1.07346591  0.15881905 1.988113 0.0148278* 1000-1175 vs. 1275-1375
#C-B 0.09210714 -0.54389575 0.728110 0.9811375
#D-B 0.79148214  0.05964368 1.523321 0.0288590* 1125-1275 vs. 1275-1375
#D-C 0.69937500 -0.04262838 1.441378 0.0719385

# Sulphur vs. Time
hist(all.thre_isotope_table_timeupdated$Meand34SV.corrected)
shapiro.test(all.thre_isotope_table_timeupdated$Meand34SV.corrected) #W = 0.94073, p-value = 0.005837 = NOT normal distribution
qqnorm(all.thre_isotope_table_timeupdated$Meand34SV.corrected)
qqline(all.thre_isotope_table_timeupdated$Meand34SV.corrected) #Normal distributed

leveneTest(Meand34SV.corrected~Layer2, data = all.thre_isotope_table_timeupdated) #F-value = 2.0731; p-value = 0.1141 = equal variance
kruskal.test(Meand34SV.corrected~Layer2, data = all.thre_isotope_table_timeupdated) #Kruskal-Wallis chi-squared = 4.6881, df = 3, p-value = 0.1961

#Plots isotope vs. time
library("Rmisc")
#Carbon
#tgc <- summarySE(C.N_table_timeupdated, measurevar="d13C", groupvars="Layer2") #To get the ci
ggplot(data = C.N_table_timeupdated, aes(Layer2, d13C, color = Layer2, fill = Layer2)) +
  #geom_point(size = 4, position=position_dodge(width=0.5)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") +
  #geom_errorbar(data = tgc,
  #              aes(ymin=d13C-ci, ymax=d13C+ci),
  #              width = 0.1, color = "black",
  #              #linetype = "dotted",
  #              position=position_dodge(width=0.5)) +
  #geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Hydrogen
#tgc <- summarySE(C.N_table_timeupdated_2, measurevar="Meand2HV.SMOW.corrected", groupvars="Layer2") #To get the ci
ggplot(data = C.N_table_timeupdated_2, aes(Layer2, Meand2HV.SMOW.corrected, color = Layer2, fill = Layer2)) +
  #geom_point(size = 4, position=position_dodge(width=0.5)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") +
  #geom_errorbar(data = tgc,
  #              aes(ymin=Meand2HV.SMOW.corrected-ci, ymax=Meand2HV.SMOW.corrected+ci),
  #              width = 0.1, color = "black",
  #              #linetype = "dotted",
  #              position=position_dodge(width=0.5)) +
  #geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Nitrogen
#tgc <- summarySE(C.N_table_timeupdated, measurevar="d15N", groupvars="Layer2") #To get the ci
ggplot(data = C.N_table_timeupdated, aes(Layer2, d15N, color = Layer2, fill = Layer2)) +
  #geom_point(size = 4, position=position_dodge(width=0.5)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") +
  #geom_errorbar(data = tgc,
  #              aes(ymin=d15N-ci, ymax=d15N+ci),
  #              width = 0.1, color = "black",
  #              #linetype = "dotted",
  #              position=position_dodge(width=0.5)) +
  #geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Sulphur
#tgc <- summarySE(all.thre_isotope_table_timeupdated, measurevar="Meand34SV.corrected", groupvars="Layer2") #To get the ci
ggplot(data = all.thre_isotope_table_timeupdated, aes(Layer2, Meand34SV.corrected, color = Layer2, fill = Layer2)) +
  #geom_point(size = 4, position=position_dodge(width=0.5)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") +
  #geom_errorbar(data = tgc,
  #              aes(ymin=Meand34SV.corrected-ci, ymax=Meand34SV.corrected+ci),
  #              width = 0.1, color = "black",
  #              #linetype = "dotted",
  #              position=position_dodge(width=0.5)) +
  #geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#######################################
# Bimodal distribution Carbon and Size#
#######################################
# Density Plot
head(C.N_table) # 86 individuals
#Subset only size 80-100 cm with the Carbon values
medsizeC <- C.N_table[C.N_table$Size == '80-100cm', ]
dim(medsizeC) #38 ind
grasizeC <- C.N_table[C.N_table$Size == '>100cm', ]
dim(grasizeC) #23 ind
chisizeC <- C.N_table[C.N_table$Size == '50-80cm', ]
dim(chisizeC) #25 ind
ggplot(chisizeC, 
       aes(x = d13C, 
           fill = Size)) +
  geom_density(alpha = 0.2) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

#########################################
# Bimodal distribution Nitrogen and Size#
#########################################
head(C.N_table) # 86 ind
medsizeN <- C.N_table[C.N_table$Size == '80-100cm', ]
dim(medsizeN) #38 ind
grasizeN <- C.N_table[C.N_table$Size == '>100cm', ]
dim(grasizeN) #23 ind
chisizeN <- C.N_table[C.N_table$Size == '50-80cm', ]
dim(chisizeN) #25 ind
ggplot(chisizeN, 
       aes(x = d15N, 
           fill = Size)) +
  geom_density(alpha = 0.4) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

#########################################
# Bimodal distribution Hydrogen and Size#
#########################################
head(C.N_table) # 64 ind
medsizeH <- C.N_table[C.N_table$Size == '80-100cm', ]
dim(medsizeH) #38 ind
grasizeH <- C.N_table[C.N_table$Size == '>100cm', ]
dim(grasizeH) #23 ind
chisizeH <- C.N_table[C.N_table$Size == '50-80cm', ]
dim(chisizeH) #25 ind
ggplot(chisizeH, 
       aes(x = Meand2HV.SMOW.corrected, 
           fill = Size)) +
  geom_density(alpha = 0.4) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())


######################
#Isotope Correlations#
######################
dim(master_table) # 93individuals - all isotope data
dim(C.N_table) #86 individuals
dim(all.thre_isotope_table) #64 individuals

# Nitrogen and Carbon
cor.test(C.N_table$d13C, C.N_table$d15N, method = "pearson") #t = 4.7221, df = 84, p-value = 9.246e-06, cor = 0.4580078

p<-C.N_table %>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  ggplot(aes(d13C, d15N)) +
  geom_point(aes(size=Size)) + ggtitle("Oslo") + 
  #scale_color_manual(values=c("blue","paleturquoise3", "yellow3"))+
  #geom_text(aes(label=aDNA_ID),size=3,hjust=1,vjust=1) +
  geom_smooth(method=lm) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p


# Hydrogen and Nitrogen
cor.test(C.N_table$d15N, C.N_table$Meand2HV.SMOW.corrected, method = "pearson") # t = 5.3676, df = 81, p-value = 7.424e-07, corr=0.5122222
p<-all.thre_isotope_table %>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  ggplot(aes(d15N,Meand2HV.SMOW.corrected)) +
  geom_point(aes(size=Size)) + ggtitle("Oslo") + 
  #scale_color_manual(values=c("blue","paleturquoise3", "yellow3"))+
  #geom_text(aes(label=aDNA_ID),size=3,hjust=1,vjust=1) +
  geom_smooth(method=lm) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

# Sulfur and Carbon
hist(all.thre_isotope_table$d13C)
shapiro.test(all.thre_isotope_table$d13C) #W = 0.95804, p-value = 0.02896 = not normal distribution
qqnorm(all.thre_isotope_table$d13C)
qqline(all.thre_isotope_table$d13C)

cor.test(all.thre_isotope_table$d13C, all.thre_isotope_table$Meand34SV.corrected, method = "spearman", exact = FALSE) # S = 51636, p-value = 0.1497, corr = -0.182145

p<-all.thre_isotope_table %>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  ggplot(aes(d13C, Meand34SV.corrected)) +
  geom_point(aes(size=Size)) + ggtitle("Oslo") + 
  #scale_color_manual(values=c("blue","paleturquoise3", "yellow3"))+
  #geom_text(aes(label=aDNA_ID),size=3,hjust=1,vjust=1) +
  geom_smooth(method=lm) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

#######################
#Isotopes vs Genotypes#
#######################
#Genotypes of these individuals were removed because low endogenous DNA and/or SNPs supporting population assignments: COD348 (2535), COD379 (2570), COD420 (2612), COD415 (2607)
dim(master_table) # 93individuals - all isotope data
dim(C.N_table) # 86 individuals - For only Carbon and Nitrogen data (passed the C.N_ratio threshold)
dim(all.thre_isotope_table) # 64individuals - those that passed all 3 ratio's thresholds (C.N_ratio, C.S_ratio, N.S_ratio)

all.thre.isotope_genotype_table1 <- subset(master_table, C.N_ratio>2.9& C.N_ratio<3.6 & C.S_ratio>125 & 
                                             C.S_ratio<225 & N.S_ratio>40 & N.S_ratio<80,
                                           select = c(aDNA_ID,Old_ID ,Size,Element,Layer,Layer2,Age,d13C,d15N,Meand2HV.SMOW.corrected,Meand34SV.corrected,logd13C,logd34S,
                                                      LG01,LofotenSkrei, LofotenCoastal, NorthSea,  Oresund, IrishSea, Iceland, IcFrontal, IcCoastal, WestCoast, northcentral, northernmost,
                                                      population.assigment, ncno, C.N_ratio)) 

all.thre.isotope_genotype_table2<-na.omit(all.thre.isotope_genotype_table1) #20 individuals
dim(all.thre.isotope_genotype_table2) #20 individuals
all.thre.isotope_genotype_table3<-subset(all.thre.isotope_genotype_table2,Old_ID!= 2535 &
                                           Old_ID!=2570 & 
                                           Old_ID!=2612 &
                                           aDNA_ID!="COD408") #COD421 does not count with isotope data. Still 20 individuals. We removed everything needed by applying the NA filter

dim(all.thre.isotope_genotype_table3) # 20 individuals with genotype - those that passed all 3 ratio's thresholds (C.N_ratio, C.S_ratio, N.S_ratio)

C.N.isotope_genotype_table1 <-subset(master_table, C.N_ratio>2.9& C.N_ratio<3.6,
                                     select = c(aDNA_ID,Old_ID ,Size,Element,Layer,Layer2,Age,d13C,d15N,Meand2HV.SMOW.corrected,Meand34SV.corrected,logd13C,logd34S,
                                                LG01,LofotenSkrei, LofotenCoastal, NorthSea,  Oresund, IrishSea, Iceland, IcFrontal, IcCoastal, WestCoast, northcentral, northernmost,
                                                population.assigment,ncno))
dim(C.N.isotope_genotype_table1) #86 individuals
C.N.isotope_genotype_table2<-na.omit(C.N.isotope_genotype_table1)
dim(C.N.isotope_genotype_table2) #25 individuals
C.N.isotope_genotype_table3<-subset(C.N.isotope_genotype_table2,Old_ID!= 2535 &
                                      Old_ID!=2570 & 
                                      Old_ID!=2612 &
                                      aDNA_ID!="COD408") #25 individuals. We removed everything needed by applying the NA filter
dim(C.N.isotope_genotype_table3) # 25 individuals with genotype - For only Carbon and Nitrogen data (passed the C.N_ratio threshold)
dim(all.thre.isotope_genotype_table3) # 20 individuals with genotype - those that passed all 3 ratio's thresholds (C.N_ratio, C.S_ratio, N.S_ratio)

#Logistically transformation of NEA and Iceland
C.N.isotope_genotype_table3$LofotenSkrei_log <- logit(C.N.isotope_genotype_table3$LofotenSkrei) #In logit(probNEA) : proportions remapped to (0.025, 0.975)
C.N.isotope_genotype_table3$Iceland_log <- logit(C.N.isotope_genotype_table3$Iceland)
C.N.isotope_genotype_table3$IcFrontal_log <- logit(C.N.isotope_genotype_table3$IcFrontal)
all.thre.isotope_genotype_table3$LofotenSkrei_log <- logit(all.thre.isotope_genotype_table3$LofotenSkrei) #In logit(probNEA) : proportions remapped to (0.025, 0.975)
all.thre.isotope_genotype_table3$Iceland_log <- logit(all.thre.isotope_genotype_table3$Iceland)
all.thre.isotope_genotype_table3$IcFrontal_log <- logit(all.thre.isotope_genotype_table3$IcFrontal)

#################################################################
#######################Linear regression#########################
#Based on: https://www.datacamp.com/tutorial/linear-regression-R#
################################################################
hist(all.thre.isotope_genotype_table3$LofotenSkrei_log)
shapiro.test(all.thre.isotope_genotype_table3$LofotenSkrei_log) #W = 0.78024, p-value = 0.0004448 :  significantly deviate from a normal distribution
qqnorm(all.thre.isotope_genotype_table3$LofotenSkrei)
qqline(all.thre.isotope_genotype_table3$LofotenSkrei)

hist(all.thre.isotope_genotype_table3$Iceland_log)
shapiro.test(all.thre.isotope_genotype_table3$Iceland_log) #W = 0.81719, p-value = 0.02388 :  significantly deviate from a normal distribution
qqnorm(all.thre.isotope_genotype_table3$Iceland_log)
qqline(all.thre.isotope_genotype_table3$Iceland_log)

hist(all.thre.isotope_genotype_table3$IcFrontal_log)
shapiro.test(all.thre.isotope_genotype_table3$IcFrontal_log) #W = 0.81784, p-value = 0.001614 :  significantly deviate from a normal distribution
qqnorm(all.thre.isotope_genotype_table3$IcFrontal_log)
qqline(all.thre.isotope_genotype_table3$IcFrontal_log)

#Model (with stats package)
#NEA
model <- lm(LofotenSkrei_log~d13C+Meand34SV.corrected+Meand2HV.SMOW.corrected+d15N, data = all.thre.isotope_genotype_table3)
summary(model)
#The sum of my residuals is = 0.0873, so it is approximately zero :) 
#Also, our model is fairly symmetrical because our median is close to zero too
#If using log values of C and S, my residuals are = 0.0098
# Residuals:
#Min      1Q  Median      3Q     Max 
#-2.2616 -1.0107  0.1158  0.8954  2.3484 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)             -26.59200   13.70426  -1.940   0.0714 .
#d13C                     -1.50430    0.55213  -2.725   0.0157 *
#  Meand34SV.corrected       0.34171    0.24087   1.419   0.1764  
#Meand2HV.SMOW.corrected  -0.01754    0.04118  -0.426   0.6762  
#d15N                     -0.02802    0.50592  -0.055   0.9566  
#Residual standard error: 1.409 on 15 degrees of freedom
#Multiple R-squared:  0.6036,	Adjusted R-squared:  0.4978 
#F-statistic: 5.709 on 4 and 15 DF,  p-value: 0.005356
plot(model$residuals, pch = 16, col = "red") #Not seeing any clear patterns on your residuals is good
#Equal variance test if using lm
ncvTest(model) #Chisquare = 0.3189053, Df = 1, p = 0.57227 = we assume that the residuals are homoscedastic.
#To avoid Multicolinearity I take away Nitrogen
#model <- lm(LofotenSkrei_log~d13C+Meand34SV.corrected+Meand2HV.SMOW.corrected, data = all.thre.isotope_genotype_table3)
#summary(model)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-2.2352 -0.9957  0.1172  0.8932  2.3530 
#Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)   
#(Intercept)             -27.23007    7.18870  -3.788  0.00161 **
#  d13C                     -1.51485    0.50183  -3.019  0.00816 **
#  Meand34SV.corrected       0.34734    0.21149   1.642  0.12003   
#Meand2HV.SMOW.corrected  -0.01780    0.03961  -0.449  0.65914   
#Residual standard error: 1.365 on 16 degrees of freedom
#Multiple R-squared:  0.6035,	Adjusted R-squared:  0.5291 
#F-statistic: 8.117 on 3 and 16 DF,  p-value: 0.001641
#plot(model$residuals, pch = 16, col = "red") #Not seeing any clear patterns on your residuals is good
#Equal variance test if using lm
#ncvTest(model) #Chisquare = 0.3122428, Df = 1, p = 0.57631 = we assume that the residuals are homoscedastic.

#Individual plot for NEA and C
p<-all.thre.isotope_genotype_table3%>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  ggplot(aes(x=LofotenSkrei_log, y=d13C, color = LofotenSkrei_log, fill = LofotenSkrei_log)) +
  geom_point(aes(size=Size)) + 
  #geom_violin(trim= FALSE, alpha=0.6) +
  #geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  #geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Iceland
model <- lm(Iceland_log~d13C+Meand34SV.corrected+Meand2HV.SMOW.corrected+d15N, data = all.thre.isotope_genotype_table3)
summary(model)
#The sum of my residuals is =-0-0041, so it is approximately zero :) 
#Also, our model is fairly symmetrical because our median is close to zero too
#Residuals:
#   Min      1Q  Median      3Q     Max 
#-6.6723 -1.3899  0.1419  2.1502  5.7531 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)              6.68147   32.41463   0.206    0.839
#d13C                    -1.22942    1.30595  -0.941    0.361
#Meand34SV.corrected     -0.01563    0.56973  -0.027    0.978
#Meand2HV.SMOW.corrected  0.06682    0.09741   0.686    0.503
#d15N                    -1.77778    1.19665  -1.486    0.158
#Residual standard error: 3.333 on 15 degrees of freedom
#Multiple R-squared:  0.2868,	Adjusted R-squared:  0.09667 
#F-statistic: 1.508 on 4 and 15 DF,  p-value: 0.2497
plot(model$residuals, pch = 16, col = "red") #Not seeing any clear patterns on your residuals is good
#Equal variance test if using lm
ncvTest(model) #CChisquare = 2.904987, Df = 1, p = 0.088306 = we assume that the residuals are homoscedastic.

#To avoid Multicolinearity I take away Nitrogen
#model <- lm(Iceland_log~d13C+Meand34SV.corrected+Meand2HV.SMOW.corrected, data = all.thre.isotope_genotype_table3)
#summary(model)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-7.3184 -1.6401  0.3287  2.1603  4.7495 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)             -33.79701   18.20958  -1.856    0.082 .
#d13C                     -1.89870    1.27118  -1.494    0.155  
#Meand34SV.corrected       0.34129    0.53572   0.637    0.533  
#Meand2HV.SMOW.corrected   0.05014    0.10034   0.500    0.624  
#Residual standard error: 3.456 on 16 degrees of freedom
#Multiple R-squared:  0.1819,	Adjusted R-squared:  0.02852 
#F-statistic: 1.186 on 3 and 16 DF,  p-value: 0.3464
#plot(model$residuals, pch = 16, col = "red") #Not seeing any clear patterns on your residuals is good
#Equal variance test if using lm
#ncvTest(model) #Chisquare = 0.5868632, Df = 1, p = 0.44363 = we assume that the residuals are homoscedastic.

#Individual plot for Iceland and C
p<-all.thre.isotope_genotype_table3%>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  ggplot(aes(x=Iceland_log, y=d13C)) +
  geom_point(aes(size=Size)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

#Iceland Frontal
model <- lm(IcFrontal_log~d13C+Meand34SV.corrected+Meand2HV.SMOW.corrected+d15N, data = all.thre.isotope_genotype_table3)
summary(model)
#The sum of my residuals is =-0-0041, so it is approximately zero :) 
#Also, our model is fairly symmetrical because our median is close to zero too
#Residuals:
#Min      1Q  Median      3Q     Max 
#-8.8537 -1.7331 -0.4159  3.0263  7.6154 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)
#(Intercept)             -7.64222   44.51853  -0.172    0.866
#d13C                    -2.38587    1.79360  -1.330    0.203
#Meand34SV.corrected      0.26672    0.78247   0.341    0.738
#Meand2HV.SMOW.corrected  0.06245    0.13378   0.467    0.647
#d15N                    -2.30440    1.64349  -1.402    0.181
#Residual standard error: 4.578 on 15 degrees of freedom
#Multiple R-squared:  0.3861,	Adjusted R-squared:  0.2224 
#F-statistic: 2.359 on 4 and 15 DF,  p-value: 0.1003
plot(model$residuals, pch = 16, col = "red") #Not seeing any clear patterns on your residuals is good
#Equal variance test if using lm
ncvTest(model) #Chisquare = 3.2636, Df = 1, p = 0.070833 = we assume that the residuals are homoscedastic.

#To avoid Multicolinearity I take away Nitrogen
#model <- lm(IcFrontal_log~d13C+Meand34SV.corrected+Meand2HV.SMOW.corrected, data = all.thre.isotope_genotype_table3)
#summary(model)
#Residuals:
#Min     1Q Median     3Q    Max 
#-9.691 -2.094  1.175  3.145  6.314 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)             -60.11133   24.83336  -2.421   0.0278 *
#  d13C                     -3.25340    1.73358  -1.877   0.0789 .
#Meand34SV.corrected       0.72938    0.73060   0.998   0.3330  
#Meand2HV.SMOW.corrected   0.04083    0.13684   0.298   0.7693  
#Residual standard error: 4.714 on 16 degrees of freedom
#Multiple R-squared:  0.3057,	Adjusted R-squared:  0.1755 
#F-statistic: 2.348 on 3 and 16 DF,  p-value: 0.1112
#plot(model$residuals, pch = 16, col = "red") #Not seeing any clear patterns on your residuals is good
#Equal variance test if using lm
#ncvTest(model) #Chisquare = 0.8666083, Df = 1, p = 0.3519 = we assume that the residuals are homoscedastic.

#Individual plot for Iceland Frontal and C
p<-all.thre.isotope_genotype_table3%>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  ggplot(aes(x=d13C, y=IcFrontal_log)) +
  geom_point(aes(size=Size)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

#Overall graph including NEA, Iceland and Iceland Frontal against Carbon values
library(lattice)
###Original values
#Option one
xyplot(d13C~LofotenSkrei_log + Iceland_log + IcFrontal_log, all.thre.isotope_genotype_table3,
#       groups = LofotenSkrei + Iceland + IcFrontal,
        pch=19, cex = 2,
        xlab = "Genetic Probability",
#       auto.key=TRUE,
        auto.key=list(x=0.75,y=0.95, text=c("LofotenSkrei","Iceland", "IcelandFrontal")))

#Option two (using this one)
all.thre.isotope_genotype_table3 %>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  pivot_longer(c(LofotenSkrei,IcFrontal, Iceland), names_to = "pop", values_to = "response") %>%
  ggplot(aes(y = d13C, x = response,
             color = pop, fill = pop)) +
  geom_point(aes(size=Size, stroke = 0)) + theme(legend.position="top") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

##Log values
#Option one
xyplot(d13C~LofotenSkrei_log + Iceland_log + IcFrontal_log,all.thre.isotope_genotype_table3,
       pch=19, cex = 2,
       xlab = "Genetic Probability",
       auto.key=list(x=0.05,y=0.95, text=c("LofotenSkrei","Iceland", "IcelandFrontal")))

#Option two (using this one)
#library(ggbreak)
all.thre.isotope_genotype_table3 %>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  pivot_longer(LofotenSkrei_log:IcFrontal_log, names_to = "pop", values_to = "response") %>%
  ggplot(aes(y = d13C, x = response,
             color = pop, fill = pop)) +
  geom_point(aes(size=Size, stroke = 0)) + theme(legend.position="top") +
  #scale_x_break(c(-15, -12), scales = 5) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


#################################
#Individual isotopes vs Genotype#
#################################
# Carbon
dim(C.N.isotope_genotype_table3) # 25 individuals
#Removing individuals with a source population with less than 70% proability
#NOTE that individuals from Oslo Mindets Tomt have no isotope data.
C.N.isotope_genotype_table4 <-subset(C.N.isotope_genotype_table3, 
                                     aDNA_ID!="COD332" & aDNA_ID!="COD334" & aDNA_ID!="COD338" &
                                       aDNA_ID!="COD339" & aDNA_ID!="COD341" & aDNA_ID!="COD082" &
                                       aDNA_ID!="COD086" & aDNA_ID!="COD343" &
                                       aDNA_ID!="COD356" & aDNA_ID!="COD367" &
                                       aDNA_ID!="COD365" &
                                       aDNA_ID!="COD377" & aDNA_ID!="COD383" & aDNA_ID!="COD389" &
                                       aDNA_ID!="COD084" & aDNA_ID!="COD398" & aDNA_ID!="COD405" &
                                       aDNA_ID!="COD421" & aDNA_ID!="COD408" & aDNA_ID!="COD411" &
                                       aDNA_ID!="COD414" & aDNA_ID!="COD409") 
#aDNA_ID!="COD372" -> outlier in Iceland
dim(C.N.isotope_genotype_table4) #10 individuals
#Including just NEA and Iceland Frontal
C.N.isotope_genotype_table5 <-subset(C.N.isotope_genotype_table4, 
                                     aDNA_ID!="COD372" & aDNA_ID!="COD394" & aDNA_ID!="COD089")
dim(C.N.isotope_genotype_table5) #8 individuals

hist(C.N.isotope_genotype_table4$d13C)
shapiro.test(C.N.isotope_genotype_table4$d13C) #W = 0.95658, p-value = 0.7463 = normal distribution
qqnorm(C.N.isotope_genotype_table4$d13C)
qqline(C.N.isotope_genotype_table4$d13C)
#leveneTest(logd13C~population.assigment, data = C.N.isotope_genotype_table4) #F-value = 1.0284; p-value = 0.3402 = equal variance
bartlett.test(d13C~population.assigment, data = C.N.isotope_genotype_table4) #Bartlett's K-squared = 0.68524, df = 1, p-value = 0.4078
#kruskal.test(logd13C~population.assigment, data = C.N.isotope_genotype_table4) #Kruskal-Wallis chi-squared = 4.6881, df = 3, p-value = 0.1961
size.model<-aov(formula = d13C~population.assigment, data = C.N.isotope_genotype_table4)
summary(size.model)
#                     Df Sum Sq Mean Sq F value Pr(>F)  
#population.assigment  1  1.154  1.1537   5.986 0.0402 *
#Residuals             8  1.542  0.1927  
TukeyHSD(size.model)
#  Tukey multiple comparisons of means
#95% family-wise confidence level
#Fit: aov(formula = d13C ~ population.assigment, data = C.N.isotope_genotype_table4)
#$population.assigment
#diff       lwr         upr     p adj
#NEA-IC -0.6933333 -1.346828 -0.03983834 0.0401502

p<-C.N.isotope_genotype_table4%>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  ggplot(aes(x=population.assigment, y=d13C, group=population.assigment)) +
  geom_point(aes(size=Size)) + 
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  #geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

#Hydrogen
hist(C.N.isotope_genotype_table4$Meand2HV.SMOW.corrected)
shapiro.test(C.N.isotope_genotype_table4$Meand2HV.SMOW.corrected) #W = 0.86935, p-value = 0.09824 = normal distribution
qqnorm(C.N.isotope_genotype_table4$Meand2HV.SMOW.corrected)
qqline(C.N.isotope_genotype_table4$Meand2HV.SMOW.corrected)
#leveneTest(Meand2HV.SMOW.corrected~population.assigment, data = C.N.isotope_genotype_table4) #F-value = 8.1514; p-value = 0.02131 * = not an equal variance
bartlett.test(Meand2HV.SMOW.corrected~population.assigment, data = C.N.isotope_genotype_table4) #Bartlett's K-squared = 4.2239, df = 1, p-value = 0.03986 = unequal variance
size.model <- oneway.test(formula = Meand2HV.SMOW.corrected~population.assigment, data = C.N.isotope_genotype_table4, var.equal = FALSE)
size.model #F = 0.61123, num df = 1.0000, denom df = 3.4151, p-value = 0.4849

p<-C.N.isotope_genotype_table4%>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  ggplot(aes(x=population.assigment, y=Meand2HV.SMOW.corrected, group=population.assigment)) +
  geom_point(aes(size=Size)) + 
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  #geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

#Nitrogen
hist(C.N.isotope_genotype_table4$d15N)
shapiro.test(C.N.isotope_genotype_table4$d15N) #W = 0.8674, p-value = 0.09321 = normal distribution
qqnorm(C.N.isotope_genotype_table4$d15N)
qqline(C.N.isotope_genotype_table4$d15N)
#leveneTest(d15N~population.assigment, data = C.N.isotope_genotype_table4) #F-value = 8.1514; p-value = 0.02131 * = not an equal variance
bartlett.test(d15N~population.assigment, data = C.N.isotope_genotype_table4) #Bartlett's K-squared = 2.1424, df = 1, p-value = 0.1433 = equal variance
size.model<-aov(formula = d15N~population.assigment, data = C.N.isotope_genotype_table4)
summary(size.model) #
#                      Df Sum Sq Mean Sq F value Pr(>F)
#population.assigment  1  0.310  0.3096   0.421  0.535
#Residuals             8  5.882  0.7352

p<-C.N.isotope_genotype_table4%>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  ggplot(aes(x=population.assigment, y=d15N, group=population.assigment)) +
  geom_point(aes(size=Size)) + 
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

#Sulphur
#Removing individuals with a source popualtions with less thab 70% probability.
all.thre.isotope_genotype_table4 <-subset(all.thre.isotope_genotype_table3, 
                                          aDNA_ID!="COD332" & aDNA_ID!="COD334" & aDNA_ID!="COD338" &
                                            aDNA_ID!="COD339" & aDNA_ID!="COD341" & aDNA_ID!="COD082" &
                                            aDNA_ID!="COD086" & aDNA_ID!="COD343" & aDNA_ID!="COD356" & 
                                            aDNA_ID!="COD365" & aDNA_ID!="COD367" &
                                            aDNA_ID!="COD377" & aDNA_ID!="COD383" & aDNA_ID!="COD389" &
                                            aDNA_ID!="COD084" & aDNA_ID!="COD398" & aDNA_ID!="COD405" &
                                            aDNA_ID!="COD421" & aDNA_ID!="COD408" & aDNA_ID!="COD411" &
                                            aDNA_ID!="COD414" & aDNA_ID!="COD409") 
dim(all.thre.isotope_genotype_table4) #9 individuals
#Only leaving NEA and Iceland Frontal
all.thre.isotope_genotype_table5 <-subset(all.thre.isotope_genotype_table4, 
                                          aDNA_ID!="COD372" & aDNA_ID!="COD394" & aDNA_ID!="COD089")
dim(all.thre.isotope_genotype_table5) #7 individuals
hist(all.thre.isotope_genotype_table4$Meand34SV.corrected)
shapiro.test(all.thre.isotope_genotype_table4$Meand34SV.corrected) #W = 0.93721, p-value = 0.5529 = normal distribution
qqnorm(all.thre.isotope_genotype_table4$Meand34SV.corrected)
qqline(all.thre.isotope_genotype_table4$Meand34SV.corrected)
#leveneTest(logd34S~population.assigment, data = all.thre.isotope_genotype_table4) #F-value = 8.1514; p-value = 0.02131 * = not an equal variance
bartlett.test(Meand34SV.corrected~population.assigment, data = all.thre.isotope_genotype_table4) #Bartlett's K-squared = 0.012509, df = 1, p-value = 0.9109 = equal variance
size.model<-aov(formula = Meand34SV.corrected~population.assigment, data = all.thre.isotope_genotype_table4)
summary(size.model)
#                     Df  Sum Sq  Mean Sq F value Pr(>F)
#population.assigment  1  4.987   4.987   2.557  0.154
#Residuals             7 13.653   1.950

p<-all.thre.isotope_genotype_table4%>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  ggplot(aes(x=population.assigment, y=Meand34SV.corrected, group=population.assigment)) +
  #geom_text(aes(label=Old_ID),size=2,hjust=0.5,vjust=1.8)+ 
  geom_point(aes(size=Size)) + 
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  #geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

##################################################
#Fisher test for Cranial vs. Postcranial evidence#
##################################################
#NOT including 50/50% combinations (COD082, COD084, COD339, COD356, COD398, COD405) and COD408 and COD421 because of they uncertain northernmost genotype
ONE <- data.frame(
  "Local" = c(3, 8),
  "Trade" = c(2, 14),
  row.names = c("Cranial", "Postcranial"),
  stringsAsFactors = FALSE)
ONE

mosaicplot(ONE,
           main = "Mosaic plot",
           color = TRUE)

testONE <- fisher.test(ONE) #p-value = 0.2806
testONE

#Plot
DATA <- c()
for (row in rownames(ONE)) {
  for (col in colnames(ONE)) {
    DATA <- rbind(DATA, matrix(rep(c(row, col), ONE[row, col]), ncol = 2, byrow = TRUE))
  }
}
DATABONES <- as.data.frame(DATA)
colnames(DATABONES) <- c("Element", "Fish_Exchange")
df
SelectedTest <- fisher.test(table(DATABONES))
SelectedTest

library(ggstatsplot)
ggbarstats(
  DATABONES, Element, Fish_Exchange,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(SelectedTest$p.value < 0.001, "< 0.001", round(SelectedTest$p.value, 3))
  )
)


###Correlation C/N ratio to d13C (with data set of 86 individuals)
head(master_table)
dim(master_table) #93 individuals
dim(C.N_table) #86 individuals - For only Carbon, Nitrogen and Hydrogen data (passed the C.N_ratio threshold)
dim(all.thre_isotope_table) #64 individuals but with all information

hist(C.N_table$C.N_ratio)
shapiro.test(C.N_table$C.N_ratio) #W = 0.93861, p-value = 0.0004916 = NOT normal distribution
qqnorm(C.N_table$C.N_ratio)
qqline(C.N_table$C.N_ratio) #Not normal distributed

cor.test(C.N_table$d13C, C.N_table$C.N_ratio, method = "spearman") #S = 126545, p-value = 0.07367, rho = -0.1938751
#cor.test(C.N_table$d15N, C.N_table$C.N_ratio, method = "spearman") #S = 73072, p-value = 0.003607*, rho = 0.3106065
#cor.test(C.N_table$d13C, C.N_table$C.N_ratio, method = "pearson") #t = -2.2576, df = 84, p-value = 0.02657*, cor = -0.2391721 

p<-C.N_table %>%
  mutate(Size = factor(Size, levels=c("50-80cm","80-100cm",">100cm"))) %>%
  ggplot(aes(C.N_ratio, d13C)) +
  geom_point(aes(size=Size)) + ggtitle("Oslo") + 
  #scale_color_manual(values=c("blue","paleturquoise3", "yellow3"))+
  #geom_text(aes(label=aDNA_ID),size=3,hjust=1,vjust=1) +
  geom_smooth(method=lm) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

min(C.N_table$C.N_ratio)
max(C.N_table$C.N_ratio)
###Differences between C/N ratio across time (4 Brann Layers, with data set of 86 individuals)
#Re opening table were we removed Brann 5,6,7,8
dim(C.N_table_timeupdated) #81 individuals

# Carbon vs. Time
hist(C.N_table_timeupdated$C.N_ratio)
shapiro.test(C.N_table_timeupdated$C.N_ratio) #W = 0.92554, p-value = 0.0001595 = NOT normal distribution
qqnorm(C.N_table_timeupdated$C.N_ratio)
qqline(C.N_table_timeupdated$C.N_ratio) #Not normal distributed

leveneTest(C.N_ratio~Layer2, data = C.N_table_timeupdated) #F-value = 3.8722; p-value = 0.01234 * = NOT equal variance
kruskal.test(C.N_ratio~Layer2, data = C.N_table_timeupdated) #Kruskal-Wallis chi-squared = 13.47, df = 3, p-value = 0.003724 *
#Dunn test
library(FSA)
i <- dunnTest(C.N_ratio~Layer2, data = C.N_table_timeupdated, method="holm")
i
#  Comparison          Z      P.unadj       P.adj
#1      A - B  1.5504164 0.1210415969 0.363124791
#2      A - C  2.3753960 0.0175301309 0.070120524
#3      B - C  1.1113637 0.2664118489 0.532823698
#4      A - D -0.5092986 0.6105429682 0.610542968
#5      B - D -2.3969448 0.0165324148 0.082662074
#6      C - D -3.3167058 0.0009108548 0.005465129 * (1250-1325 vs. 1275-1375)

ggplot(data = C.N_table_timeupdated, aes(Layer2, C.N_ratio, color = Layer2, fill = Layer2)) +
  #geom_point(size = 4, position=position_dodge(width=0.5)) +
  geom_violin(trim= FALSE, alpha=0.6) +
  geom_boxplot(width=0.1, color= "black", fill = "white", alpha=0.8) +
  geom_jitter(position = position_jitter(height = 0, width = .1), size=6.5, colour = "black") +
  #geom_errorbar(data = tgc,
  #              aes(ymin=d13C-ci, ymax=d13C+ci),
  #              width = 0.1, color = "black",
  #              #linetype = "dotted",
  #              position=position_dodge(width=0.5)) +
  #geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


