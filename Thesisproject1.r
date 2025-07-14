#loaded the datasets
library(DESeq2)
#CIBORSORT DATASET
test = read.table(file='/N/project/mmge_pah/Analysis/PAHB_RNASeq/CIBERSORTx_Job5_Adjusted.txt' , header = T, sep = "\t")
ciber = read.table(file="/N/project/mmge_pah/Analysis/PAHB_RNASeq/CIBERSORTx_Job5_Adjusted.txt",header=T,sep="\t")
library(survival)
library(biomaRt)
install.package("gtools")
library(gtools)
library(dplyr) 
library(tidyverse)
install.packages("survminer")
library(survminer)
install.packages("corrplot")
library(corrplot)

#Phenotype Dataset (pah)
load(file="/N/project/mmge_pah/Analysis/PAHB_RNASeq/clinical_treatment_phenotype_plus_omics_data_availability_07_08_2020.Rda")
#RNA- Seq dataset 
RDSfile = readRDS("/N/project/mmge_pah/Analysis/PAHB_RNASeq/BGIvALL_GENE.Filtered.NormalizedTMM.RDS") #RNA- Seq dataset 
colnames(pah)
#loading dataset with EV's 
EV = read.table(file='/N/project/mmge_pah/Analysis/analyses/new_survival_analysis/uchl1_phenotype_covars_v8_plus_key_metabolites_v2.txt', head= T, sep="\t")
colnames(EV)
#extracting the file (counts data) 

gene.counts= counts(RDSfile)
dim(gene.counts) #number of rows and columns
dim(pah) 

ciber$pheno.id=gsub("_","-",ciber$Mixture) #adds an extra column to ciber named pheno.i
key1=match(EV$phenoID,ciber$pheno.id) #matches up IDs (phynotypes/genes to persons covariants)
key1[1:30]
#matches the Neutrophil cells to proper patients in the phenotype dataset
EV$cell.frac=ciber$Neutrophils[key1] 

# rename genes from ENSBL to gene IDs 

# Gene name conversion
ensembl.names=rownames(RDSfile) # remove decimal version information
# Assigning what data to use 
mart=useMart("ENSEMBL_MART_ENSEMBL")
mart=useDataset("hsapiens_gene_ensembl",mart)
# Annotating ENSG ids to gene names
annotLookup=getBM(mart=mart,attributes=c("ensembl_transcript_id", "ensembl_gene_id", "gene_biotype", "external_gene_name"),filter="ensembl_gene_id",values=ensembl.names,uniqueRows=T)

# Replace those ENSG that mapped to a gene name
key=match(ensembl.names,annotLookup$ensembl_gene_id)
# Gene names have duplicated genes (probably pseudogenes or alternatively spliced products as well as non-matched genes
# These will stay as Ensembl ID
gene.names=annotLookup$external_gene_name[key]
gene.names[which(gene.names %in% c("",NA))]=NA
rownames(RDSfile)[which(!is.na(gene.names))]=gene.names[which(!is.na(gene.names))]
#transposing /extracting  
gene.values=t(counts(RDSfile, normalize=T))


#converted the large matrix into a data.frame. 
test3=data.frame(gene.values)






#putting the patient ID in the right format
test3$id.interim=unlist(lapply(rownames(test3),function(x) unlist(strsplit(x,"[.]"))[1]))
test3$id.interim2=gsub("/N/project/mmge_pah/BGI/BGI_data_download_08102020/salmon/quants/quants/","",test3$id.interim)
test3$id.interim3=gsub("_","-",test3$id.interim2)

#merged datasets!! #delete test.merged (the ORIGNIAL LARGE DATASET)
ev.test=merge(EV,test3,by.y="id.interim3",by.x="phenoID")
View(ev.test)


#matched phenoID with their correlating RNAseq data with patient ID  
#key2=match(EV$phenoID,test3$id.interim3)
#key2[1:50]




#Creating a Function (Model #1) 
enroll.run("ev.test","TSPAN6","sex,age.diagnosis,EV1,EV2,EV3,EV4,EV5,cell.frac"),data=pv[which(pah$ancestry=="EA"),])
enroll.run<-function(EV,gene,covars=NA){
  #complete.cases filters the dataset to keep only rows with complete data for the specified gene and covariates 
  test.set=get(EV)[complete.cases(get(EV)[,c(gene,unlist(strsplit(covars,",")))]),] # Getting only complete records for testing
  #test.set=get(EV)[complete.cases(get(EV)[,c(gene,unlist(strsplit(covars,",")))]),]

  test.set2=test.set
  #test.set2 is a compy of the filtered dataset, where the gene expression values are log-trnasformed using 'log10'.
  test.set2[,gene]=log10(test.set2[,gene]+1)
  test.set[,gene]=quantcut(test.set[,gene],3)
  levels(test.set[,gene])=c("Low","Mid","High")
  
  
  #'cox.full' fits the Cox model using the binned gene expression levels and covariates.
  cox.full=coxph(as.formula(paste("Surv(primary.yr,primary.status)~",gene,"+",gsub(",","+",covars),sep="")),data=test.set)
  #'cox.summ' summarizes the results of the full Cox model.
  cox.summ=summary(cox.full)
  #cox.base fits a base Cox model using only the covariates (no gene expression).
  cox.base=coxph(as.formula(paste("Surv(primary.yr,primary.status)~",gsub(",","+",covars),sep="")),data=test.set)
  #cox.phv checks the proportional hazards assumption of the full model using cox.zph.
  cox.phv=cox.zph(cox.full)
  #cox.cont fits a continuous Cox model using the log-transformed gene expression values.
  cox.cont=coxph(as.formula(paste("Surv(primary.yr,primary.status)~",gene,"+",gsub(",","+",covars),sep="")),data=test.set2)
  #cox.sumt summarizes the results of the continuous Cox model.
  cox.sumt=summary(cox.cont)
  
  #storing the results
  output=data.frame(matrix(nrow=1,ncol=19))
  output[1,1]=gene #gene name
  output[1,2]=cox.summ$n #number of observations in the model 
  output[1,3]=cox.summ$nevent # Number of events (e.g. deaths) in the model  
  output[1,4]=cox.summ$coefficients[paste(gene,"Mid",sep=""),1] 
  output[1,5]=cox.summ$coefficients[paste(gene,"Mid",sep=""),3]
  output[1,6]=cox.summ$coefficients[paste(gene,"Mid",sep=""),4]
  output[1,7]=cox.summ$coefficients[paste(gene,"Mid",sep=""),5]
  output[1,8]=cox.summ$coefficients[paste(gene,"High",sep=""),1] 
  output[1,9]=cox.summ$coefficients[paste(gene,"High",sep=""),3]
  output[1,10]=cox.summ$coefficients[paste(gene,"High",sep=""),4]
  output[1,11]=cox.summ$coefficients[paste(gene,"High",sep=""),5]
  output[1,12]=cox.phv$table[gene,3] #p-value from the proportional hazards assumption test
  output[1,13]=summary(cox.base)$rsq[1] #R-squared value for the base model 
  output[1,14]=cox.summ$rsq[1] #R-squared value for the full model 
  output[1,15]=anova(cox.full,cox.base)$P[2] #is the sub model better than the full model? p-value from the likelihood ratio test comparing the full model to the base model
  output[1,16]=cox.sumt$coefficients[gene,1] #continuous coefficient 
  output[1,17]=cox.sumt$coefficients[gene,3] #Standard Deviation 
  output[1,18]=cox.sumt$coefficients[gene,4] #z-statistic
  output[1,19]=cox.sumt$coefficients[gene,5] #p-value
  
  colnames(output)=c("Gene","N","Event","Mid.Coef","Mid.SD","Mid.Z.stat","Mid.pvalue",
                                     "High.Coef","High.SD","High.Z.stat","High.pvalue",
                                     "ZPH.pvalue","base.r2","full.r2","loglike.pvalue",
                                     "Cont.Coef","Cont.SD","Cont.Z.stat","Cont.pvalue")
  return(output)
}

#Generating Results for Multiple Genes 
#A data frame 'output.enroll.survival' is initialized to store the results for each gene.
output.enroll.survival=data.frame(matrix(nrow=length(ev.test),ncol=19)) # subtract the number of added columns

colnames(output.enroll.survival)=c("Gene","N","Event","Mid.Coef","Mid.SD","Mid.Z.stat","Mid.pvalue",
                                   "High.Coef","High.SD","High.Z.stat","High.pvalue",
                                   "ZPH.pvalue","base.r2","full.r2","loglike.pvalue",
                                   "Cont.Coef","Cont.SD","Cont.Z.stat","Cont.pvalue")


#how to delete a column
ev.test.2=ev.test[,-c(23586)]


#give you all the genes all ancestrys: 
output=enroll.run("ev.test",colnames(ev.test)[124],"sex,age.diagnosis,EV1,EV2,EV3,EV4,EV5,cell.frac")
for(i in c(125:23585,23587:(dim(ev.test)[2]-2))){
 output=rbind(output,enroll.run("ev.test",colnames(ev.test)[i],"sex,age.diagnosis,EV1,EV2,EV3,EV4,EV5,cell.frac"))
}

ev.test.EA=ev.test[which(ev.test$ancestry=="EA"),]
ev.test.AA=ev.test[which(ev.test$ancestry=="AA"),]
ev.test.AMR=ev.test[which(ev.test$ancestry=="AMR"),]

write.csv(ev.test.EA, "EURlarge_file.csv", row.names=TRUE) 
output22 <- read.csv("EURlarge_file.csv")

write.csv(ev.test.AA, "AFRLarge_file.csv", row.names=TRUE)
output23 <- read.csv("AFRLarge_file.csv")

write.csv(ev.test.AMR, "AMRLarge_file.csv", row.names=TRUE)
output24 <- read.csv("AMRLarge_file.csv")

#load a CSV file: 
output1 <- read.csv("output1_file.csv")
#delete column one
output1 <- output1[,-1] 

#Running all genes for only the European Americans: 

output1=enroll.run("ev.test.EA",colnames(ev.test.EA)[124],"sex,age.diagnosis,EV1,EV2,EV3,EV3,EV4,EV5,cell.frac")
for(i in c(125:10346,10348:16232,16234:16974,16976:19438,19440:19965,19967:21057,21059:21391,21393:21667,21669:21753,21755:
           22248,22250:22321,22323:22423,22425:22841,22843:23585,23587:23979,23981:24518,24520:24844,24846:25120,25122:25822,25824:27017,27019:
           28672,28674:29033,29035:29578,29580:31302,31304:31582,31584:32159,32161:32623,32625:33162,33164:37031,37033:37881,37883:38768,38770:39261,39263:
           41236,41238:(dim(ev.test.EA)[2]-2))){ 
    output1=rbind(output1,enroll.run("ev.test.EA",colnames(ev.test.EA)[i],"sex,age.diagnosis,EV1,EV2,EV3,EV4,EV5,cell.frac"))
}

#output1=enroll.run("ev.test.EA",colnames(ev.test.EA)[124],"sex,age.diagnosis,EV1,EV2,EV3,EV3,EV4,EV5,cell.frac")
#for(i in c(125::(dim(ev.test.EA)[2]-2))){ 
#  output1=rbind(output1,try(enroll.run("ev.test.EA",colnames(ev.test.EA)[i],"sex,age.diagnosis,EV1,EV2,EV3,EV4,EV5,cell.frac"))
#}

#or this way! 
output1=enroll.run("ev.test.EA",colnames(ev.test.EA)[124],"sex,age.diagnosis,EV1,EV2,EV3,EV3,EV4,EV5,cell.frac")
for(i in 125:46037){#dim(ev.test.AA)[2]){ 
  output1=try(rbind(output1,try(enroll.run("ev.test.EA",colnames(ev.test.EA)[i],"sex,age.diagnosis,EV1,EV2,EV3,EV4,EV5,cell.frac"))))
}

#converting output1 dataframe to a csv
write.csv(output1, "output1_file.csv", row.names=TRUE) 
#load a CSV file: 
output1 <- read.csv("output1_file.csv")
#delete column one
output1 <- output1[,-1] 


#converting ZPH.pvalue and High.pvalue to numeric 
output1t <- output1

#Converted all the rows in to numeric values except for Genes column
output1t <- output1t %>% 
           mutate(across(-Gene, as.numeric))

# got rid on all the NA in High.pvalues 
output1t <- output1t %>% 
           filter(!is.na(output1t$High.pvalue & output1t$ZPH.pvalue))


filtered1 <- output1t %>% 
  filter(High.pvalue < 0.000000336 & ZPH.pvalue>=0.05)




#filtering 
rm(filtere_values)
hist(newXa)
hist(filtered1$ZPH.pvalue)



ggplot(filtered1, aes(x=High.pvalue, na.rm=TRUE))+geom_histogram(binwidth = 200)

#Running all genes for only the African Americans: 

output2=enroll.run("ev.test.AA",colnames(ev.test.AA)[124],"sex,age.diagnosis,EV1,EV2,EV3,EV3,EV4,EV5,cell.frac")
for(i in 125:46037){#dim(ev.test.AA)[2]){ 
  output2=try(rbind(output2,try(enroll.run("ev.test.AA",colnames(ev.test.AA)[i],"sex,age.diagnosis,EV1,EV2,EV3,EV4,EV5,cell.frac"))))
}
write.csv(output2, "output2_file.csv", row.names=TRUE) 

#converting ZPH.pvalue and High.pvalue to numeric 
output2t <- output2

#Converted all the rows in to numeric values except for Genes column
output2t <- output2t %>% 
  mutate(across(-Gene, as.numeric))

# got rid on all the NA in High.pvalues 
output2t <- output2t %>% 
  filter(!is.na(output2t$High.pvalue & output2t$ZPH.pvalue))


filtered2 <- output2t %>% 
  filter(High.pvalue < 0.05 & ZPH.pvalue>=0.05)

#Running all genes for only the Hispanics: 

output3=enroll.run("ev.test.AMR",colnames(ev.test.AMR)[124],"sex,age.diagnosis,EV1,EV2,EV3,EV3,EV4,EV5,cell.frac")
for(i in 125:46037){#dim(ev.test.AMR)[2]){ 
  output3=try(rbind(output3,try(enroll.run("ev.test.AMR",colnames(ev.test.AMR)[i],"sex,age.diagnosis,EV1,EV2,EV3,EV4,EV5,cell.frac"))))
}

#load a CSV file: 
output3 <- read.csv("output3_file.csv")
#delete column one
output3 <- output3[,-1] 

#converting ZPH.pvalue and High.pvalue to numeric 
output3t <- output3

#Converted all the rows in to numeric values except for Genes column
output3t <- output3t %>% 
  mutate(across(-Gene, as.numeric))

# got rid on all the NA in High.pvalues 
output3t <- output3t %>% 
  filter(!is.na(output3t$High.pvalue & output3t$ZPH.pvalue))


filtered3 <- output3t %>% 
  filter(High.pvalue < 0.05 & ZPH.pvalue>=0.05)

#Second Model (including covariants (PAH and PAH)

output4=enroll.run("ev.test.EA",colnames(ev.test.EA)[124],"sex,age.diagnosis,PAH.x,PVR.x,EV1,EV2,EV3,EV3,EV4,EV5,cell.frac")
for(i in 125:46037){#dim(ev.test.EA)[2]){ 
  output4=try(rbind(output4,try(enroll.run("ev.test.EA",colnames(ev.test.EA)[i],"sex,age.diagnosis,PAH.x,PVR.x,EV1,EV2,EV3,EV4,EV5,cell.frac"))))
}

#converting output4 dataframe to a csv
write.csv(output4, "output4_file.csv", row.names=TRUE) 
#load a CSV file: 
output4 <- read.csv("output4_file.csv")
#delete column one
output4 <- output4[,-1] 


#converting ZPH.pvalue and High.pvalue to numeric 
output4t <- output4

#Converted all the rows in to numeric values except for Genes column
output4t <- output4t %>% 
  mutate(across(-Gene, as.numeric))

# got rid on all the NA in High.pvalues 
output4t <- output4t %>% 
  filter(!is.na(output4t$High.pvalue & output4t$ZPH.pvalue))


filtered4 <- output4t %>% 
  filter(High.pvalue < 0.00000186 & ZPH.pvalue >=0.05)

write.csv(filtered4, "filtered4_file.csv", row.names=TRUE) 
output4 <- as.data.frame(read.csv("filtered4_file.csv"))

filtered4.exp<-exp(output4$High.Coef)
filtered4.exp<- cbind(output4$Gene, output4$N,output4$High.SD, filtered4.exp) 
colnames(filtered4.exp)<-c("Gene","N","SD", "HR")
SEeu <- output4$High.SD/sqrt(output4$N)
filtered4.exp = as.data.frame(cbind(filtered4.exp, SEeu))
filtered4.exp$HR = as.numeric(filtered4.exp$HR)
filtered4.exp$SEeu = as.numeric(filtered4.exp$SEeu)
filtered4.exp$lowcieu = exp(log(filtered4.exp$HR)-1.96*filtered4.exp$SEeu)
filtered4.exp$upcieu = exp(log(filtered4.exp$HR)+1.96*filtered4.exp$SEeu) 
filtered4.exp$pvalue= output4$High.pvalue

write.csv(filtered4.exp, "EUsecondtable.csv", row.names=TRUE) 

  
#Kapalan Miere Survival graphs for the 11 genes for EU! 
#plot(survfit(Surv(primary.yr,primary.status)~quantcut(REXO2,3),data=ev.test.EA))

tiff(file="EUREXO2.tif", width=8, height =6, unit="in",res=600, compression="lzw")
EUREXO2 <- survfit(Surv(primary.yr,primary.status)~quantcut(REXO2,3),data=ev.test.EA)
ggsurvplot(EUREXO2, data =ev.test.EA,title="A1", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="EUTMC5.tif", width=8, height =6, unit="in",res=600, compression="lzw")
EUTMC5 <- survfit(Surv(primary.yr,primary.status)~quantcut(TMC5,3),data=ev.test.EA)
ggsurvplot(EUTMC5, data =ev.test.EA, risk.table = TRUE,title="KM for EUTMC5", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="EUCPE.tif", width=8, height =6, unit="in",res=600, compression="lzw")
EUCPE <- survfit(Surv(primary.yr,primary.status)~quantcut(CPE,3),data=ev.test.EA)
ggsurvplot(EUCPE, data =ev.test.EA, risk.table = TRUE,title="KM for EUCPE", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="EUFHL2.tif", width=8, height =6, unit="in",res=600, compression="lzw")
EUFHL2 <- survfit(Surv(primary.yr,primary.status)~quantcut(FHL2,3),data=ev.test.EA)
ggsurvplot(EUFHL2, data =ev.test.EA, risk.table = FALSE,title="B1", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="EURAB13.tif", width=8, height =6, unit="in",res=600, compression="lzw")
EURAB13 <- survfit(Surv(primary.yr,primary.status)~quantcut(RAB13,3),data=ev.test.EA)
ggsurvplot(EURAB13, data =ev.test.EA, risk.table = TRUE,title="KM for EURAB13", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="EUH4C8.tif", width=8, height =6, unit="in",res=600, compression="lzw")
EUH4C8 <- survfit(Surv(primary.yr,primary.status)~quantcut(H4C8,3),data=ev.test.EA)
ggsurvplot(EUH4C8, data =ev.test.EA, risk.table = TRUE,title="KM for EUH4C8", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="EUTHOC7.tif", width=8, height =6, unit="in",res=600, compression="lzw")
EUTHOC7 <- survfit(Surv(primary.yr,primary.status)~quantcut(THOC7,3),data=ev.test.EA)
ggsurvplot(EUTHOC7, data =ev.test.EA, risk.table = TRUE,title="KM for EUTHOC7", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="EUCABP4.tif", width=8, height =6, unit="in",res=600, compression="lzw")
EUCABP4 <- survfit(Surv(primary.yr,primary.status)~quantcut(CABP4,3),data=ev.test.EA)
ggsurvplot(EUCABP4, data =ev.test.EA, risk.table = FALSE,title="C1", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="EUZBTB7C.tif", width=8, height =6, unit="in",res=600, compression="lzw")
EUZBTB7C <- survfit(Surv(primary.yr,primary.status)~quantcut(ZBTB7C,3),data=ev.test.EA)
ggsurvplot(EUZBTB7C, data =ev.test.EA, risk.table = TRUE,title="KM for EUZBTB7C", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="EUH12.tif", width=8, height =6, unit="in",res=600, compression="lzw")
EUH12 <- survfit(Surv(primary.yr,primary.status)~quantcut(H1.2,3),data=ev.test.EA)
ggsurvplot(EUH12, data =ev.test.EA, risk.table = TRUE,title="KM for EUH1.2", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="EULINC02345.tif", width=8, height =6, unit="in",res=600, compression="lzw")
EULINC02345 <- survfit(Surv(primary.yr,primary.status)~quantcut(LINC02345,3),data=ev.test.EA)
ggsurvplot(EULINC02345, data =ev.test.EA, risk.table = TRUE,title="KM for LINC02345", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

dev.off()

#Boxplots: 
library(ggplot2)
ev.testbox <- ev.test
ev.testbox$primary.status <- factor(ev.testbox$primary.status)

ggplot(ev.testbox, aes(x=ancestry, y=REXO2, fill=sex))+
  geom_boxplot()

ggplot(ev.testbox, aes(x=ancestry, y=REXO2, fill=primary.status))+
  geom_boxplot()+
  labs(title="REX02 gene", x="Ancestry", y=" REX02 Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=REXO2, fill= interaction(primary.status, sex)))+
  geom_boxplot()+
  labs(fill="Primary Status & Sex", title="REX02 Gene", x="Ancestry", y="Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=TMC5, fill=primary.status))+
  geom_boxplot()+
  labs(title="TMC5 gene", x="Ancestry", y=" TMC5 Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=TMC5, fill= interaction(primary.status, sex)))+
  geom_boxplot()+
  labs(fill="Primary Status & Sex", title="TMC5 Gene", x="Ancestry", y=" TMC5 Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=CPE, fill=primary.status))+
  geom_boxplot()+
  labs(title="CPE gene", x="Ancestry", y="CPE Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=CPE, fill= interaction(primary.status, sex)))+
  geom_boxplot()+
  labs(fill="Primary Status & Sex", title="CPE Gene", x="Ancestry", y=" CPE Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=FHL2, fill=primary.status))+
  geom_boxplot()+
  labs(title="FHL2 gene", x="Ancestry", y="FHL2 Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=FHL2, fill= interaction(primary.status, sex)))+
  geom_boxplot()+
  labs(fill="Primary Status & Sex", title="FHL2 Gene", x="Ancestry", y="FHL2 Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=RAB13, fill=primary.status))+
  geom_boxplot()+
  labs(title="RAB13 gene", x="Ancestry", y="RAB13 Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=RAB13, fill= interaction(primary.status, sex)))+
  geom_boxplot()+
  labs(fill="Primary Status & Sex", title="RAB13 Gene", x="Ancestry", y="RAB13 Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=H4C8, fill=primary.status))+
  geom_boxplot()+
  labs(title="H4C8 gene", x="Ancestry", y="H4C8 Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=H4C8 , fill= interaction(primary.status, sex)))+
  geom_boxplot()+
  labs(fill="Primary Status & Sex", title="H4C8  Gene", x="Ancestry", y="H4C8  Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=THOC7, fill=primary.status))+
  geom_boxplot()+
  labs(title="THOC7 gene", x="Ancestry", y="THOC7 Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=THOC7, fill= interaction(primary.status, sex)))+
  geom_boxplot()+
  labs(fill="Primary Status & Sex", title="THOC7 Gene", x="Ancestry", y="THOC7 Gene Expression")

#thesis Boxplot 
#got rid of the NA's 
ev.testbox_clean <- ev.testbox %>%
  filter(!is.na(primary.status))

ggplot(ev.testbox_clean, aes(x=ancestry, y=CABP4, fill=primary.status))+
  geom_boxplot()+
  labs(title="", x="Ancestry", y="CABP4 Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=CABP4, fill= interaction(primary.status, sex)))+
  geom_boxplot()+
  labs(fill="Primary Status & Sex", title="CABP4 Gene", x="Ancestry", y="CABP4 Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=ZBTB7C, fill=primary.status))+
  geom_boxplot()+
  labs(title="ZBTB7C gene", x="Ancestry", y="ZBTB7C Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=ZBTB7C, fill= interaction(primary.status, sex)))+
  geom_boxplot()+
  labs(fill="Primary Status & Sex", title="ZBTB7C Gene", x="Ancestry", y="ZBTB7C Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=H1.2, fill=primary.status))+
  geom_boxplot()+
  labs(title="H1.2 gene", x="Ancestry", y="H1.2 Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=H1.2, fill= interaction(primary.status, sex)))+
  geom_boxplot()+
  labs(fill="Primary Status & Sex", title="H1.2 Gene", x="Ancestry", y="H1.2 Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=LINC02345, fill=primary.status))+
  geom_boxplot()+
  labs(title="LINC02345 gene", x="Ancestry", y="LINC02345 Gene Expression")

ggplot(ev.testbox, aes(x=ancestry, y=LINC02345, fill= interaction(primary.status, sex)))+
  geom_boxplot()+
  labs(fill="Primary Status & Sex", title="LINC02345 Gene", x="Ancestry", y="LINC02345 Gene Expression")

#Kapalan Miere Survival graphs for the 11 genes for AA! 
tiff(file="AAREXO2.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AAREXO2 <- survfit(Surv(primary.yr,primary.status)~quantcut(REXO2,3),data=ev.test.AA)
ggsurvplot(AAREXO2, data =ev.test.AA, risk.table = FALSE,title="A2", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AATMC5.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AATMC5 <- survfit(Surv(primary.yr,primary.status)~quantcut(TMC5,3),data=ev.test.AA)
ggsurvplot(AATMC5, data =ev.test.AA, risk.table = TRUE,title="KM for AATMC5", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AACPE.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AACPE <- survfit(Surv(primary.yr,primary.status)~quantcut(CPE,3),data=ev.test.AA)
ggsurvplot(AACPE, data =ev.test.AA, risk.table = TRUE,title="KM for AACPE")
dev.off()

tiff(file="AAFHL2.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AAFHL2 <- survfit(Surv(primary.yr,primary.status)~quantcut(FHL2,3),data=ev.test.AA)
ggsurvplot(AAFHL2, data =ev.test.AA, risk.table = FALSE,title="B2", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AARAB13.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AARAB13 <- survfit(Surv(primary.yr,primary.status)~quantcut(RAB13,3),data=ev.test.AA)
ggsurvplot(AARAB13, data =ev.test.AA, risk.table = TRUE,title="KM for AARAB13", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AAH4C8.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AAH4C8 <- survfit(Surv(primary.yr,primary.status)~quantcut(H4C8,3),data=ev.test.AA)
ggsurvplot(AAH4C8, data =ev.test.AA, risk.table = TRUE,title="KM for AAH4C8", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AATHOC7.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AATHOC7 <- survfit(Surv(primary.yr,primary.status)~quantcut(THOC7,3),data=ev.test.AA)
ggsurvplot(AATHOC7, data =ev.test.AA, risk.table = TRUE,title="KM for AATHOC7", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AACABP4.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AACABP4 <- survfit(Surv(primary.yr,primary.status)~quantcut(CABP4,3),data=ev.test.AA)
ggsurvplot(AACABP4, data =ev.test.AA, risk.table = FALSE,title="C2", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AAZBTB7C.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AAZBTB7C <- survfit(Surv(primary.yr,primary.status)~quantcut(ZBTB7C,3),data=ev.test.AA)
ggsurvplot(AAZBTB7C, data =ev.test.AA, risk.table = TRUE,title="KM for AAZBTB7C", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AAH12.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AAH12 <- survfit(Surv(primary.yr,primary.status)~quantcut(H1.2,3),data=ev.test.AA)
ggsurvplot(AAH12, data =ev.test.AA, risk.table = TRUE,title="KM for AAH1.2", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AALINC02345.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AALINC02345 <- survfit(Surv(primary.yr,primary.status)~quantcut(LINC02345,3),data=ev.test.AA)
ggsurvplot(AALINC02345, data =ev.test.AA, risk.table = TRUE,title="KM for AALINC02345", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()


#Kapalan Miere Survival graphs for the 11 genes for AMR! 

tiff(file="AMRREXO2.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AMRREXO2 <- survfit(Surv(primary.yr,primary.status)~quantcut(REXO2,3),data=ev.test.AMR)
ggsurvplot(AMRREXO2, data =ev.test.AMR, risk.table = FALSE,title="A3", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AMRTMC5.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AMRTMC5 <- survfit(Surv(primary.yr,primary.status)~quantcut(TMC5,3),data=ev.test.AMR)
ggsurvplot(AMRTMC5, data =ev.test.AMR, risk.table = TRUE,title="KM for AMRTMC5", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AMRCPE.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AMRCPE <- survfit(Surv(primary.yr,primary.status)~quantcut(CPE,4),data=ev.test.AMR)
ggsurvplot(AMRCPE, data =ev.test.AMR, risk.table = TRUE,title="KM for AMRCPE")
dev.off()

tiff(file="AMRFHL2.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AMRFHL2 <- survfit(Surv(primary.yr,primary.status)~quantcut(FHL2,3),data=ev.test.AMR)
ggsurvplot(AMRFHL2, data =ev.test.AMR, risk.table = FALSE,title="B3", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AMRRAB13.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AMRRAB13 <- survfit(Surv(primary.yr,primary.status)~quantcut(RAB13,3),data=ev.test.AMR)
ggsurvplot(AMRRAB13, data =ev.test.AMR, risk.table = TRUE,title="KM for AMRRAB13", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AMRH4C8.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AMRH4C8 <- survfit(Surv(primary.yr,primary.status)~quantcut(H4C8,3),data=ev.test.AMR)
ggsurvplot(AMRH4C8, data =ev.test.AMR, risk.table = TRUE,title="KM for AMRH4C8", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AMRTHOC7.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AMRTHOC7 <- survfit(Surv(primary.yr,primary.status)~quantcut(THOC7,3),data=ev.test.AMR)
ggsurvplot(AMRTHOC7, data =ev.test.AMR, risk.table = TRUE,title="KM for AMRTHOC7", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AMRCABP4.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AMRCABP4 <- survfit(Surv(primary.yr,primary.status)~quantcut(CABP4,3),data=ev.test.AMR)
ggsurvplot(AMRCABP4, data =ev.test.AMR, risk.table = FALSE,title="C3", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AMRZBTB7C.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AMRZBTB7C <- survfit(Surv(primary.yr,primary.status)~quantcut(ZBTB7C,3),data=ev.test.AMR)
ggsurvplot(AMRZBTB7C, data =ev.test.AMR, risk.table = TRUE,title="KM for AMRZBTB7C", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AMRH12.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AMRH12 <- survfit(Surv(primary.yr,primary.status)~quantcut(H1.2,3),data=ev.test.AMR)
ggsurvplot(AMRH12, data =ev.test.AMR, risk.table = TRUE,title="KM for AMRH1.2", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()

tiff(file="AMRLINC02345.tif", width=8, height =6, unit="in",res=600, compression="lzw")
AMRLINC02345 <- survfit(Surv(primary.yr,primary.status)~quantcut(LINC02345,3),data=ev.test.AMR)
ggsurvplot(AMRLINC02345, data =ev.test.AMR, risk.table = TRUE,title="KM for AMRLINC02345", legend.labs = c("Bottom Tertile", "Middle Tertile", "Top Tertile"))
dev.off()



output5=enroll.run("ev.test.AA",colnames(ev.test.AA)[124],"sex,age.diagnosis,PVR.x,EV1,EV2,EV3,EV3,EV4,EV5,cell.frac")
for(i in 125:46037){#dim(ev.test.AA)[2]){ 
  output5=try(rbind(output5,try(enroll.run("ev.test.AA",colnames(ev.test.AA)[i],"sex,age.diagnosis,PVR.x,EV1,EV2,EV3,EV4,EV5,cell.frac"))))
}

#converting output5 dataframe to a csv
write.csv(output5, "output5_file.csv", row.names=TRUE) 

#load a CSV file: 
output5 <- read.csv("output5_file.csv")
#delete column one
output5 <- output5[,-1] 

#converting ZPH.pvalue and High.pvalue to numeric 
output5t <- output5

#Converted all the rows in to numeric values except for Genes column
output5t <- output5t %>% 
  mutate(across(-Gene, as.numeric))

# got rid on all the NA in High.pvalues 
output5t <- output5t %>% 
  filter(!is.na(output5t$High.pvalue & output5t$ZPH.pvalue))


filtered5 <- output5t %>% 
  filter(ZPH.pvalue >=.05 & High.pvalue < 0.05)



filtered_df <- output5t[output5t$Gene %in% c("REXO2", "TMC5", "CPE", "FHL2", "RAB13", "H4C8", "THOC7", "CABP4","ZBTB7C", "H1.2","LINC02345"), ]
write.csv(filtered_df, "filteredAFR11_file.csv", row.names=TRUE) 
output5 <- read.csv("filteredAFR11_file.csv")

#table 2 for poster!!!
table2AFR = exp(output5$High.Coef)
table2AFR = cbind(output5$Gene, output5$N,output5$High.SD, table2AFR) 
colnames(table2AFR) = c("Gene","N","SD", "HR")
SEAA = output5$High.SD/sqrt(output5$N)
table2AFR = as.data.frame(cbind(table2AFR, SEAA))
table2AFR$HR = as.numeric(table2AFR$HR)
table2AFR$SEAA = as.numeric(table2AFR$SEAA)
table2AFR$lowcieu = exp(log(table2AFR$HR)-1.96*table2AFR$SEAA)
table2AFR$upcieu = exp(log(table2AFR$HR)+1.96*table2AFR$SEAA) 
table2AFR$pvalue= output5$High.pvalue
write.csv(table2AFR, "AFRsecondtable.csv", row.names=TRUE)
  
  
  
output6=enroll.run("ev.test.AMR",colnames(ev.test.AMR)[124],"sex,age.diagnosis,PVR.x,EV1,EV2,EV3,EV3,EV4,EV5,cell.frac")
for(i in 125:46037){#dim(ev.test.AMR)[2]){ 
  output6=try(rbind(output6,try(enroll.run("ev.test.AMR",colnames(ev.test.AMR)[i],"sex,age.diagnosis,PVR.x,EV1,EV2,EV3,EV4,EV5,cell.frac"))))
}

#converting output6 dataframe to a csv
write.csv(output6, "output6_file.csv", row.names=TRUE) 

#load a CSV file: 
output6 <- read.csv("output6_file.csv")
#delete column one
output6 <- output6[,-1] 

#converting ZPH.pvalue and High.pvalue to numeric 
output6t <- output6

#Converted all the rows in to numeric values except for Genes column
output6t <- output6t %>% 
  mutate(across(-Gene, as.numeric))

# got rid on all the NA in High.pvalues 
output6t <- output6t %>% 
  filter(!is.na(output6t$High.pvalue & output6t$ZPH.pvalue))

filtered6 <- output6t %>% 
  filter(ZPH.pvalue >=.05 & High.pvalue < 0.05)

filtered_11 <- output6t[output6t$Gene %in% c("REXO2", "TMC5", "CPE", "FHL2", "RAB13", "H4C8", "THOC7", "CABP4","ZBTB7C", "H1.2","LINC02345"), ]
write.csv(filtered_11, "filteredAMR11_file.csv", row.names=TRUE) 

output6 <- read.csv("filteredAMR11_file.csv")

#table 2 for poster!!!
table2AMR = exp(output6$High.Coef)
table2AMR = cbind(output6$Gene, output6$N,output6$High.SD, table2AMR) 
colnames(table2AMR) = c("Gene","N","SD", "HR")
SE = output6$High.SD/sqrt(output6$N)
table2AMR = as.data.frame(cbind(table2AMR, SE))
table2AMR$HR = as.numeric(table2AMR$HR)
table2AMR$SE = as.numeric(table2AMR$SE)
table2AMR$lowcieu = exp(log(table2AMR$HR)-1.96*table2AMR$SE)
table2AMR$upcieu = exp(log(table2AMR$HR)+1.96*table2AMR$SE) 
table2AMR$pvalue= output6$High.pvalue
write.csv(table2AMR, "AMRsecondtable.csv", row.names=TRUE)





library(ggplot2)
library(MASS)
library(dplyr)

#correlation between the 11 genes 
ev.test.EV11 <- ev.test.EA |>
       select(REXO2, TMC5, CPE, FHL2, RAB13, H4C8, THOC7, CABP4,ZBTB7C,H1.2,LINC02345)
EAcor <- cor(ev.test.EV11)
corrplot(EAcor, method="circle")
          title("Correlation for the EUR 11 genes")
          write.csv(EAcor, "EUcorrelation_file.csv", row.names=TRUE)

#correlations between the 8 Genes        
ev.test.EV8 <- ev.test.EA[, c("REXO2", "TMC5", "CPE", "FHL2", "H4C8", "THOC7", "CABP4","ZBTB7C")]    #selecting specific columns   
EAcor <- cor(ev.test.EV8)
corrplot(EAcor, method="circle")
title("EUR")          
          
ev.test.AA11 <- ev.test.AA |>
        select(REXO2, TMC5, CPE, FHL2, RAB13, H4C8, THOC7, CABP4,ZBTB7C,H1.2,LINC02345)
AAcor <- cor(ev.test.AA11)
corrplot(AAcor, method="circle")
             title("Correlation for the AFR 11 genes")
             write.csv(AAcor, "AFRcorrelation_file.csv", row.names=TRUE)

#correlations between the 8 Genes        
ev.test.AA8 <- ev.test.AA[, c("REXO2", "TMC5", "CPE", "FHL2", "H4C8", "THOC7", "CABP4","ZBTB7C")]    #selecting specific columns   
AAcor <- cor(ev.test.AA8)
corrplot(AAcor, method="circle")
title("AFR")             
             
ev.test.AMR11 <- ev.test.AMR |>
  select(REXO2, TMC5, CPE, FHL2, RAB13, H4C8, THOC7, CABP4,ZBTB7C,H1.2,LINC02345)
AMRcor <- cor(ev.test.AMR11)
corrplot(AMRcor, method="circle")
title("Correlation for the AMR 11 genes")
write.csv(AMRcor, "AMRcorrelation_file.csv", row.names=TRUE)

ev.test.AMR8 <- ev.test.AMR[, c("REXO2", "TMC5", "CPE", "FHL2", "H4C8", "THOC7", "CABP4","ZBTB7C")]    #selecting specific columns   
AMRcor <- cor(ev.test.AMR8)
corrplot(AMRcor, method="circle")
title("AMR")       

#VennDiagram
install.packages("ggVennDiagram")
library(ggVennDiagram)

list_venn <- list(
  EUR1=filtered4$Gene,
  AFR2=filtered5$Gene,
  AMR3=filtered6$Gene
)

ggVennDiagram(list_venn) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")






#practice 3 

# Function to find shared values and return specific columns from the data frames
find_shared_rows_specific_columns <- function(filtered4, filtered5, filtered6, column_name, cols_df1, cols_df2, cols_df3) {
  
  # Extract the values from the specified column in each data frame
  values_df1 <- filtered4[[column_name]]
  values_df2 <- filtered5[[column_name]]
  values_df3 <- filtered6[[column_name]]
  
  # Find the intersection of the values across the three data frames
  shared_values <- Reduce(intersect, list(values_df1, values_df2, values_df3))
  
  # Filter the rows in each data frame where the shared values are found
  df1_filtered <- filtered4[filtered4[[column_name]] %in% shared_values, c(column_name, cols_df1)]
  df2_filtered <- filtered5[filtered5[[column_name]] %in% shared_values, c(column_name, cols_df2)]
  df3_filtered <- filtered6[filtered6[[column_name]] %in% shared_values, c(column_name, cols_df3)]
  
  # Merge the filtered data frames based on the common column (e.g., 'Name')
  merged_df <- merge(df1_filtered, df2_filtered, by = column_name, suffixes = c("_EU", "_AA"))
  merged_df <- merge(merged_df, df3_filtered, by = column_name, suffixes = c("_AMR"))
  
  # Return the merged data frame with only the selected columns
  return(merged_df)
}

# Example data frames


# Use the function to find shared rows based on the "Name" column
# Select specific columns from each data frame (e.g., "Age", "City" from df1, df2, and df3)
allsharedgenes <- find_shared_rows_specific_columns(filtered4, filtered5, filtered6, "Gene",
                                                    cols_df1 = c("High.Coef", "High.pvalue"),
                                                    cols_df2 = c("High.Coef", "High.pvalue"),
                                                    cols_df3 = c("High.Coef", "High.pvalue", "ZPH.pvalue"))

#adding AMR to last three columns. 
colnames(allsharedgenes)[6:8]=paste(colnames(shared_rows_df)[6:8],"",sep="")
colnames(allsharedgenes)

#
write.csv(allsharedgenes, "allsharedgenesfile.csv", row.names=TRUE) 


#all genes that have a negative coefficient 
alldecrease <- allsharedgenes %>% 
  filter(High.Coef_EU <0 & High.Coef_AA <0 & High.Coef_AMR <0)
write.csv(alldecrease, "alldecreasefile.csv", row.names=TRUE) 

#all genes that have a positive coefficient 
allincrease <- allsharedgenes %>% 
  filter(High.Coef_EU >0 & High.Coef_AA >0 & High.Coef_AMR >0)

write.csv(allincrease, "allincreasefile.csv", row.names=TRUE) 

#sorting parts of the Venn Diagram 
EUR=filtered4[,c("Gene","High.pvalue","High.Coef")] 
AFR =filtered5[,c("Gene","High.pvalue","High.Coef")] 
AMR =filtered6[,c("Gene","High.pvalue","High.Coef")] 

colnames(EUR)[2:3]=paste(colnames(EUR)[2:3],"",sep="_EU")
colnames(AFR)[2:3]=paste(colnames(AFR)[2:3],"",sep="_AF")
colnames(AMR)[2:3]=paste(colnames(AMR)[2:3],"",sep="_HIS")

EUR$estatus=TRUE
AFR$astatus=TRUE
AMR$hstatus=TRUE


popall <- merge(EUR, AMR, by="Gene", all=TRUE)
popall2 <- merge(popall, AFR, by="Gene", all=TRUE)

onlyEURAMR <- popall2[which(popall2$estatus==TRUE&popall2$hstatus==TRUE&is.na(popall2$astatus)),]
onlyEURnotAMR <- popall2[which(popall2$estatus==TRUE&is.na(popall2$astatus)),]


filtered_11 <- output6t[output6t$Gene %in% c("REXO2", "TMC5", "CPE", "FHL2", "RAB13", "H4C8", "THOC7", "CABP4","ZBTB7C", "H1.2","LINC02345"), ]
write.csv(filtered_11, "filteredAMR11_file.csv", row.names=TRUE) 





all_zero_columns <- which(colSums(ev.test.EA !=0) ==0)
column_number <- which(colnames(ev.test.EA) == "column_name")
#Single Gene Analysis 
enroll.run("ev.test","TSPAN6","sex,age.diagnosis,EV1,EV2,EV3,EV4,EV5,cell.frac") data=ev.test[which(ev.test$ancestry=="EA"),])



#The gene names are cleaned by replacing periods (".") with hyphens ("-").
output.enroll.survival$Gene=unlist(lapply(output.enroll.survival$Gene,function(x) gsub("[.]","-",x))) 
output.enroll.survival(output)
output(output.enroll.survival)

#Whinney Mann Test 

df1_dead <- ev.test.EA[ev.test.EA$primary.status == 1, ]
df2_dead <- ev.test.AMR[ev.test.AMR$primary.status ==1, ]

results <- c() 
for(gene in c("REXO2","TMC5", "CPE", "FHL2", "RAB13","H4C8", "THOC7", "CABP4", "ZBTB7C", "H1.2", "LINC02345")){
  geneEA <- df1_dead[[gene]]
  geneAMR <- df2_dead[[gene]]
  
  test_result <- wilcox.test(geneEA, geneAMR, alternative ="two.sided")
  print(gene)
  print(test_result)
  results <- c(results, test_result)
}



print(test_result)


#calculating the p-value between two genes that are positively corrrelated
REXO2 <- ev.test$REXO2
FHL2 <- ev.test$FHL2
correlation_result <- cor.test(REXO2, FHL2, method = "pearson")
print(correlation_result)

#calculating the p-value between three genes that are positively corrrelated
REXO2 <- ev.test$REXO2
FHL2 <- ev.test$FHL2
CABP4 <- ev.test$CABP4
correlation_result3 <- cor.test(REXO2, CABP4, method ="pearson")
print(correlation_result3)

correlation_result4 <- cor.test(FHL2, CABP4, method ="pearson")
print(correlation_result4)

# Ran a Partial Correlation
install.packages("ppcor")
library(ppcor)
ev.test.3 <- ev.test[, c("REXO2", "FHL2", "CABP4")] 
# Compute partial correlation
partial_cor <- pcor(ev.test.3)

# View the results
print(partial_cor)


#extracting columns from ev.test for ID, FHL2, REXO2, and CABP4
library(dplyr)
threegenes <- ev.test %>% select(IID,phenoID, ancestry, FHL2, CABP4, REXO2)
print(threegenes)
write.csv(threegenes, "threegenes.csv", row.names = FALSE)
write.table(threegenes, "threegenes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Extracted the 3 genes and stratified them by EUR, AMR, and AFR 
threegenes.EUR=threegenes[which(threegenes$ancestry=="EA"),]
threegenes.AMR=threegenes[which(threegenes$ancestry=="AMR"),]
threegenes.AFR=threegenes[which(threegenes$ancestry=="AA"),]

write.table(threegenes.EUR, "threegenesEUR.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(threegenes.AMR, "threegenesAMR.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(threegenes.AFR, "threegenesAFR.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#log10 for EUR the gene expressions
threegenes.log.EUR = threegenes.EUR
threegenes.log.EUR$FHL2 = log10(threegenes.log.EUR$FHL2)
threegenes.log.EUR$CABP4 = log10(threegenes.log.EUR$CABP4)
threegenes.log.EUR$REXO2 = log10(threegenes.log.EUR$REXO2)

write.table(threegenes.log.EUR, "/N/project/mmge_pah/Analysis/analyses/SOX18/threegenes.log.EUR.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#log10 for AFR the gene expressions
threegenes.log.AFR = threegenes.AFR
threegenes.log.AFR$FHL2 = log10(threegenes.log.AFR$FHL2)
threegenes.log.AFR$CABP4 = log10(threegenes.log.AFR$CABP4)
threegenes.log.AFR$REXO2 = log10(threegenes.log.AFR$REXO2)

write.table(threegenes.log.AFR, "/N/project/mmge_pah/Analysis/analyses/SOX18/threegenes.log.AFR.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#log10 for AMR the gene expressions
threegenes.log.AMR = threegenes.AMR
threegenes.log.AMR$FHL2 = log10(threegenes.log.AMR$FHL2)
threegenes.log.AMR$CABP4 = log10(threegenes.log.AMR$CABP4)
threegenes.log.AMR$REXO2 = log10(threegenes.log.AMR$REXO2)

write.table(threegenes.log.AMR, "/N/project/mmge_pah/Analysis/analyses/SOX18/threegenes.log.AMR.txt", sep = "\t", row.names = FALSE, quote = FALSE)


### Callie ruining your code
metabolite = read.table(file='/N/project/mmge_pah/Analysis/analyses/SOX18/uchl1_phenotype_covars_v8_plus_key_metabolites_v2.txt' , header = T, sep = "\t")
keyx=match(metabolite$phenoID, ciber$pheno.id)
metabolite$Neutrophils<-ciber$Neutrophils[keyx]
temp=metabolite[which(!is.na(metabolite$sex)),]
write.table(temp, "/N/project/mmge_pah/Analysis/analyses/SOX18/covariates.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#removed all the NA's from Neutrophils column
metabolite1 <- metabolite %>% drop_na(Neutrophils)
metabolite1EUR=metabolite1[which(metabolite1$ancestry=="EA"),]
write.table(metabolite1EUR, "covariatesEUR.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#threegenesAFR and AMR  extracting/editing the IID

newfile = read.table(file='/N/project/mmge_pah/Analysis/analyses/new_survival_analysis/uchl1_phenotype_covars_v8_plus_key_metabolites_v2.txt', head= T, sep="\t")
newfile.AFR = newfile=newfile[which(newfile$ancestry=="AA"),]
library(dplyr)
new.AFR = newfile.AFR %>% select(IID, phenoID, self.report.ancestry, ancestry)
IID.AFR = new.AFR %>% select(IID)
write.table(IID.AFR, "IID.AFR.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#AMR extracting
newfile.AMR = newfile=newfile[which(newfile$ancestry=="AMR"),]
new.AMR = newfile.AMR %>% select(IID, phenoID, self.report.ancestry, ancestry)
IID.AMR = new.AMR %>% select(IID)
write.table(IID.AMR, "IID.AMR.txt", sep = "\t", row.names = FALSE, quote = FALSE)


library(gtools)
library(dplyr) 
library(tidyverse)
#boxplot for mpap
testmpap <- ev.test

#removing -inf for CABP4
testmpap<-  %>%
  filter(is.finite(CABP4))

boxplot(mpap)~ quantcut(CABP4, 3), ancestry, data=ev.test)

boxplot(mPAP ~ quantcut(CABP4, 3), data = testmpap, subset = ancestry)
table(quantcut(ev.test$CABP4, 3))

ev.test$CABP4.tert<-as.numeric(quantcut(ev.test$CABP4, 3))

boxplot(mPAP~ CABP4.tert,data=ev.test)

boxplot(mPAP ~ CABP4.tert, data = ev.test, subset = gender)

ggplot(ev.test, aes(y=mPAP, x=as.factor(CABP4.tert), fill=ancestry)) + 
  geom_boxplot() 
labs(fill="Ancestry", title="mPAP Severity", x=" CABP4 Tertiles", y="mPAP mmHg")

ggplot(ev.test, aes(x=CABP4.tert, y=mPAP, fill= interaction(ancestry)))+
  geom_boxplot()+
  labs(fill="ancestry", title="mPAP Severity", x=" CABP4 Tertiles", y="mPAP")
