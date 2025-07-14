library(gtools)
library(dplyr) 
library(tidyverse)

#locate the column # given the gene name
which(colnames(ev.test.EUR) == "PGAM1")

#loaded the ev.phenotypes for ancestries 
ev.test.EUR <- read.csv("EURlarge_file.csv")
ev.test.AA<- read.csv("AFRLarge_file.csv")
ev.test.AMR <- read.csv("AMRLarge_file.csv")

#log10 for datasets 
ev.test.EUR.log <- ev.test.EUR
ev.test.AA.log <- ev.test.AA
ev.test.AMR.log <- ev.test.AMR 

#location column number
which(colnames(ev.test.AA) == "CABP4")

#log10 for EUR the gene expressions
ev.test.EUR.log$FHL2 = log10(ev.test.EUR$FHL2)
ev.test.EUR.log$REXO2 = log10(ev.test.EUR$REXO2)
ev.test.EUR.log$CABP4 = log10(ev.test.EUR$CABP4)

#log10 for AA the gene expressions
ev.test.AA.log$FHL2 = log10(ev.test.AA$FHL2)
ev.test.AA.log$REXO2 = log10(ev.test.AA$REXO2)
ev.test.AA.log$CABP4 = log10(ev.test.AA$CABP4)

#log10 for AA the gene expressions
ev.test.AMR.log$FHL2 = log10(ev.test.AMR$FHL2)
ev.test.AMR.log$REXO2 = log10(ev.test.AMR$REXO2)
ev.test.AMR.log$CABP4 = log10(ev.test.AMR$CABP4)

#linear model regression for EUR 
#sixMWD.in.meter for EUR 
lm(sixMWD.in.meter~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log)
summary(lm(sixMWD.in.meter~age.diagnosis+sex+PAH.x+cell.frac+FHL2, data=ev.test.EUR.log))

#mRAP_EUR
lm(mRAP~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log)
summary(lm(mRAP~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log))

#mPAP_ CI 
lm(mPAP~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log)
summary(lm(mPAP~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log))

#CI EUR
lm(CI~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log)
summary(lm(CI~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log))

#CO EUR
lm(CO~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log)
summary(lm(CO~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log))

#CI.2 EUR
lm(CI.2~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log)
summary(lm(CI.2~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log))

#PVRi_Thermo EUR
lm(PVRi_Thermo~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log)
summary(lm(PVRi_Thermo~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log))

#RVSW_Thermo EUR
lm(RVSW_Thermo~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log)
summary(lm(RVSW_Thermo~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log))

#RVSWi_Thermo EUR
lm(RVSWi_Thermo~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log)
summary(lm(RVSWi_Thermo~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log))

#linear model regression for AFR
#sixMWD.in.meter for AA
lm(sixMWD.in.meter~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AA.log)
summary(lm(sixMWD.in.meter~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AA.log))

#mRAP_AA
lm(mRAP~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AA.log)
summary(lm(mRAP~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AA.log))

#mPAP_ AA
lm(mPAP~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AA.log)
summary(lm(mPAP~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AA.log))

#CI AA
lm(CI~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AA.log)
summary(lm(CI~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AA.log))

#CO AA
lm(CO~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AA.log)
summary(lm(CO~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AA.log))

#CI.2 AA
lm(CI.2~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AA.log)
summary(lm(CI.2~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AA.log))

#PVRi_Thermo AA
lm(PVRi_Thermo~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AA.log)
summary(lm(PVRi_Thermo~age.diagnosis+EV1+EV2+EV3+EV4+EV5+sex+cell.frac+FHL2, data=ev.test.AA.log))
 
#RVSW_Thermo AA
lm(RVSW_Thermo~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AA.log)
summary(lm(RVSW_Thermo~age.diagnosis+EV1+EV2+EV3+EV4+EV5+sex+cell.frac+FHL2, data=ev.test.AA.log))

#RVSWi_Thermo AA
lm(RVSWi_Thermo~age.diagnosis+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AA.log)
summary(lm(RVSWi_Thermo~age.diagnosis+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AA.log))

#linear model regression for AMR
#sixMWD.in.meter for AMR
lm(sixMWD.in.meter~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AMR.log)
summary(lm(sixMWD.in.meter~age.diagnosis+EV1+EV2+EV3+EV4+EV5+sex+cell.frac+FHL2, data=ev.test.AMR.log))

#mRAP_AMR
lm(mRAP~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AMR.log)
summary(lm(mRAP~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AMR.log))

#mPAP_ CI AMR
lm(mPAP~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AMR.log)
summary(lm(mPAP~age.diagnosis+EV1+EV2+EV3+EV4+EV5+sex+cell.frac+FHL2, data=ev.test.AMR.log))

#CI AMR
lm(CI~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AMR.log)
summary(lm(CI~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AMR.log))

#CO AMR
lm(CO~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AMR.log)
summary(lm(CO~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AMR.log))

#CI.2 AMR
lm(CI.2~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AMR.log)
summary(lm(CI.2~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AMR.log))

#PVRi_Thermo AMR
lm(PVRi_Thermo~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AMR.log)
summary(lm(PVRi_Thermo~age.diagnosis+EV1+EV2+EV3+EV4+EV5+sex+cell.frac+FHL2, data=ev.test.AMR.log))

#RVSW_Thermo AMR
lm(RVSW_Thermo~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AMR.log)
summary(lm(RVSW_Thermo~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AMR.log))

#RVSWi_Thermo AMR
lm(RVSWi_Thermo~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AMR.log)
summary(lm(RVSWi_Thermo~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.AMR.log))



#linear model regression for EUR REXO2
#sixMWD.in.meter for EUR 
lm(sixMWD.in.meter~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.EUR.log)
summary(lm(sixMWD.in.meter~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.EUR.log))

#mRAP_EUR REXO2
lm(mRAP~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+FHL2, data=ev.test.EUR.log)
summary(lm(mRAP~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.EUR.log))

#mPAP_ EUR REXO2
lm(mPAP~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.EUR.log)
summary(lm(mPAP~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.EUR.log))

#CI EUR REXO2
lm(CI~age.diagnosis+sex+PAH.x+cell.frac+REXO2, data=ev.test.EUR.log)
summary(lm(CI~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.EUR.log))

#CO EUR REXO2
lm(CO~age.diagnosis+sex+PAH.x+cell.frac+REXO2, data=ev.test.EUR.log)
summary(lm(CO~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.EUR.log))

#CI.2 EUR REXO2
lm(CI.2~age.diagnosis+sex+PAH.x+cell.frac+REXO2, data=ev.test.EUR.log)
summary(lm(CI.2~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.EUR.log))

#PVRi_Thermo EUR REXO2
lm(PVRi_Thermo~age.diagnosis+sex+PAH.x+cell.frac+REXO2, data=ev.test.EUR.log)
summary(lm(PVRi_Thermo~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.EUR.log))

#RVSW_Thermo EUR REXO2
lm(RVSW_Thermo~age.diagnosis+sex+PAH.x+cell.frac+REXO2, data=ev.test.EUR.log)
summary(lm(RVSW_Thermo~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.EUR.log))

#RVSWi_Thermo EUR  REXO2
lm(RVSWi_Thermo~age.diagnosis+sex+PAH.x+cell.frac+REXO2, data=ev.test.EUR.log)
summary(lm(RVSWi_Thermo~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.EUR.log))

#linear model regression for AFR
#sixMWD.in.meter for AA REXO2 
lm(sixMWD.in.meter~age.diagnosis+sex+cell.frac+REXO2 , data=ev.test.AA.log)
summary(lm(sixMWD.in.meter~age.diagnosis+EV1+EV2+EV3+EV4+EV5+sex+cell.frac+REXO2 , data=ev.test.AA.log))

#mRAP_AA REXO2 
lm(mRAP~age.diagnosis+sex+cell.frac+REXO2 , data=ev.test.AA.log)
summary(lm(mRAP~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2 , data=ev.test.AA.log))

#mPAP_ CI AA REXO2 
lm(mPAP~age.diagnosis+sex+cell.frac+REXO2 , data=ev.test.AA.log)
summary(lm(mPAP~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2 , data=ev.test.AA.log))

#CI AA REXO2 
lm(CI~age.diagnosis+sex+cell.frac+REXO2, data=ev.test.AA.log)
summary(lm(CI~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.AA.log))

#CO AA REXO2 
lm(CO~age.diagnosis+sex+cell.frac+REXO2, data=ev.test.AA.log)
summary(lm(CO~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.AA.log))

#CI.2 AA REXO2 
lm(CI.2~age.diagnosis+sex+cell.frac+REXO2, data=ev.test.AA.log)
summary(lm(CI.2~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.AA.log))

#PVRi_Thermo AA REXO2 
lm(PVRi_Thermo~age.diagnosis+sex+cell.frac+REXO2, data=ev.test.AA.log)
summary(lm(PVRi_Thermo~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.AA.log))

#RVSW_Thermo AA REXO2 
lm(RVSW_Thermo~age.diagnosis+sex+cell.frac+REXO2, data=ev.test.AA.log)
summary(lm(RVSW_Thermo~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.AA.log))

#RVSWi_Thermo AA REXO2 
lm(RVSWi_Thermo~age.diagnosis+cell.frac+REXO2, data=ev.test.AA.log)
summary(lm(RVSWi_Thermo~age.diagnosis+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.AA.log))

#linear model regression for AMR REXO2
#sixMWD.in.meter for AMR REXO2
lm(sixMWD.in.meter~age.diagnosis+sex+cell.frac+REXO2, data=ev.test.AMR.log)
summary(lm(sixMWD.in.meter~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.AMR.log))

#mRAP_AMR REXO2
lm(mRAP~age.diagnosis+sex+cell.frac+REXO2, data=ev.test.AMR.log)
summary(lm(mRAP~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.AMR.log))

#mPAP_ CI AMR REXO2
lm(mPAP~age.diagnosis+sex+cell.frac+REXO2, data=ev.test.AMR.log)
summary(lm(mPAP~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.AMR.log))

#CI AMR REXO2
lm(CI~age.diagnosis+sex+cell.frac+REXO2, data=ev.test.AMR.log)
summary(lm(CI~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.AMR.log))

#CO AMR REXO2
lm(CO~age.diagnosis+sex+cell.frac+REXO2, data=ev.test.AMR.log)
summary(lm(CO~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.AMR.log))

#CI.2 AMR REXO2
lm(CI.2~age.diagnosis+sex+cell.frac+REXO2, data=ev.test.AMR.log)
summary(lm(CI.2~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5++cell.frac+REXO2, data=ev.test.AMR.log))

#PVRi_Thermo AMR REXO2
lm(PVRi_Thermo~age.diagnosis+sex+cell.frac+REXO2, data=ev.test.AMR.log)
summary(lm(PVRi_Thermo~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.AMR.log))

#RVSW_Thermo AMR REXO2
lm(RVSW_Thermo~age.diagnosis+sex+cell.frac+REXO2, data=ev.test.AMR.log)
summary(lm(RVSW_Thermo~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.AMR.log))

#RVSWi_Thermo AMR REXO2
lm(RVSWi_Thermo~age.diagnosis+sex+cell.frac+REXO2, data=ev.test.AMR.log)
summary(lm(RVSWi_Thermo~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+REXO2, data=ev.test.AMR.log))

#removing -inf for CABP4
ev.test.AA.log <- ev.test.AA.log %>%
  filter(is.finite(CABP4))
ev.test.EUR.log <- ev.test.EUR.log %>%
  filter(is.finite(CABP4))
ev.test.AMR.log <- ev.test.AMR.log %>%
  filter(is.finite(CABP4))


#linear model regression for EUR CABP4
#sixMWD.in.meter for EUR 
lm(sixMWD.in.meter~age.diagnosis+sex+PAH.x+cell.frac+CABP4, data=ev.test.EUR.log)
summary(lm(sixMWD.in.meter~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.EUR.log))

#mRAP_EUR CABP4
lm(mRAP~age.diagnosis+sex+PAH.x+cell.frac+CABP4, data=ev.test.EUR.log)
summary(lm(mRAP~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.EUR.log))

#mPAP_ EUR CABP4
lm(mPAP~age.diagnosis+sex+PAH.x+cell.frac+CABP4, data=ev.test.EUR.log)
summary(lm(mPAP~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.EUR.log))

#CI EUR CABP4
lm(CI~age.diagnosis+sex+PAH.x+cell.frac+CABP4, data=ev.test.EUR.log)
summary(lm(CI~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.EUR.log))

#CO EUR CABP4
lm(CO~age.diagnosis+sex+PAH.x+cell.frac+CABP4, data=ev.test.EUR.log)
summary(lm(CO~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.EUR.log))

#CI.2 EUR CABP4
lm(CI.2~age.diagnosis+sex+PAH.x+cell.frac+CABP4, data=ev.test.EUR.log)
summary(lm(CI.2~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.EUR.log))

#PVRi_Thermo EUR CABP4
lm(PVRi_Thermo~age.diagnosis+sex+PAH.x+cell.frac+CABP4, data=ev.test.EUR.log)
summary(lm(PVRi_Thermo~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.EUR.log))

#RVSW_Thermo EUR CABP4
lm(RVSW_Thermo~age.diagnosis+sex+PAH.x+cell.frac+CABP4, data=ev.test.EUR.log)
summary(lm(RVSW_Thermo~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.EUR.log))

#RVSWi_Thermo EUR  CABP4
lm(RVSWi_Thermo~age.diagnosis+sex+PAH.x+cell.frac+CABP4, data=ev.test.EUR.log)
summary(lm(RVSWi_Thermo~age.diagnosis+sex+PAH.x+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.EUR.log))

#linear model regression for CABP4
#sixMWD.in.meter for AA CABP4 
lm(sixMWD.in.meter~age.diagnosis+sex+cell.frac+CABP4 , data=ev.test.AA.log)
summary(lm(sixMWD.in.meter~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4 , data=ev.test.AA.log))

#mRAP_AA CABP4
lm(mRAP~age.diagnosis+sex+cell.frac+CABP4 , data=ev.test.AA.log)
summary(lm(mRAP~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4 , data=ev.test.AA.log))

#mPAP_AA CABP4
lm(mPAP~age.diagnosis+sex+cell.frac+CABP4 , data=ev.test.AA.log)
summary(lm(mPAP~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4 , data=ev.test.AA.log))

#CI AA CABP4
lm(CI~age.diagnosis+sex+cell.frac+CABP4, data=ev.test.AA.log)
summary(lm(CI~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.AA.log))

#CO AA CABP4 
lm(CO~age.diagnosis+sex+cell.frac+CABP4, data=ev.test.AA.log)
summary(lm(CO~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.AA.log))

#CI.2 AA CABP4 
lm(CI.2~age.diagnosis+sex+cell.frac+CABP4, data=ev.test.AA.log)
summary(lm(CI.2~age.diagnosis+EV1+EV2+EV3+EV4+EV5+sex+cell.frac+CABP4, data=ev.test.AA.log))

#PVRi_Thermo AA CABP4
lm(PVRi_Thermo~age.diagnosis+sex+cell.frac+CABP4, data=ev.test.AA.log)
summary(lm(PVRi_Thermo~age.diagnosis+EV1+EV2+EV3+EV4+EV5+sex+cell.frac+CABP4, data=ev.test.AA.log))

#RVSW_Thermo AA CABP4
lm(RVSW_Thermo~age.diagnosis+sex+cell.frac+CABP4, data=ev.test.AA.log)
summary(lm(RVSW_Thermo~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.AA.log))

#RVSWi_Thermo AA CABP4
lm(RVSWi_Thermo~age.diagnosis+cell.frac+CABP4, data=ev.test.AA.log)
summary(lm(RVSWi_Thermo~age.diagnosis+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.AA.log))

#linear model regression for AMR CABP4
#sixMWD.in.meter for AMR CABP4
lm(sixMWD.in.meter~age.diagnosis+sex+cell.frac+CABP4, data=ev.test.AMR.log)
summary(lm(sixMWD.in.meter~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.AMR.log))

#mRAP_AMR CABP4
lm(mRAP~age.diagnosis+sex+cell.frac+CABP4, data=ev.test.AMR.log)
summary(lm(mRAP~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.AMR.log))

#mPAP_ AMR CABP4
lm(mPAP~age.diagnosis+sex+cell.frac+CABP4, data=ev.test.AMR.log)
summary(lm(mPAP~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.AMR.log))

#CI AMR CABP4
lm(CI~age.diagnosis+sex+cell.frac+CABP4, data=ev.test.AMR.log)
summary(lm(CI~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.AMR.log))

#CO AMR CABP4
lm(CO~age.diagnosis+sex+cell.frac+CABP4, data=ev.test.AMR.log)
summary(lm(CO~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.AMR.log))

#CI.2 AMR CABP4
lm(CI.2~age.diagnosis+sex+cell.frac+CABP4, data=ev.test.AMR.log)
summary(lm(CI.2~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.AMR.log))

#PVRi_Thermo AMR CABP4
lm(PVRi_Thermo~age.diagnosis+sex+cell.frac+CABP4, data=ev.test.AMR.log)
summary(lm(PVRi_Thermo~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.AMR.log))

#RVSW_Thermo AMR CABP4
lm(RVSW_Thermo~age.diagnosis+sex+cell.frac+CABP4, data=ev.test.AMR.log)
summary(lm(RVSW_Thermo~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.AMR.log))

#RVSWi_Thermo AMR CABP4
lm(RVSWi_Thermo~age.diagnosis+sex+cell.frac+CABP4, data=ev.test.AMR.log)
summary(lm(RVSWi_Thermo~age.diagnosis+sex+EV1+EV2+EV3+EV4+EV5+cell.frac+CABP4, data=ev.test.AMR.log))





#Metabolite 
summary(lm(MALATE..M.H.~age.diagnosis+sex+PVR.x+CABP4, data=ev.test.EUR.log))

# Severity analysis in the EA
severe=function(in.dat,outcome,covars){
  test.set=get(in.dat)
  test.sev=lm(as.formula(paste(outcome,"~",gsub(",","+",covars),"+",sep="")),data=test.set)
  output=data.frame(matrix(nrow=1,ncol=7,NA))
  output[1,1]=outcome
  output[1,2]=length(test.sev$residuals)
  output[1,3]=NA
  output[1,4]=summary(test.sev)$coefficient[covars,1]
  output[1,5]=summary(test.sev)$coefficient[covars,2]
  output[1,6]=summary(test.sev)$coefficient[covars,3]
  output[1,7]=summary(test.sev)$coefficient[covars,4]
  return(output)
}


summary(lm(PYRUVATE..M.H.~age.diagnosis+sex+PVR.x+FHL2, data=ev.test.EUR.log))
# Define the function
summarize_pvalues <- function(data, cols, predictor, covariates) {
  # Create an empty results table
  results <- data.frame(
    Outcome = character(),  # Column name (outcome)
    P_Value_FHL2 = numeric() # P-value for FHL2
  )
  
  # Loop through the specified columns
  for (col_index in cols) {
    # Get the column name for the outcome
    outcome <- colnames(data)[col_index]
    
    # Create the formula dynamically
    formula <- as.formula(paste(outcome, "~", paste(covariates, collapse = "+"), "+", predictor))
    
    # Try to fit the model and extract the p-value
    tryCatch({
      model <- lm(formula, data = data)
      p_value <- summary(model)$coefficients[predictor, "Pr(>|t|)"] # Extract p-value for FHL2
      # Append the results to the table
      results <- rbind(results, data.frame(Outcome = outcome, P_Value_FHL2 = p_value))
    }, error = function(e) {
      # Handle cases where the model fails (e.g., NA or singular matrix)
      results <- rbind(results, data.frame(Outcome = outcome, P_Value_REXO2 = NA))
    })
  }
  
  # Return the results table
  return(results)
}

# Example usage
# Assuming `ev.test.EUR.log` is your dataset
cols_to_test <- 93:113
covariates <- c("age.diagnosis", "sex", "PAH.x", "PVR.x") # List of covariates
covariates2 <- c("age.diagnosis", "sex", "PVR.x")
predictor <- "FHL2" # The predictor of interest
predictor2 <- "REXO2"
predictor3 <- "CABP4"

# Run the function
#EUR results 
results_table <- summarize_pvalues(ev.test.EUR.log, cols_to_test, predictor, covariates)
results_table2 <- summarize_pvalues(ev.test.EUR.log, cols_to_test, predictor2, covariates)
results_table3 <- summarize_pvalues(ev.test.EUR.log, cols_to_test, predictor3, covariates)
#AFR results 
results_table4 <- summarize_pvalues(ev.test.AA.log, cols_to_test, predictor, covariates2)
results_table5 <- summarize_pvalues(ev.test.AA.log, cols_to_test, predictor2, covariates2)
results_table6 <- summarize_pvalues(ev.test.AA.log, cols_to_test, predictor3, covariates2)
#AMR Results 
results_table7 <- summarize_pvalues(ev.test.AMR.log, cols_to_test, predictor, covariates2)
results_table8 <- summarize_pvalues(ev.test.AMR.log, cols_to_test, predictor2, covariates2)
results_table9 <- summarize_pvalues(ev.test.AMR.log, cols_to_test, predictor3, covariates2)

# Print the results
print(results_table)
print(results_table2)
print(results_table7)


# Optionally save to a file
write.table(results_table, file = "EUR_FHL2_pvalues_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(results_table2, file = "EUR_REXO2_pvalues_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(results_table3, file = "EUR_CABP4_pvalues_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

write.table(results_table4, file = "AFR_FHL2_pvalues_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(results_table5, file = "AFR_REXO2_pvalues_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(results_table6, file = "AFR_CABP4_pvalues_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

write.table(results_table7, file = "AMR_FHL2_pvalues_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(results_table8, file = "AMR_REXO2_pvalues_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(results_table9, file = "AMR_CABP4_pvalues_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


#test 
#Metabolite 
summary(lm(PYRUVATE..M.H.~age.diagnosis+sex+PVR.x+FHL2, data=ev.test.AMR.log))



summarize_fhl2_result <- function(data, cols, predictor, covariates) {
  # Create an empty results data frame
  results <- data.frame(
    Outcome = character(),
    Estimate_FHL2 = numeric(),
    StdError_FHL2 = numeric(),
    P_Value_FHL2 = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through the specified columns
  for (col_index in cols) {
    # Get the column name for the outcome
    outcome <- colnames(data)[col_index]
    
    # Create the formula dynamically
    formula <- as.formula(paste(outcome, "~", paste(covariates, collapse = "+"), "+", predictor))
    
    # Try to fit the model and extract results
    tryCatch({
      model <- lm(formula, data = data)
      summary_model <- summary(model)
      
      # Extract coefficients for FHL2
      estimate <- summary_model$coefficients[predictor, "Estimate"]
      std_error <- summary_model$coefficients[predictor, "Std. Error"]
      p_value <- summary_model$coefficients[predictor, "Pr(>|t|)"]
       
      # Append the results to the data frame
      results <- rbind(
        results,
        data.frame(Outcome = outcome, Estimate_FHL2 = estimate, StdError_FHL2 = std_error, P_Value_FHL2 = p_value)
      )
    }, error = function(e) {
      # Handle errors (e.g., NA values or singular matrix)
      results <- rbind(
        results,
        data.frame(Outcome = outcome, Estimate_FHL2 = NA, StdError_FHL2 = NA, P_Value_FHL2 = NA)
      )
    })
  }
  
  # Return the results
  return(results)
}

# Example Usage
# Define the range of columns for outcomes (93 to 113)
cols_to_test <- 93:113

# Define the covariates
covariates <- c("age.diagnosis", "sex", "PAH.x", "PVR.x")
covariates1 <- c("age.diagnosis", "sex", "PVR.x")

# Define the predictor of interest
predictor <- "FHL2"
predictor2 <- "CABP4"
predictor3 <- "REXO2"

# Assuming `ev.test.EUR.log` is your dataset
results_table <- summarize_fhl2_result(ev.test.EUR.log, cols_to_test, predictor, covariates)
results_table1 <- summarize_fhl2_result(ev.test.EUR.log, cols_to_test, predictor2, covariates)
results_table2 <- summarize_fhl2_result(ev.test.EUR.log, cols_to_test, predictor3, covariates)

# Assuming `ev.test.AA.log` is your dataset
results_table3 <- summarize_fhl2_result(ev.test.AA.log, cols_to_test, predictor, covariates1)
results_table4 <- summarize_fhl2_result(ev.test.AA.log, cols_to_test, predictor2, covariates1)
results_table5 <- summarize_fhl2_result(ev.test.AA.log, cols_to_test, predictor3, covariates1)

# Assuming `ev.test.AMR.log` is your dataset
results_table6 <- summarize_fhl2_result(ev.test.AMR.log, cols_to_test, predictor, covariates1)
results_table7<- summarize_fhl2_result(ev.test.AMR.log, cols_to_test, predictor2, covariates1)
results_table8<- summarize_fhl2_result(ev.test.AMR.log, cols_to_test, predictor3, covariates1)

# Print the results
print(results_table)
print(results_table8)

# Optionally save the results to a file
write.table(results_table, file = "FHL2_results_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(results_table1, file = "CABP4_results_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(results_table2, file = "REXO2_results_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

write.table(results_table3, file = "FHL2_AFR_results_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(results_table4, file = "CAPB4_AFR_results_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(results_table5, file = "REXO2_AFR_results_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

write.table(results_table6, file = "FHL2_AMR_results_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(results_table7, file = "CABP4_AMR_results_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(results_table8, file = "REXO2_AMR_results_table.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


# total number of NA'S: sum(is.na(data$column_name))
sum(is.na(ev.test.EUR$IID))

#makes a table for transplant including NA's
table(ev.test.EUR$transplant, useNA="ifany")

#includes all NA's in Transplant and primary.status 
table(ev.test.AA$transplant, ev.test.AA$primary.status, useNA = "ifany")

library(dplyr)
new_df_AA <- ev.test.AA %>% select(CI, primary.status, Alive, transplant)

#total amount of zeros in a df second column
total_zeros <- sum(new_df_AA[, 2] == 0, na.rm = TRUE)

table(new_df_AA$transplant, useNA="ifany")

#new dr for EUR on CI, primary.status, Alive, transplant
new_df_EUR <- ev.test.EUR %>% select(CI, primary.status, Alive, transplant)

new_df_EUR$CI2 <- ifelse(new_df_EUR$CI>0, 1) 

new_df_AMR <- ev.test.AMR %>% select(CI, primary.status, Alive, transplant)

new_df_AMR$CI2 <- ifelse(new_df_AMR$CI>0, 1) 




#old function

library(survival)
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

#Running a survival analysis for gene PGAM1 for EUR
EUR_PGAM1=enroll.run("ev.test.EUR",colnames(ev.test.EUR)[125],"sex,age.diagnosis,PAH.x,PVR.x,EV1,EV2,EV3,EV3,EV4,EV5,cell.frac")
for(i in 13008){#dim(ev.test.EUR)[2]){ 
  EUR_PGAM1=try(rbind(EUR_PGAM1,try(enroll.run("ev.test.EUR",colnames(ev.test.EUR)[i],"sex,age.diagnosis,PAH.x,PVR.x,EV1,EV2,EV3,EV4,EV5,cell.frac"))))
}

#how to delete a row from a df
EUR_PGAM1 <- EUR_PGAM1[-1, ]

write.csv(EUR_PGAM1, "EUR_PGAM1.csv", row.names=TRUE) 

#Running a survival analysis for gene PGAM1 for AFR
AFR_PGAM1=enroll.run("ev.test.AA",colnames(ev.test.AA)[124],"sex,age.diagnosis,PVR.x,EV1,EV2,EV3,EV3,EV4,EV5,cell.frac")
for(i in 13008){#dim(ev.test.AA)[2]){ 
  AFR_PGAM1=try(rbind(AFR_PGAM1,try(enroll.run("ev.test.AA",colnames(ev.test.AA)[i],"sex,age.diagnosis,PVR.x,EV1,EV2,EV3,EV4,EV5,cell.frac"))))
}

#how to delete a row from a df
AFR_PGAM1 <- AFR_PGAM1[-1, ]

write.csv(AFR_PGAM1, "AFR_PGAM1.csv", row.names=TRUE) 

AMR_PGAM1=enroll.run("ev.test.AMR",colnames(ev.test.AMR)[124],"sex,age.diagnosis,PVR.x,EV1,EV2,EV3,EV3,EV4,EV5,cell.frac")
for(i in 13008){#dim(ev.test.AMR)[2]){ 
  AMR_PGAM1=try(rbind(AMR_PGAM1,try(enroll.run("ev.test.AMR",colnames(ev.test.AMR)[i],"sex,age.diagnosis,PVR.x,EV1,EV2,EV3,EV4,EV5,cell.frac"))))
}

#how to delete a row from a df
AMR_PGAM1 <- AMR_PGAM1[-1, ]

write.csv(AMR_PGAM1, "AMR_PGAM1.csv", row.names=TRUE) 
