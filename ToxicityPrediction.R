library(readxl)
library(dplyr)
library(readr)
library(tidyr)
library(ROCit)

#Grade >= 3 vs others (compare aSMA and MHU)
#PS vs sarcopenia correlation

setwd("D:/OneDrive - Imperial College London/Documents/Sarcopenia/Toxicity Study/New_workspace")

Milan <- read_excel("Milan - Worksheet.xlsx")
Imperial <- read_excel("Imperial - Worksheet.xlsx")
Imperial_sarc <- read_excel("Sarcopenia-measurements.xltx")

#Milan <- read_excel("/home/mitchchen/Downloads/Sarcopenia_workspace/Milan - Worksheet.xlsx")
#Imperial <- read_excel("/home/mitchchen/Downloads/Sarcopenia_workspace/Imperial - Worksheet.xlsx")
#Imperial_sarc <- read_excel("/home/mitchchen/Downloads/Sarcopenia_workspace/Sarcopenia-measurements.xltx")

Imperial$ID <- parse_number(Imperial$...1)

##############Processing sarcopenia measurements ##################################

SMA <- aggregate(`MuscleArea (cmÂ²)`~ ID + StudyDateTime, Imperial_sarc, FUN=c)
MHU <- aggregate(`MuscleMeanAttenuation (HU)`~ ID + StudyDateTime, Imperial_sarc, FUN=c)

sarc_reading <- data.frame()

for (i in 1:dim(SMA)[1]){
  sarc_reading [i,1:2] <- MHU[i,1:2]
  index <- which.min(unlist(MHU[i,3]))
  sarc_reading [i,3] <- unlist(SMA[i,3])[index]  
  sarc_reading [i,4] <- unlist(MHU[i,3])[index] 
}

sarc_reading <- sarc_reading %>% arrange(ID, StudyDateTime)

baseline_sarc <- data.frame()
fu_sarc <- data.frame()
j <- 1
k <- 1
current_id <- 0

#Split into baseline and followup vectors
for (i in 1:dim(sarc_reading)[1]){
  if (sarc_reading[i,1] != current_id){
    current_id <- sarc_reading[i,1]
    baseline_sarc[j,1:4] <- sarc_reading[i,1:4]
    j <- j + 1
  }
  else {
    fu_sarc[k,1:4] <- sarc_reading[i,1:4]
    k <- k + 1
  }
}

#Imperial <- subset(Imperial, select = -c(SMA, MHU, baselinedate))
Imperial$SMA <- NA
Imperial$MHU <- NA
Imperial$baselinedate <- NA

for (i in 1:dim(baseline_sarc)[1]){
  idex <- which(Imperial$ID == baseline_sarc[i,1]) 
  Imperial$baselinedate[idex] <- baseline_sarc[i,2]
  Imperial$SMA[idex] <- baseline_sarc[i,3]
  Imperial$MHU[idex] <- baseline_sarc[i,4]
}

#Up to three followup scans
Imperial$fudate1 <- NA
Imperial$fuMHU1 <- NA
Imperial$fuSMA1 <- NA
Imperial$fudate2 <- NA
Imperial$fuMHU2 <- NA
Imperial$fuSMA2 <- NA
Imperial$fudate3 <- NA
Imperial$fuMHU3 <- NA
Imperial$fuSMA3 <- NA

current_id <- 0

for (i in 1:dim(fu_sarc)[1]){
  if (current_id != fu_sarc[i,1]){
    current_id <- fu_sarc[i,1]
    idex <- which(Imperial$ID == fu_sarc[i,1]) 
    Imperial$fudate1[idex] <- fu_sarc[i,2]
    Imperial$fuSMA1[idex] <- fu_sarc[i,3]
    Imperial$fuMHU1[idex] <- fu_sarc[i,4]
  }
  else{
    idex <- which(Imperial$ID == fu_sarc[i,1])
    if (is.na(Imperial$fudate2[idex])){
      Imperial$fudate2[idex] <- fu_sarc[i,2]
      Imperial$fuSMA2[idex] <- fu_sarc[i,3]
      Imperial$fuMHU2[idex] <- fu_sarc[i,4]
    }
    else{
      Imperial$fudate3[idex] <- fu_sarc[i,2]
      Imperial$fuSMA3[idex] <- fu_sarc[i,3]
      Imperial$fuMHU3[idex] <- fu_sarc[i,4]  
    }
  }
}

#Calculates BMI, SMA adjusted by height a
Milan$BMI <- Milan$Weight/((Milan$Height/100)^2)
Milan$aSMA <- Milan$SMA/((Milan$Height/100)^2)
Milan$BSA <- 0
Milan$BSA <- ((Milan$Weight*Milan$Height)/3600)^(1/2)

Imperial$BMI <- Imperial$Weight/((Imperial$Height/100)^2)
Imperial$aSMA <- Imperial$SMA/((Imperial$Height/100)^2)
Imperial$aSMAfu1 <- Imperial$fuSMA1/((Imperial$Height/100)^2)
Imperial$aSMAfu2 <- Imperial$fuSMA2/((Imperial$Height/100)^2)
Imperial$aSMAfu3 <- Imperial$fuSMA3/((Imperial$Height/100)^2)

Imperial$BSA <- 0
Imperial$BSA <- ((Imperial$Weight*Milan$Height)/3600)^(1/2)


#Replace NA with 0s
Milan[,20:42] <- Milan[,20:42] %>% replace(is.na(.), 0)
Imperial[,14:91] <- Imperial[,14:91] %>% replace(is.na(.), 0)

#Replace \\ with NAs
Imperial$PS[Imperial$PS == "\\"] <- NA

#Replace genders in Milan
Milan$Sex[Milan$Sex == "M"] <- 1
Milan$Sex[Milan$Sex == "F"] <- 2


#############Cohort splits########################

#Split Imperial into Imperial_IO
ImperialIO <- Imperial[Imperial$Immunotherapy==1,]

#Split both cohorts into male and female for aSMA use
ImperialM <- Imperial[Imperial$Sex==1,]
ImperialF <- Imperial[Imperial$Sex==2,]

ImperialIOM <- ImperialIO[ImperialIO$Sex==1,]
ImperialIOF <- ImperialIO[ImperialIO$Sex==2,]

MilanM <- Milan[Milan$Sex==1,]
MilanF <- Milan[Milan$Sex==2,]


######Check case matching with propensity score##############
install.packages("MatchIt")       # Only if not already installed
install.packages("tableone")      # For balance diagnostics

library(MatchIt)
library(tableone)

#Match and merge
ImperialIO_df <- ImperialIO
names(ImperialIO_df)[names(ImperialIO_df) == "...1"] <- "Patient ID"
ImperialIO_df$cohort <- 0

Milan_df <- Milan
names(Milan_df) <- gsub("Performance Status", "PS", names(Milan_df))
Milan_df$cohort <- 1

common_cols <- intersect(colnames(ImperialIO_df), colnames(Milan_df))
merged_df <- rbind(ImperialIO_df[, common_cols], Milan_df[, common_cols])

merged_df$Age <- as.factor(merged_df$Age)
merged_df$Sex <- as.factor(merged_df$Sex)
merged_df$PS <- as.factor(merged_df$PS)
merged_df$Height <- as.factor(merged_df$Height)
merged_df$Weight <- as.factor(merged_df$Weight)
merged_df$cohort <- as.factor(merged_df$cohort)

merged_df <- merged_df[complete.cases(merged_df[c("Age", "Sex", "PS", "Height", "Weight")]), ]

# Estimate propensity score and match
match_model <- matchit(cohort ~ Age + Sex + PS + Height + Weight,
                       data = merged_df,
                       method = "nearest",        # Other options: "optimal", "genetic", "full"
                       ratio = 1)                 

# View a summary
summary(match_model)
matched_data <- match.data(match_model)

table1 <- CreateTableOne(vars = c("Age", "Sex", "PS", "Height", "Weight"),
                         strata = "cohort",
                         data = matched_data,
                         test = FALSE)
print(table1, smd = TRUE)

t.test(outcome ~ gro, data = matched_data)

################Multivariable Predictor###############################

#Do LASSO With aSMA, MHU, BMI, and BSA

install.packages("survival")
library(survival)

#Composite1Y metric

fit1 <- coxph(Surv(NewOS, Event) ~ aSMA + MHU + BMI,
              data=MilanF)

fit2 <- coxph(Surv(NewOS, Event) ~ aSMA + MHU + BMI,
              data=MilanM)

MilanF$Composite <- predict(fit1, newdata=MilanF, type = "lp")
MilanM$Composite <- predict(fit2, newdata=MilanM, type = "lp")

#Imperial
ImperialF$Composite <- predict(fit1, newdata=ImperialF, type = "lp")
ImperialM$Composite <- predict(fit2, newdata=ImperialM, type = "lp")

#ImperialIO
ImperialIOF$Composite <- predict(fit1, newdata=ImperialIOF, type = "lp")
ImperialIOM$Composite <- predict(fit2, newdata=ImperialIOM, type = "lp")

Milan <- rbind(MilanF, MilanM)
Imperial <- rbind (ImperialF, ImperialM)
ImperialIO <- rbind (ImperialIOF, ImperialIOM)

#Composite90D metric - don't use for now

fit1 <- coxph(Surv(NewOS, Event) ~ aSMA + MHU + BMI + BSA,
              data=MilanF)

fit2 <- coxph(Surv(NewOS, Event) ~ aSMA + MHU + BMI + BSA,
              data=MilanM)

MilanF$Composite2 <- predict(fit1, newdata=MilanF, type = "lp")
MilanM$Composite2 <- predict(fit2, newdata=MilanM, type = "lp")

#Imperial
ImperialF$Composite2 <- predict(fit1, newdata=ImperialF, type = "lp")
ImperialM$Composite2 <- predict(fit2, newdata=ImperialM, type = "lp")

#ImperialIO
ImperialIOF$Composite2 <- predict(fit1, newdata=ImperialIOF, type = "lp")
ImperialIOM$Composite2 <- predict(fit2, newdata=ImperialIOM, type = "lp")



##########Assign NewOS & Hospital Stay################################

#Milan

#Specify event up to 36 months (3 year OS)
Milan$Event <- 0 
Milan$Event <- ifelse(Milan$`Overall Survival (months)` < 36, 1, 0)

Milan$NewOS <- 0
Milan$NewOS <- ifelse(Milan$`Overall Survival (months)` < 36, Milan$`Overall Survival (months)`, 36) * 30.44

#Specify event up to 12 months (1 year OS)
Milan$Event <- 0 
Milan$Event <- ifelse(Milan$`Overall Survival (months)` < 12, 1, 0)

Milan$NewOS <- 0
Milan$NewOS <- ifelse(Milan$`Overall Survival (months)` < 12, Milan$`Overall Survival (months)`, 12) * 30.44

#Specify event up to 90 days
Milan$Event <- 0 
Milan$Event <- ifelse(Milan$`Overall Survival (months)` < 3, 1, 0)

Milan$NewOS <- 0
Milan$NewOS <- ifelse(Milan$`Overall Survival (months)` < 3, Milan$`Overall Survival (months)`, 3) * 30.44

#Imperial

Imperial$Event <- 0 
Imperial$Event <- ifelse(Imperial$OS <365, 1, 0)

Imperial$NewOS <- 0
Imperial$NewOS <- ifelse(Imperial$OS < 365, Imperial$OS, 365)

Imperial$Event <- 0 
Imperial$Event <- ifelse(Imperial$OS < 90, 1, 0)

Imperial$NewOS <- 0
Imperial$NewOS <- ifelse(Imperial$OS < 90, Imperial$OS, 90)

Imperial$HospitalStayB <-  0
Imperial$LongStay <- 0 

Imperial$HospitalStayB <- ifelse(Imperial$HospitalStay>0, 1, 0)
Imperial$LongStay <- ifelse(Imperial$HospitalStay>7, 1, 0)


#ImperialIO

ImperialIO$Event <- 0 
ImperialIO$Event <- ifelse(ImperialIO$OS < 365, 1, 0)

ImperialIO$NewOS <- 0
ImperialIO$NewOS <- ifelse(ImperialIO$OS < 365, ImperialIO$OS, 365)

ImperialIO$Event <- 0 
ImperialIO$Event <- ifelse(ImperialIO$OS < 90, 1, 0)

ImperialIO$NewOS <- 0
ImperialIO$NewOS <- ifelse(ImperialIO$OS < 90, Imperial$OS, 90)


##########Patient characteristics############################

#Summary statistics
mean (as.numeric(Imperial$Age))
sd (as.numeric(Imperial$Age))
max (as.numeric(Imperial$Age))
min (as.numeric(Imperial$Age))

tbl <- table(Imperial$Sex)
cbind(tbl,prop.table(tbl))

#Summary statistics
mean (as.numeric(!is.na(ImperialIO$PS)))
sd (as.numeric(!is.na(ImperialIO$PS)))
max (as.numeric(!is.na(ImperialIO$PS)))
min (as.numeric(!is.na(ImperialIO$PS)))

tbl <- table(!is.na(ImperialIO$PS))
cbind(tbl,prop.table(tbl))



Imperial_M <- Imperial[Imperial$Sex==1,]
Imperial_F <- Imperial[Imperial$Sex==2,]

ICL_temp_F <- Imperial_F[!is.na(Imperial_F$Height),]
ICL_temp_M <- Imperial_M[!is.na(Imperial_M$Height),]

mean (as.numeric(ICL_temp$Height))
sd (as.numeric(ICL_temp$Height))
max (as.numeric(ICL_temp$Height))
min (as.numeric(ICL_temp$Height))

wilcox.test(x = ICL_temp_M$Height, y = as.numeric(ICL_temp_F$Height),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

mean (as.numeric(ICL_temp_M$BMI))
sd (as.numeric(ICL_temp_M$BMI))
max (as.numeric(ICL_temp_M$BMI))
min (as.numeric(ICL_temp_M$BMI))

wilcox.test(x = ICL_temp_M$BMI, y = as.numeric(ICL_temp_F$BMI),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

ICL_temp_F <- Imperial_F[!is.na(Imperial_F$PS),]
ICL_temp_M <- Imperial_M[!is.na(Imperial_M$PS),]

tbl <- table(ICL_temp_M$PS)
cbind(tbl,prop.table(tbl))

wilcox.test(x = as.numeric(ICL_temp_M$PS), y = as.numeric(ICL_temp_F$PS),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

ICL_temp_F <- Imperial_F[!is.na(Imperial_F$aSMA),]
ICL_temp_M <- Imperial_M[!is.na(Imperial_M$aSMA),]

mean (as.numeric(ICL_temp_M$SMA))
sd (as.numeric(ICL_temp_M$SMA))
max (as.numeric(ICL_temp_M$SMA))
min (as.numeric(ICL_temp_M$SMA))

wilcox.test(x = as.numeric(ICL_temp_M$SMA), y = as.numeric(ICL_temp_F$SMA),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

mean (as.numeric(ICL_temp_M$aSMA))
sd (as.numeric(ICL_temp_M$aSMA))
max (as.numeric(ICL_temp_M$aSMA))
min (as.numeric(ICL_temp_M$aSMA))

wilcox.test(x = as.numeric(ICL_temp_M$aSMA), y = as.numeric(ICL_temp_F$aSMA),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

mean (as.numeric(ICL_temp_M$MHU))
sd (as.numeric(ICL_temp_M$MHU))
max (as.numeric(ICL_temp_M$MHU))
min (as.numeric(ICL_temp_M$MHU))

wilcox.test(x = as.numeric(ICL_temp_M$MHU), y = as.numeric(ICL_temp_F$MHU),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

wilcox.test(x = as.numeric(Imperial$Age), y = as.numeric(Milan$Age),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)


#ImperialIO

#Summary statistics
mean (as.numeric(ImperialIO$Age))
sd (as.numeric(ImperialIO$Age))
max (as.numeric(ImperialIO$Age))
min (as.numeric(ImperialIO$Age))

wilcox.test(x = as.numeric(Imperial$Age), y = as.numeric(ImperialIO$Age),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

tbl <- table(ImperialIO$Sex)
cbind(tbl,prop.table(tbl))

wilcox.test(x = as.numeric(Imperial$Sex), y = as.numeric(ImperialIO$Sex),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)


ICL_temp_F <- ImperialIOF[!is.na(ImperialIOF$Height),]
ICL_temp_M <- ImperialIOM[!is.na(ImperialIOM$Height),]

ICL_temp <- ICL_temp_F

mean (as.numeric(ICL_temp$Height))
sd (as.numeric(ICL_temp$Height))
max (as.numeric(ICL_temp$Height))
min (as.numeric(ICL_temp$Height))

wilcox.test(x = ICL_temp_M$Height, y = as.numeric(ICL_temp_F$Height),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

mean (as.numeric(ICL_temp$BMI))
sd (as.numeric(ICL_temp$BMI))
max (as.numeric(ICL_temp$BMI))
min (as.numeric(ICL_temp$BMI))

wilcox.test(x = ICL_temp_M$BMI, y = as.numeric(ICL_temp_F$BMI),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

ICL_temp_F <- ImperialIOF[!is.na(ImperialIOF$PS),]
ICL_temp_M <- ImperialIOM[!is.na(ImperialIOM$PS),]

tbl <- table(ICL_temp_M$PS)
cbind(tbl,prop.table(tbl))

tbl <- table(ICL_temp_F$PS)
cbind(tbl,prop.table(tbl))

wilcox.test(x = as.numeric(ICL_temp_M$PS), y = as.numeric(ICL_temp_F$PS),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

ICL_temp_F <- ImperialIOF[!is.na(ImperialIOF$aSMA),]
ICL_temp_M <- ImperialIOM[!is.na(ImperialIOM$aSMA),]

mean (as.numeric(ICL_temp_M$SMA))
sd (as.numeric(ICL_temp_M$SMA))
max (as.numeric(ICL_temp_M$SMA))
min (as.numeric(ICL_temp_M$SMA))

mean (as.numeric(ICL_temp_F$SMA))
sd (as.numeric(ICL_temp_F$SMA))
max (as.numeric(ICL_temp_F$SMA))
min (as.numeric(ICL_temp_F$SMA))

wilcox.test(x = as.numeric(ICL_temp_M$SMA), y = as.numeric(ICL_temp_F$SMA),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

mean (as.numeric(ICL_temp_M$aSMA))
sd (as.numeric(ICL_temp_M$aSMA))
max (as.numeric(ICL_temp_M$aSMA))
min (as.numeric(ICL_temp_M$aSMA))

mean (as.numeric(ICL_temp_F$aSMA))
sd (as.numeric(ICL_temp_F$aSMA))
max (as.numeric(ICL_temp_F$aSMA))
min (as.numeric(ICL_temp_F$aSMA))

wilcox.test(x = as.numeric(ICL_temp_M$aSMA), y = as.numeric(ICL_temp_F$aSMA),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

mean (as.numeric(ICL_temp_M$MHU))
sd (as.numeric(ICL_temp_M$MHU))
max (as.numeric(ICL_temp_M$MHU))
min (as.numeric(ICL_temp_M$MHU))

mean (as.numeric(ICL_temp_F$MHU))
sd (as.numeric(ICL_temp_F$MHU))
max (as.numeric(ICL_temp_F$MHU))
min (as.numeric(ICL_temp_F$MHU))

wilcox.test(x = as.numeric(ICL_temp_M$MHU), y = as.numeric(ICL_temp_F$MHU),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)


#Milan

MI <- Milan[!is.na(Milan$Age),]

wilcox.test(x = as.numeric(Imperial$Sex), y = as.numeric(MI$Sex),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)
  

mean (as.numeric(MI$Age))
sd (as.numeric(MI$Age))
max (as.numeric(MI$Age))
min (as.numeric(MI$Age))

tbl <- table(MI$Sex)
cbind(tbl,prop.table(tbl))

Milan_M <- MI[MI$Sex==1,]
Milan_F <- MI[MI$Sex==2,]

mean (as.numeric(Milan_M$Height))
sd (as.numeric(Milan_M$Height))
max (as.numeric(Milan_M$Height))
min (as.numeric(Milan_M$Height))

mean (as.numeric(Milan_F$Height))
sd (as.numeric(Milan_F$Height))
max (as.numeric(Milan_F$Height))
min (as.numeric(Milan_F$Height))

wilcox.test(x = Milan_M$Height, y = as.numeric(Milan_F$Height),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

mean (as.numeric(Milan_M$BMI))
sd (as.numeric(Milan_M$BMI))
max (as.numeric(Milan_M$BMI))
min (as.numeric(Milan_M$BMI))

mean (as.numeric(Milan_F$BMI))
sd (as.numeric(Milan_F$BMI))
max (as.numeric(Milan_F$BMI))
min (as.numeric(Milan_F$BMI))

wilcox.test(x = Milan_M$BMI, y = as.numeric(Milan_F$BMI),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

ICL_temp_F <- Imperial_F[!is.na(Imperial_F$PS),]
ICL_temp_M <- Imperial_M[!is.na(Imperial_M$PS),]

tbl <- table(Milan_M$`Performance Status`)
cbind(tbl,prop.table(tbl))

tbl <- table(Milan_F$`Performance Status`)
cbind(tbl,prop.table(tbl))

tbl <- table(Milan$`Tumour Stage`)
cbind(tbl,prop.table(tbl))

wilcox.test(x = as.numeric(Milan_M$`Performance Status`), y = as.numeric(Milan_F$`Performance Status`),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

Milan_temp_M <- Milan_M[!is.na(Milan_M$aSMA),]
Milan_temp_F <- Milan_F[!is.na(Milan_F$aSMA),]

mean (as.numeric(Milan_temp_M$SMA))
sd (as.numeric(Milan_temp_M$SMA))
max (as.numeric(Milan_temp_M$SMA))
min (as.numeric(Milan_temp_M$SMA))

Milan_temp_F <- Milan_temp_F[-2,]

mean (as.numeric(Milan_temp_F$SMA))
sd (as.numeric(Milan_temp_F$SMA))
max (as.numeric(Milan_temp_F$SMA))
min (as.numeric(Milan_temp_F$SMA))

wilcox.test(x = as.numeric(Milan_temp_M$SMA), y = as.numeric(Milan_temp_F$SMA),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

mean (as.numeric(Milan_temp_M$aSMA))
sd (as.numeric(Milan_temp_M$aSMA))
max (as.numeric(Milan_temp_M$aSMA))
min (as.numeric(Milan_temp_M$aSMA))

mean (as.numeric(Milan_temp_F$aSMA))
sd (as.numeric(Milan_temp_F$aSMA))
max (as.numeric(Milan_temp_F$aSMA))
min (as.numeric(Milan_temp_F$aSMA))

wilcox.test(x = as.numeric(Milan_temp_M$aSMA), y = as.numeric(Milan_temp_F$aSMA),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

mean (as.numeric(Milan_temp_M$MHU))
sd (as.numeric(Milan_temp_M$MHU))
max (as.numeric(Milan_temp_M$MHU))
min (as.numeric(Milan_temp_M$MHU))

mean (as.numeric(Milan_temp_F$MHU))
sd (as.numeric(Milan_temp_F$MHU))
max (as.numeric(Milan_temp_F$MHU))
min (as.numeric(Milan_temp_F$MHU))

wilcox.test(x = as.numeric(Milan_temp_M$MHU), y = as.numeric(Milan_temp_F$MHU),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)


wilcox.test(x = as.numeric(Milan_temp_M$MHU), y = as.numeric(ICL_temp_M$MHU),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

wilcox.test(x = as.numeric(Milan_temp_F$MHU), y = as.numeric(ICL_temp_F$MHU),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)


wilcox.test(x = as.numeric(Milan_temp_M$aSMA), y = as.numeric(ICL_temp_M$aSMA),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

wilcox.test(x = as.numeric(Milan_temp_F$aSMA), y = as.numeric(ICL_temp_F$aSMA),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

wilcox.test(x = as.numeric(Milan_temp_M$SMA), y = as.numeric(ICL_temp_M$SMA),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

wilcox.test(x = as.numeric(Milan_temp_F$SMA), y = as.numeric(ICL_temp_F$SMA),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

wilcox.test(x = as.numeric(Milan_M$`Performance Status`), y = as.numeric(ICL_temp_M$PS),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

wilcox.test(x = as.numeric(Milan_F$`Performance Status`), y = as.numeric(ICL_temp_F$PS),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

wilcox.test(x = as.numeric(Milan_M$BMI), y = as.numeric(ICL_temp_M$BMI),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

wilcox.test(x = as.numeric(Milan_F$BMI), y = as.numeric(ICL_temp_F$BMI),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

wilcox.test(x = as.numeric(Milan_M$Height), y = as.numeric(ICL_temp_M$Height),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

wilcox.test(x = as.numeric(Milan_F$Height), y = as.numeric(ICL_temp_F$Height),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)


#Survival, trial and disease information

MI <- Milan[!is.na(Milan$`Overall Survival (months)`),]

mean (as.numeric(MI$`Overall Survival (months)`))
sd (as.numeric(MI$`Overall Survival (months)`))
max (as.numeric(MI$`Overall Survival (months)`))
min (as.numeric(MI$`Overall Survival (months)`))

mean (as.numeric(MI$`Progression-free Survival (months)`))
sd (as.numeric(MI$`Progression-free Survival (months)`))
max (as.numeric(MI$`Progression-free Survival (months)`))
min (as.numeric(MI$`Progression-free Survival (months)`))

ICL <- Imperial[!is.na(Imperial$OS),]

mean (as.numeric(ICL$OS/30.42))
sd (as.numeric(ICL$OS/30.42))
max (as.numeric(ICL$OS/30.42))
min (as.numeric(ICL$OS/30.42))

#NAE
ICL <- Imperial[!is.na(Imperial$NAE),]

ICL <- ImperialIO[!is.na(ImperialIO$NAE),]

MI <- Milan[!is.na(Milan$NAE),]

wilcox.test(x = as.numeric(ICL$HighestGrade), y = as.numeric(MI$HighestGrade),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

mean (as.numeric(ICL$NAE))
sd (as.numeric(ICL$NAE))
max (as.numeric(ICL$NAE))
min (as.numeric(ICL$NAE))

tbl <- table(Imperial$NAE)
cbind(tbl,prop.table(tbl))

tbl <- table(ImperialIO$NAE)
cbind(tbl,prop.table(tbl))

tbl <- table(ImperialIO$HighestGrade)
cbind(tbl,prop.table(tbl))

tbl <- table(Milan$HighestGrade)
cbind(tbl,prop.table(tbl))


wilcox.test(x = as.numeric(ICL$OS/30.42), y = as.numeric(MI$`Overall Survival (months)`),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)

ICL <- ImperialIO[!is.na(ImperialIO$OS),]

mean (as.numeric(ICL$OS/30.42))
sd (as.numeric(ICL$OS/30.42))
max (as.numeric(ICL$OS/30.42))
min (as.numeric(ICL$OS/30.42))

wilcox.test(x = as.numeric(ICL$OS/30.42), y = as.numeric(MI$`Overall Survival (months)`),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)


tbl <- table(ICL$`Cancer type`)
cbind(tbl,prop.table(tbl))

tbl <- table(ICL$HighestGrade)
cbind(tbl,prop.table(tbl))

tbl <- table(Milan$HighestGrade)
cbind(tbl,prop.table(tbl))

wilcox.test(x = as.numeric(ICL$HighestGrade), y = as.numeric(Milan$HighestGrade),
            mu=0, alt="two.sided", conf.int=T, conf.level=0.95,
            paired=FALSE, exact=T, correct=T)


tbl <- table(Milan_df$`Trial Drug`)
cbind(tbl,prop.table(tbl))

tbl <- table(ImperialIO_df$`Trial Drug`)
cbind(tbl,prop.table(tbl))

#################Draw Forest Plot for Multivariable Regression##################
#install.packages("xfun")
#install.packages("forestplot")
#install.packages("forestmodel")
library(forestplot)
library(caret)

library(tidyverse)
library(forestmodel)


set.seed(100)

Milan$`Performance Status`

logitmod1 <- coxph(Surv(NewOS, Event) ~ aSMA + MHU + BMI + BSA + Performance Status`, data=Milan)

logitmod1 <- coxph(Surv(NewOS, Event) ~ aSMA + MHU + BMI + BSA, data=Imperial)

logitmod1 <- coxph(Surv(NewOS, Event) ~ MHU, data=Imperial)

logitmod1 <- glm(SevereAE ~ aSMA + MHU + BMI + BSA + PS, data = Milan, family = binomial)

logitmod1 <- glm(SevereAE ~ aSMA + MHU + BMI + BSA, data = Imperial, family = binomial)

summary(logitmod1)
confint(logitmod1)


plot1 <- forest_model(logitmod1)
plot1

#Every pair analysis 
############Investigative Correlation Plots######################################
install.packages("corrplot")
library(corrplot)

install.packages("GGally")
library(GGally)

#Milan - Variables
aSMA <- as.numeric(Milan$aSMA)
MHU <- as.numeric(Milan$MHU)
BMI <- as.numeric(Milan$BMI)
BSA <- as.numeric(Milan$BSA)
Composite <- as.numeric(Milan$Composite)
NAE <- as.numeric(Milan$NAE)

Survival <- as.numeric(Milan$`Overall Survival (months)`)
HighestGrade <- as.numeric(Milan$HighestGrade)

#data <- data.frame(aSMA, MHU, BMI, BSA, Composite, Survival, HighestGrade)

x_vars <- data.frame(aSMA, MHU, BMI, BSA, Composite)
y_vars <- data.frame(NAE, HighestGrade, Survival)


corr_matrix <- cor(x_vars, y_vars, use="pairwise.complete.obs")
# Plot heatmap
corrplot(corr_matrix, method="color", col=colorRampPalette(c("blue", "white", "red"))(200), tl.col="black", tl.cex=0.8, addCoef.col="black", number.cex=0.7)


# Compute correlation matrix
corr_matrix <- cor(x_vars, y_vars, use="pairwise.complete.obs")

# Plot heatmap
corrplot(corr_matrix, method="color", col=colorRampPalette(c("blue", "white", "red"))(200), 
             +          tl.col="black", tl.cex=0.8, addCoef.col="black", number.cex=0.7)



#Imperial - Variables
aSMA <- as.numeric(Imperial $aSMA)
MHU <- as.numeric(Imperial $MHU)
BMI <- as.numeric(Imperial $BMI)
BSA <- as.numeric(Imperial $BSA)
Composite <- as.numeric(Imperial $Composite)
Survival <- as.numeric(Imperial$Event)
HighestGrade <- as.numeric(Imperial$HighestGrade)
HospitalStay <- as.numeric(Imperial$HospitalStay)

data <- data.frame(aSMA, MHU, BMI, BSA, Composite, Survival, HighestGrade, HospitalStay)

ggpairs(data,
        lower = list(continuous = wrap("smooth", method = "lm", se = FALSE)),  # Regression lines
        upper = list(continuous = wrap("cor", size = 5)),  # Correlation coefficients
        diag = list(continuous = wrap("densityDiag")))  # Density plots on diagonal

#ImperialIO - Variables
aSMA <- as.numeric(ImperialIO$aSMA)
MHU <- as.numeric(ImperialIO$MHU)
BMI <- as.numeric(ImperialIO$BMI)
BSA <- as.numeric(ImperialIO$BSA)
Composite <- as.numeric(ImperialIO$Composite)
Survival <- as.numeric(ImperialIO$Event)
HighestGrade <- as.numeric(ImperialIO$HighestGrade)
HospitalStay <- as.numeric(ImperialIO$HospitalStay)

data <- data.frame(aSMA, MHU, BMI, BSA, Composite, Survival, HighestGrade, HospitalStay)

ggpairs(data,
        lower = list(continuous = wrap("smooth", method = "lm", se = FALSE)),  # Regression lines
        upper = list(continuous = wrap("cor", size = 5)),  # Correlation coefficients
        diag = list(continuous = wrap("densityDiag")))  # Density plots on diagonal


###############As a heatmap####################################
install.packages("corrplot")
library(corrplot)

# Compute correlation matrix
data <- na.omit(data)
corr_matrix <- cor(data)

# Plot heatmap
corrplot(corr_matrix, method="color", col=colorRampPalette(c("blue", "white", "red"))(200), 
         tl.col="black", tl.cex=0.8, addCoef.col="black", number.cex=0.7)

#####Milan - Survival and SAE prediction#############

Milan_old <- Milan

Milan <- MilanM

Milan <- Milan[!is.na(Milan$Event),]

Milan <- Milan[!is.na(Milan$SMA),]

#Survival prediction
cate_dis <- as.numeric(Milan$Event)

#AE prediction
Milan <- Milan[!is.na(Milan$SevereAE),]
Milan <- Milan[!is.na(Milan$SMA),]

cate_dis <- as.numeric(Milan$SevereAE)

for (i in 1:length(cate_dis)) {
  if (cate_dis[i] == "1"){
    cate_dis[i] <- '+'
  } else {
    cate_dis[i] <- '-'
  }
}


#SMA
score_dis <- Milan$SMA
#Adjusted SMA to height
score_dis <- Milan$aSMA
#Prado formula for fat-free mass estimate
Milan$pSMA <- Milan$aSMA * 0.30 + 6.06
#MHU
score_dis <- Milan$MHU
#BMI - predictive
score_dis <- Milan$BMI
#Body total surface area - predictive
score_dis <- Milan$BSA
#Composite score
score_dis <- Milan$Composite
#Performance score
score_dis <- Milan$`Performance Status`


ROCit_obj <- rocit(score=score_dis, class=cate_dis, negref = "-", method ="bin")

par(cex.axis=1.0)

AUC_obj <- ciAUC(ROCit_obj, level = 0.95)
p <- plot(ROCit_obj)
text(0.8, 0.4, paste0("AUC=", round(AUC_obj$AUC, 2), ", 95% CI [", round(AUC_obj$lower, 2), ",", round(AUC_obj$upper, 2), "]"), adj = 1, font = 4, cex=1.0)
title("SAE Prediction - Milan")
title("MHU for 90D OS - Milan")
title("aSMA for 3Y OS - Milan")
title("Composite1Y for 1Y OS Event - Milan")

########Imperial - Survival and SAE Prediction########### 

Imperial_old <- Imperial

Imperial_old <- ImperialIO

Imperial <- ImperialIO

Imperial <- ImperialIOF


Imperial <- Imperial[!is.na(Imperial$aSMA),]

cate_dis <- as.numeric(Imperial$Event)

#AE prediction
Imperial <- Imperial[!is.na(Imperial$SevereAE),]
Imperial <- Imperial[!is.na(Imperial$SMA),]
Imperial <- Imperial[!is.na(Imperial$BMI),]
Imperial <- Imperial[!is.na(Imperial$BSA),]
Imperial <- Imperial[!is.na(Imperial$Composite),]

cate_dis <- as.numeric(Imperial$SevereAE)

cate_dis <- as.numeric(Imperial$HospitalStayB)

cate_dis <- as.numeric(Imperial$LongStay)

for (i in 1:length(cate_dis)) {
  if (cate_dis[i] == "1"){
    cate_dis[i] <- '+'
  } else {
    cate_dis[i] <- '-'
  }
}

score_dis <- Imperial$SMA
score_dis <- Imperial$aSMA
score_dis <- Imperial$MHU

score_dis <- Imperial$BMI
score_dis <- Imperial$BSA
score_dis <- Imperial$Composite


ROCit_obj <- rocit(score=score_dis, class=cate_dis, negref = "-", method ="bin")

par(cex.axis=1.0)

AUC_obj <- ciAUC(ROCit_obj, level = 0.95)
p <- plot(ROCit_obj)
text(0.8, 0.4, paste0("AUC=", round(AUC_obj$AUC, 2), ", 95% CI [", round(AUC_obj$lower, 2), ",", round(AUC_obj$upper, 2), "]"), adj = 1, font = 4, cex=1.0)
title("MHU for 90d Event - ImperialIO")
title("aSMA for 1Y OS - Imperial")
title("MHU for 1Y OS - Imperial")
title("aSMA for 90d Event - Imperial")
title("Composite1Y for 1Y OS Event - ImperialIO")

#######C-index#####################

library(survival)
#install.packages("survAUC")
library(survAUC)

# Compute C-index 
Milan$NewOS<- as.numeric(Milan$NewOS)
Milan$Event <- as.numeric(Milan$Event)

cox_model <- coxph(Surv(NewOS, Event) ~ SMA, data = Milan)
cox_model <- coxph(Surv(NewOS, Event) ~ aSMA, data = Milan)
cox_model <- coxph(Surv(NewOS, Event) ~ MHU, data = Milan)
cox_model <- coxph(Surv(NewOS, Event) ~ BMI, data = Milan)
cox_model <- coxph(Surv(NewOS, Event) ~ BSA, data = Milan)
cox_model <- coxph(Surv(NewOS, Event) ~ `Performance Status`, data = Milan)

# Extract predicted risk scores
risk_scores <- predict(cox_model, type = "lp")  # Linear predictor

c_index <- summary(cox_model)$concordance[1]  # Concordance index
se <- summary(cox_model)$concordance[2]       # Standard error

# Compute 95% Confidence Interval
ci_lower <- c_index - 1.96 * se
ci_upper <- c_index + 1.96 * se

# Print results
print(paste("C-index:", round(c_index, 3)))
print(paste("95% CI:", round(ci_lower, 3), "-", round(ci_upper, 3)))

#Imperial

Imperial$OS <- as.numeric(Imperial$OS)
Imperial$Event <- as.numeric(Milan$Event)

cox_model <- coxph(Surv(NewOS, Event) ~ SMA, data = Imperial)
cox_model <- coxph(Surv(NewOS, Event) ~ aSMA, data = Imperial)
cox_model <- coxph(Surv(NewOS, Event) ~ MHU, data = Imperial)
cox_model <- coxph(Surv(NewOS, Event) ~ BMI, data = Imperial)
cox_model <- coxph(Surv(NewOS, Event) ~ BSA, data = Imperial)

# Extract predicted risk scores
risk_scores <- predict(cox_model, type = "lp")  # Linear predictor

c_index <- summary(cox_model)$concordance[1]  # Concordance index
se <- summary(cox_model)$concordance[2]       # Standard error

# Compute 95% Confidence Interval
ci_lower <- c_index - 1.96 * se
ci_upper <- c_index + 1.96 * se

# Print results
print(paste("C-index:", round(c_index, 3)))
print(paste("95% CI:", round(ci_lower, 3), "-", round(ci_upper, 3)))

########Survival Plots################################

library(survival)
library(survminer)

install.packages("glue")
install.packages("dplyr")
library(dplyr)

#################Milan - Survival Curves################
# Automatically find the best cutoff for Milan
cut <- surv_cutpoint(Milan, time = "NewOS", event = "Event", variables = "MHU")

# Apply the cutoff to create two groups
Milan$MHU_group <- 0
Milan$MHU_group <- ifelse(Milan$MHU > cut$cutpoint$cutpoint, "High MHU", "Low MHU")

# Fit Kaplan-Meier survival curves
fit <- survfit(Surv(NewOS, Event) ~ MHU_group, data = Milan)

#Try aSMA
cut <- surv_cutpoint(Milan, time = "NewOS", event = "Event", variables = "BMI")

# Apply the cutoff to create two groups
Milan$MHU_group <- 0
Milan$MHU_group <- ifelse(Milan$BMI> cut$cutpoint$cutpoint, "High BMI", "Low BMI")

#By composite score
cut <- surv_cutpoint(Milan, time = "NewOS", event = "Event", variables = "Composite")

# Apply the cutoff to create two groups
Milan$MHU_group <- 0
Milan$MHU_group <- ifelse(Milan$Composite> cut$cutpoint$cutpoint, "High Composite Score", "Low Composite Score")

# Fit Kaplan-Meier survival curves
fit <- survfit(Surv(NewOS, Event) ~ MHU_group, data = Milan)

# Plot survival curves
ggsurvplot(fit, data = Milan, pval = TRUE, risk.table = TRUE)

#Fancy KM Plot
a<-ggsurvplot(
  survfit(Surv(NewOS,Event) ~ MHU_group, data= Milan),     # survfit object with calculated statistics.
  data = Milan,               # data used to fit survival curves. 
  pval = T,             # show p-value of log-rank test.
  conf.int = T,         # show confidence intervals for point estimates of survival curves.
  xlim = c(0,90),        # present narrower X axis, but not affect survival estimates.
  break.time.by = 10,     # break X axis in time intervals by 100.
  ggtheme = theme,         # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T, # show bars instead of names in text annotations in legend of risk table
  legend.labs=c("Low Risk", "High Risk"),
  risk.table = T,
  palette="jco",
  tables.theme = theme,
  #title = "Kaplan-Meier Plots of Stratified Groups By MHU - Milan",
  title = "Kaplan-Meier Plots of Stratified Groups By Composite Score - Milan",
  xlab="Time in Days",
  ylab="Probability of Overall Survival",
  surv.median.line = "v",
  ylim=c(0,1),
  cumevents=F,
  surv.scale="percent"
)

# extract ggplot object from ggsurvplot
p <- a$plot 
p <- p + scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 90))

# extract table object from ggsurvplot
tab <- a$table
tab$layers = NULL # clear labels
tab <- tab + 
  geom_text(aes(x = time, y = rev(strata), label = llabels), data = tab$data[tab$data$time %in% c(0, 20, 40, 60, 80, 90),]) +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 90))

# extract cumevents object from ggsurvplot
tab2 <- a$cumevents
tab2$layers = NULL # clear labels
tab2 <- tab2 + 
  geom_text(aes(x = time, y = rev(strata), label = cum.n.event), data = tab$data[tab$data$time %in% c(0, 20, 40, 60, 80, 90),]) +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 90))

# Add plots back
a$plot <- p
a$table <- tab
a$cumevents <- tab2

a

coxph(Surv(NewOS,Event) ~ MHU_group, data= Milan) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

theme <- theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_line(colour = "white"),
               panel.grid.minor = element_line(colour = "white"),
               panel.border = element_blank(),
               panel.background = element_blank()) 

#Fancy KM Plot
a<-ggsurvplot(
  survfit(Surv(NewOS,Event) ~ MHU_group, data= Milan),     # survfit object with calculated statistics.
  data = Milan,               # data used to fit survival curves. 
  pval = T,             # show p-value of log-rank test.
  conf.int = T,         # show confidence intervals for point estimates of survival curves.
  xlim = c(0,1095),        # present narrower X axis, but not affect survival estimates.
  break.time.by = 100,     # break X axis in time intervals by 100.
  ggtheme = theme,         # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T, # show bars instead of names in text annotations in legend of risk table
  legend.labs=c("Low Risk", "High Risk"),
  risk.table = T,
  palette="jco",
  tables.theme = theme,
  title = "Kaplan-Meier Plots of Stratified Groups By MHU - Milan",
  xlab="Time in Days",
  ylab="Probability of Overall Survival",
  surv.median.line = "v",
  ylim=c(0,1),
  cumevents=F,
  surv.scale="percent"
)

# extract ggplot object from ggsurvplot
p <- a$plot 
p <- p + scale_x_continuous(breaks = c(0, 200, 400, 600, 800, 1000, 1095))

# extract table object from ggsurvplot
tab <- a$table
tab$layers = NULL # clear labels
tab <- tab + 
  geom_text(aes(x = time, y = rev(strata), label = llabels), data = tab$data[tab$data$time %in% c(0, 200, 400, 600, 800, 1000),]) +
  scale_x_continuous(breaks = c(0, 200, 400, 600, 800, 1000))

# extract cumevents object from ggsurvplot
tab2 <- a$cumevents
tab2$layers = NULL # clear labels
tab2 <- tab2 + 
  geom_text(aes(x = time, y = rev(strata), label = cum.n.event), data = tab$data[tab$data$time %in% c(0, 200, 400, 600, 800, 1000),]) +
  scale_x_continuous(breaks = c(0, 200, 400, 600, 800, 1000))

# Add plots back
a$plot <- p
a$table <- tab
a$cumevents <- tab2

a

###############Imperial - Survival Curves#######################################

# Automatically find the best cutoff for Imperial - SMA Female
cut <- surv_cutpoint(Imperial, time = "NewOS", event = "Event", variables = "aSMA")

# Apply the cutoff to create two groups
Imperial$SMA_group <- 0
Imperial$SMA_group <- ifelse(Imperial$aSMA > cut$cutpoint$cutpoint, "High SMA", "Low SMA")

# Fit Kaplan-Meier survival curves
fit <- survfit(Surv(NewOS, Event) ~ SMA_group, data = Imperial)

# Plot survival curves
ggsurvplot(fit, data = Imperial, pval = TRUE, risk.table = TRUE)

coxph(Surv(NewOS,Event) ~ SMA_group, data= Imperial) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

theme <- theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_line(colour = "white"),
               panel.grid.minor = element_line(colour = "white"),
               panel.border = element_blank(),
               panel.background = element_blank()) 

#Fancy KM Plot
a<-ggsurvplot(
  survfit(Surv(NewOS,Event) ~ MHU_group, data= Imperial),     # survfit object with calculated statistics.
  data = Imperial,               # data used to fit survival curves. 
  pval = T,             # show p-value of log-rank test.
  conf.int = T,         # show confidence intervals for point estimates of survival curves.
  xlim = c(0,90),        # present narrower X axis, but not affect survival estimates.
  break.time.by = 10,     # break X axis in time intervals by 100.
  ggtheme = theme,         # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T, # show bars instead of names in text annotations in legend of risk table
  legend.labs=c("Low Risk", "High Risk"),
  risk.table = T,
  palette="jco",
  tables.theme = theme,
  title = "Kaplan-Meier Plots of Stratified Groups By MHU - Imperial",
  xlab="Time in Days",
  ylab="Probability of Overall Survival",
  surv.median.line = "v",
  ylim=c(0,1),
  cumevents=F,
  surv.scale="percent"
)

# extract ggplot object from ggsurvplot
p <- a$plot 
p <- p + scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 90))

# extract table object from ggsurvplot
tab <- a$table
tab$layers = NULL # clear labels
tab <- tab + 
  geom_text(aes(x = time, y = rev(strata), label = llabels), data = tab$data[tab$data$time %in% c(0, 20, 40, 60, 80, 90),]) +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 90))

# extract cumevents object from ggsurvplot
tab2 <- a$cumevents
tab2$layers = NULL # clear labels
tab2 <- tab2 + 
  geom_text(aes(x = time, y = rev(strata), label = cum.n.event), data = tab$data[tab$data$time %in% c(0, 20, 40, 60, 80, 90),]) +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 90))

# Add plots back
a$plot <- p
a$table <- tab
a$cumevents <- tab2

a

# Automatically find the best cutoff for Imperial
cut <- surv_cutpoint(Imperial, time = "NewOS", event = "Event", variables = "MHU")

# Apply the cutoff to create two groups
Imperial$MHU_group <- 0
Imperial$MHU_group <- ifelse(Imperial$MHU > cut$cutpoint$cutpoint, "High MHU", "Low MHU")

# Fit Kaplan-Meier survival curves
fit <- survfit(Surv(NewOS, Event) ~ MHU_group, data = Imperial)

# Plot survival curves
ggsurvplot(fit, data = Imperial, pval = TRUE, risk.table = TRUE)

coxph(Surv(NewOS,Event) ~ MHU_group, data= Imperial) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

theme <- theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_line(colour = "white"),
               panel.grid.minor = element_line(colour = "white"),
               panel.border = element_blank(),
               panel.background = element_blank()) 

#Fancy KM Plot
a<-ggsurvplot(
  survfit(Surv(NewOS,Event) ~ MHU_group, data= Imperial),     # survfit object with calculated statistics.
  data = Imperial,               # data used to fit survival curves. 
  pval = T,             # show p-value of log-rank test.
  conf.int = T,         # show confidence intervals for point estimates of survival curves.
  xlim = c(0,90),        # present narrower X axis, but not affect survival estimates.
  break.time.by = 10,     # break X axis in time intervals by 100.
  ggtheme = theme,         # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T, # show bars instead of names in text annotations in legend of risk table
  legend.labs=c("Low Risk", "High Risk"),
  risk.table = T,
  palette="jco",
  tables.theme = theme,
  title = "Kaplan-Meier Plots of Stratified Groups By MHU - ImperialIO",
  xlab="Time in Days",
  ylab="Probability of Overall Survival",
  surv.median.line = "v",
  ylim=c(0,1),
  cumevents=F,
  surv.scale="percent"
)

# extract ggplot object from ggsurvplot
p <- a$plot 
p <- p + scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 90))

# extract table object from ggsurvplot
tab <- a$table
tab$layers = NULL # clear labels
tab <- tab + 
  geom_text(aes(x = time, y = rev(strata), label = llabels), data = tab$data[tab$data$time %in% c(0, 20, 40, 60, 80, 90),]) +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 90))

# extract cumevents object from ggsurvplot
tab2 <- a$cumevents
tab2$layers = NULL # clear labels
tab2 <- tab2 + 
  geom_text(aes(x = time, y = rev(strata), label = cum.n.event), data = tab$data[tab$data$time %in% c(0, 20, 40, 60, 80, 90),]) +
  scale_x_continuous(breaks = c(0, 20, 40, 60, 80, 90))

# Add plots back
a$plot <- p
a$table <- tab
a$cumevents <- tab2

a

theme <- theme(axis.line = element_line(colour = "black"),
               panel.grid.major = element_line(colour = "white"),
               panel.grid.minor = element_line(colour = "white"),
               panel.border = element_blank(),
               panel.background = element_blank()) 

#Fancy KM Plot
a<-ggsurvplot(
  survfit(Surv(NewOS,Event) ~ MHU_group, data= Imperial),     # survfit object with calculated statistics.
  data = Imperial,               # data used to fit survival curves. 
  pval = T,             # show p-value of log-rank test.
  conf.int = T,         # show confidence intervals for point estimates of survival curves.
  xlim = c(0,365),        # present narrower X axis, but not affect survival estimates.
  break.time.by = 50,     # break X axis in time intervals by 100.
  ggtheme = theme,         # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = T, # show bars instead of names in text annotations in legend of risk table
  legend.labs=c("Low Risk", "High Risk"),
  risk.table = T,
  palette="jco",
  tables.theme = theme,
  title = "Kaplan-Meier Plots of Stratified Groups By MHU - ImperialIO",
  xlab="Time in Days",
  ylab="Probability of Overall Survival",
  surv.median.line = "v",
  ylim=c(0,1),
  cumevents=F,
  surv.scale="percent"
)

# extract ggplot object from ggsurvplot
p <- a$plot 
p <- p + scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300, 350, 365))

# extract table object from ggsurvplot
tab <- a$table
tab$layers = NULL # clear labels
tab <- tab + 
  geom_text(aes(x = time, y = rev(strata), label = llabels), data = tab$data[tab$data$time %in% c(0, 50, 100, 150, 200, 250, 300, 350, 365),]) +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300, 350, 365))

# extract cumevents object from ggsurvplot
tab2 <- a$cumevents
tab2$layers = NULL # clear labels
tab2 <- tab2 + 
  geom_text(aes(x = time, y = rev(strata), label = cum.n.event), data = tab$data[tab$data$time %in% c(0, 50, 100, 150, 200, 250, 300, 350, 365),]) +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300, 350, 365))

# Add plots back
a$plot <- p
a$table <- tab
a$cumevents <- tab2

a

##############Nomogram Creation#######################
#install.packages(c("PASWR","rms"))
library(PASWR)
library("rms")

T_roc <- Dis_combined[ -c(10,11,12) ]
T_roc <- subset(T_roc, T_roc$Histology<3)

T_roc <- T_roc[complete.cases(T_roc$His),]

T_roc$`EGFR-RPV` <- abs(T_roc$s0)
t.data <- datadist(T_roc)
options(datadist = 't.data')

fit <- lrm(formula = His ~ Gender + Ethnicity + `EGFR-RPV`, data = T_roc)
plot(nomogram(fit, fun = function(x)plogis(x)))


#########Severity of toxicity (grade 3 or 4) Prediction############

cate_dis <- as.numeric(Milan$SevereAE)

for (i in 1:length(cate_dis)) {
  if (cate_dis[i] == "1"){
    cate_dis[i] <- '+'
  } else {
    cate_dis[i] <- '-'
  }
}

score_dis <- Milan$aSMA
score_dis <- Milan$MHU
score_dis <- Milan$SMAPrado

ROCit_obj <- rocit(score=score_dis, class=cate_dis, negref = "-", method ="bin")

par(cex.axis=1.0)

AUC_obj <- ciAUC(ROCit_obj, level = 0.95)
p <- plot(ROCit_obj)
text(0.95, 0.2, paste0("AUC=", round(AUC_obj$AUC, 2), ", 95% CI [", round(AUC_obj$lower, 2), ",", round(AUC_obj$upper, 2), "]"), adj = 1, font = 4, cex=1.0)
title("SAE - Milan Cohort aSMA")

#Imperial
Imperial <- Imperial[!is.na(Imperial$SevereAE),]

cate_dis <- as.numeric(Imperial$SevereAE)

for (i in 1:length(cate_dis)) {
  if (cate_dis[i] == "1"){
    cate_dis[i] <- '+'
  } else {
    cate_dis[i] <- '-'
  }
}

score_dis <- Imperial$MHU

ROCit_obj <- rocit(score=score_dis, class=cate_dis, negref = "-", method ="bin")

par(cex.axis=1.0)

AUC_obj <- ciAUC(ROCit_obj, level = 0.95)
p <- plot(ROCit_obj)
text(0.95, 0.2, paste0("AUC=", round(AUC_obj$AUC, 2), ", 95% CI [", round(AUC_obj$lower, 2), ",", round(AUC_obj$upper, 2), "]"), adj = 1, font = 4, cex=1.0)
title("SAE - Imperial Cohort SMA")


#################Grouped by system#################
#Cohort 0 

#Rearrange the IrAE groupings into systems (Pulmonary, GI, CVS, Endocrine, Neurological, MSK, Others)
colSums(data.frame(Milan$Nausea, Milan$Vomiting, Milan$Diarrhoea, Milan$Disgeusia, Milan$Colitis) != 0)
colSums(data.frame(Milan$Pneumonitis)!= 0)
colSums(data.frame(Milan$Oedema, Milan$Hypertension) != 0)
colSums(data.frame(Milan$Hyperlipasaemia, Milan$Hypothyroidism, Milan$Hyperthyroidism, Milan$Hypophysitis, Milan$Transaminitis) != 0)
colSums(data.frame(Milan$Neuropathy, Milan$HFS)!=0)
colSums(data.frame(Milan$Rash, Milan$Arthralgia, Milan$HFS, Milan$Myositis)!=0)
colSums(data.frame(Milan$Fatigue, Milan$Chelitis, Milan$Stomatitis, Milan$Mucositis, Milan$Proteinuria)!= 0)

sum(colSums(data.frame(Milan$Nausea, Milan$Vomiting, Milan$Diarrhoea, Milan$Disgeusia, Milan$Colitis)!=0))
sum(colSums(data.frame(Milan$Pneumonitis)!=0))
sum(colSums(data.frame(Milan$Oedema, Milan$Hypertension)!=0))
sum(colSums(data.frame(Milan$Hyperlipasaemia, Milan$Hypothyroidism, Milan$Hyperthyroidism, Milan$Hypophysitis, Milan$Transaminitis)!=0))
sum(colSums(data.frame(Milan$Neuropathy, Milan$HFS)!=0))
sum(colSums(data.frame(Milan$Rash, Milan$Arthralgia, Milan$HFS, Milan$Myositis)!=0))
sum(colSums(data.frame(Milan$Proteinuria)!=0))
sum(colSums(data.frame(Milan$Fatigue, Milan$Chelitis, Milan$Stomatitis, Milan$Mucositis)))

#GI
sum(colSums(data.frame(Imperial$Proctitis,	Imperial$Abdopain, Imperial$Nausea, Imperial$GORD, Imperial$Dysphagia, Imperial$Anorexia, Imperial$Vomiting, Imperial$Diarrhoea, Imperial$Constipation, Imperial$Colitis, Imperial$Abdopain, Imperial$Disgeusia, Imperial$SBO, Imperial$UGIBleeding, Imperial$LGIBleeding, Imperial$Ascites)!=0))
#HPB
sum(colSums(data.frame(Imperial$Oedema,	Imperial$Hyperlipasaemia,	Imperial$Hypoalbuminaemia,	Imperial$Bilirubinrise,	Imperial$Transaminitis)!=0))
#Pulmonary
sum(colSums(data.frame(Imperial$Pneumonitis, Imperial$SOB, 	Imperial$Effusion, Imperial$Cough, Imperial$Lunginfection,	Imperial$Pneumonitis)!=0))
#CVS
sum(colSums(data.frame(Imperial$Oedema, Imperial$Hypertension, Imperial$Palpitation,	Imperial$AF, Imperial$Hypotension,	Imperial$Hypertension)!=0))
#Endocrine
sum(colSums(data.frame(Imperial$Hyperlipasaemia, Imperial$Hypothyroidism, Imperial$Hyperthyroidism, Imperial$Hypophysitis, Imperial$Hyperglycaemia, Imperial$Hotflush)!=0))
#Neurological
sum(colSums(data.frame(  Imperial$Insomnia,	Imperial$Seizure,	Imperial$Syncope,	Imperial$Vertigo,	Imperial$Cordcompression,	Imperial$Dizziness,	Imperial$Neuropathy, Imperial$Neuropathy, Imperial$Handfootsyndrome, Imperial$Headache, Imperial$Neckpain)!=0))
#SkinJoints
sum(colSums(data.frame(Imperial$Ulweakness, Imperial$Chestwallpain,	Imperial$Pelvispain,	Imperial$Backpain, Imperial$Rash, Imperial$Arthralgia, Imperial$Myositis, Imperial$Injectionsite, Imperial$Skinulceration)!=0))
#Renal
sum(colSums(data.frame(Imperial$Haematuria,	Imperial$Proteinuria,	Imperial$UrinaryRetention,	Imperial$UTI,	Imperial$Electrolytedisturbance,	Imperial$Creatininerise)!=0))    
#Haem
sum(colSums(data.frame(Imperial$Lymphocytonpenia,	Imperial$Neutropenia, Imperial$Thrombocytopenia,	Imperial$Anaemia,	Imperial$INRrise,	Imperial$Thromboticevent)!=0))
#Others
sum(colSums(data.frame(Imperial$Fever,	Imperial$Sepsis,	Imperial$Cytokinerelease,Imperial$Fatigue, Imperial$Insomnia, Imperial$Chelitis, Imperial$Stomatitis, Imperial$Mucositis, Imperial$Proteinuria, Imperial$Vaginalbleed, Imperial$Anaphylaxis,	Imperial$Coryza, Imperial$Pruritis, Imperial$Nospecificpain,	Imperial$Lymphnodepain)!=0)) 


#GI
sum(colSums(data.frame(ImperialIO$Proctitis,	ImperialIO$Abdopain, ImperialIO$Nausea, ImperialIO$GORD, ImperialIO$Dysphagia, ImperialIO$Anorexia, ImperialIO$Vomiting, ImperialIO$Diarrhoea, ImperialIO$Constipation, ImperialIO$Colitis, ImperialIO$Abdopain, ImperialIO$Disgeusia, ImperialIO$SBO, ImperialIO$UGIBleeding, ImperialIO$LGIBleeding, ImperialIO$Ascites)!=0))
#HPB
sum(colSums(data.frame(ImperialIO$Oedema,	ImperialIO$Hyperlipasaemia,	ImperialIO$Hypoalbuminaemia,	ImperialIO$Bilirubinrise,	ImperialIO$Transaminitis)!=0))
#Pulmonary
sum(colSums(data.frame(ImperialIO$Pneumonitis, ImperialIO$SOB, 	ImperialIO$Effusion, ImperialIO$Cough, ImperialIO$Lunginfection,	ImperialIO$Pneumonitis)!=0))
#CVS
sum(colSums(data.frame(ImperialIO$Oedema, ImperialIO$Hypertension, ImperialIO$Palpitation,	ImperialIO$AF, ImperialIO$Hypotension,	ImperialIO$Hypertension)!=0))
#Endocrine
sum(colSums(data.frame(ImperialIO$Hyperlipasaemia, ImperialIO$Hypothyroidism, ImperialIO$Hyperthyroidism, ImperialIO$Hypophysitis, ImperialIO$Hyperglycaemia, ImperialIO$Hotflush)!=0))
#Neurological
sum(colSums(data.frame(  ImperialIO$Insomnia,	ImperialIO$Seizure,	ImperialIO$Syncope,	ImperialIO$Vertigo,	ImperialIO$Cordcompression,	ImperialIO$Dizziness,	ImperialIO$Neuropathy, ImperialIO$Neuropathy, ImperialIO$Handfootsyndrome, ImperialIO$Headache, ImperialIO$Neckpain)!=0))
#SkinJoints
sum(colSums(data.frame(ImperialIO$Ulweakness, ImperialIO$Chestwallpain,	ImperialIO$Pelvispain,	ImperialIO$Backpain, ImperialIO$Rash, ImperialIO$Arthralgia, ImperialIO$Myositis, ImperialIO$Injectionsite, ImperialIO$Skinulceration)!=0))
#Renal
sum(colSums(data.frame(ImperialIO$Haematuria,	ImperialIO$Proteinuria,	ImperialIO$UrinaryRetention,	ImperialIO$UTI,	ImperialIO$Electrolytedisturbance,	ImperialIO$Creatininerise)!=0))    
#Haem
sum(colSums(data.frame(ImperialIO$Lymphocytonpenia,	ImperialIO$Neutropenia, ImperialIO$Thrombocytopenia,	ImperialIO$Anaemia,	ImperialIO$INRrise,	ImperialIO$Thromboticevent)!=0))
#Others
sum(colSums(data.frame(ImperialIO$Fever,	ImperialIO$Sepsis,	ImperialIO$Cytokinerelease,ImperialIO$Fatigue, ImperialIO$Insomnia, ImperialIO$Chelitis, ImperialIO$Stomatitis, ImperialIO$Mucositis, ImperialIO$Proteinuria, ImperialIO$Vaginalbleed, ImperialIO$Anaphylaxis,	ImperialIO$Coryza, ImperialIO$Pruritis, ImperialIO$Nospecificpain,	ImperialIO$Lymphnodepain)!=0))


#Severity
for (x in 1:nrow(Milan)) {
  Milan$GI[x] <- max(Milan$Nausea[x], Milan$Vomiting[x], Milan$Diarrhoea[x], Milan$Disgeusia[x], Milan$Colitis[x])
  Milan$Pulmonary[x] <- max(Milan$Pneumonitis[x])
  Milan$CVS[x] <- max(Milan$Oedema[x], Milan$Hypertension[x])
  Milan$Endocrine[x] <- max(Milan$Hyperlipasaemia[x], Milan$Hypothyroidism[x], Milan$Hyperthyroidism[x], Milan$Hypophysitis[x], Milan$Transaminitis[x])
  Milan$Neurological[x] <- max(Milan$Neuropathy[x], Milan$HFS[x])
  Milan$SkinJoints[x] <- max(Milan$Rash[x], Milan$Arthralgia[x], Milan$HFS[x], Milan$Myositis[x])
  Milan$Others[x] <- max(Milan$Fatigue[x], Milan$Chelitis[x], Milan$Stomatitis[x], Milan$Mucositis[x])
  Milan$Renal[x] <- max(Milan$Proteinuria[x])
}

tbl <- table(Milan$GI)
tbl <- table(Milan$Pulmonary)
tbl <- table(Milan$CVS)
tbl <- table(Milan$Renal)
tbl <- table(Milan$Endocrine)
tbl <- table(Milan$Neurological)
tbl <- table(Milan$SkinJoints)
tbl <- table(Milan$Others)
cbind(tbl,prop.table(tbl))

Milan$GIb <- as.numeric(Milan$GI > 0)
Milan$Pulmonaryb <- as.numeric(Milan$Pulmonary > 0)
Milan$CVSb <- as.numeric(Milan$CVS > 0)
Milan$Endocrineb <- as.numeric(Milan$Endocrine > 0)
Milan$Neurologicalb <- as.numeric(Milan$Neurological > 0)
Milan$SkinJointsb <- as.numeric(Milan$SkinJoints > 0)
Milan$Othersb <- as.numeric(Milan$Others > 0)

for (x in 1:nrow(Imperial)) {
  Imperial$GI[x] <- max(Imperial$Proctitis[x],	Imperial$Abdopain[x], Imperial$Nausea[x], Imperial$GORD[x], Imperial$Dysphagia[x], Imperial$Anorexia[x], Imperial$Vomiting[x], Imperial$Diarrhoea[x], Imperial$Constipation[x], Imperial$Colitis[x], Imperial$Abdopain[x], Imperial$Disgeusia[x], Imperial$SBO[x], Imperial$UGIBleeding[x], Imperial$LGIBleeding[x], Imperial$Ascites[x])
  Imperial$HPB[x] <- max(Imperial$Oedema[x],	Imperial$Hyperlipasaemia[x],	Imperial$Hypoalbuminaemia[x],	Imperial$Bilirubinrise[x],	Imperial$Transaminitis[x])
  Imperial$Pulmonary[x] <- max(Imperial$Pneumonitis[x], Imperial$SOB[x], 	Imperial$Effusion[x], Imperial$Cough[x], Imperial$Lunginfection[x],	Imperial$Pneumonitis[x])
  Imperial$CVS[x] <- max(Imperial$Oedema[x], Imperial$Hypertension[x], Imperial$Palpitation[x],	Imperial$AF[x], Imperial$Hypotension[x],	Imperial$Hypertension[x])
  Imperial$Endocrine[x] <- max(Imperial$Hyperlipasaemia[x], Imperial$Hypothyroidism[x], Imperial$Hyperthyroidism[x], Imperial$Hypophysitis[x], Imperial$Hyperglycaemia[x], 	Imperial$Hotflush[x])
  Imperial$Neurological[x] <- max(Imperial$Insomnia[x],	Imperial$Seizure[x],	Imperial$Syncope[x],	Imperial$Vertigo[x],	Imperial$Cordcompression[x],	Imperial$Dizziness[x],	Imperial$Neuropathy[x], Imperial$Neuropathy[x], Imperial$Handfootsyndrome[x], Imperial$Headache[x], Imperial$Neckpain[x])
  Imperial$SkinJoints[x] <- max(Imperial$Ulweakness[x], Imperial$Chestwallpain[x],	Imperial$Pelvispain[x],	Imperial$Backpain[x], Imperial$Rash[x], Imperial$Arthralgia[x], Imperial$Myositis[x], Imperial$Injectionsite[x], Imperial$Skinulceration[x])
  Imperial$Renal[x] <- max(Imperial$Haematuria[x],	Imperial$Proteinuria[x],	Imperial$UrinaryRetention[x],	Imperial$UTI[x],	Imperial$Electrolytedisturbance[x],	Imperial$Creatininerise[x])
  Imperial$Haem[x] <- max(Imperial$Lymphocytonpenia[x],	Imperial$Neutropenia[x], Imperial$Thrombocytopenia[x],	Imperial$Anaemia[x],	Imperial$INRrise[x],	Imperial$Thromboticevent[x])
  Imperial$Others[x] <- max(Imperial$Fever[x],	Imperial$Sepsis[x],	Imperial$Cytokinerelease[x],Imperial$Fatigue[x], Imperial$Insomnia[x], Imperial$Chelitis[x], Imperial$Stomatitis[x], Imperial$Mucositis[x], Imperial$Proteinuria[x], Imperial$Vaginalbleed[x], Imperial$Anaphylaxis[x],	Imperial$Coryza[x], Imperial$Pruritis[x], Imperial$Nospecificpain[x],	Imperial$Lymphnodepain[x])
}

tbl <- table(Imperial$GI)
tbl <- table(Imperial$HPB)
tbl <- table(Imperial$Pulmonary)
tbl <- table(Imperial$CVS)
tbl <- table(Imperial$Renal)
tbl <- table(Imperial$Endocrine)
tbl <- table(Imperial$Neurological)
tbl <- table(Imperial$SkinJoints)
tbl <- table(Imperial$Haem)
tbl <- table(Imperial$Others)
cbind(tbl,prop.table(tbl))

Imperial$GIb <- as.numeric(Imperial$GI > 0)
Imperial$Pulmonaryb <- as.numeric(Imperial$Pulmonary > 0)
Imperial$CVSb <- as.numeric(Imperial$CVS > 0)
Imperial$Endocrineb <- as.numeric(Imperial$Endocrine > 0)
Imperial$Neurologicalb <- as.numeric(Imperial$Neurological > 0)
Imperial$SkinJointsb <- as.numeric(Imperial$SkinJoints > 0)
Imperial$Othersb <- as.numeric(Imperial$Others > 0)


for (x in 1:nrow(ImperialIO)) {
  ImperialIO$GI[x] <- max(ImperialIO$Proctitis[x],	ImperialIO$Abdopain[x], ImperialIO$Nausea[x], ImperialIO$GORD[x], ImperialIO$Dysphagia[x], ImperialIO$Anorexia[x], ImperialIO$Vomiting[x], ImperialIO$Diarrhoea[x], ImperialIO$Constipation[x], ImperialIO$Colitis[x], ImperialIO$Abdopain[x], ImperialIO$Disgeusia[x], ImperialIO$SBO[x], ImperialIO$UGIBleeding[x], ImperialIO$LGIBleeding[x], ImperialIO$Ascites[x])
  ImperialIO$HPB[x] <- max(ImperialIO$Oedema[x],	ImperialIO$Hyperlipasaemia[x],	ImperialIO$Hypoalbuminaemia[x],	ImperialIO$Bilirubinrise[x],	ImperialIO$Transaminitis[x])
  ImperialIO$Pulmonary[x] <- max(ImperialIO$Pneumonitis[x], ImperialIO$SOB[x], 	ImperialIO$Effusion[x], ImperialIO$Cough[x], ImperialIO$Lunginfection[x],	ImperialIO$Pneumonitis[x])
  ImperialIO$CVS[x] <- max(ImperialIO$Oedema[x], ImperialIO$Hypertension[x], ImperialIO$Palpitation[x],	ImperialIO$AF[x], ImperialIO$Hypotension[x],	ImperialIO$Hypertension[x])
  ImperialIO$Endocrine[x] <- max(ImperialIO$Hyperlipasaemia[x], ImperialIO$Hypothyroidism[x], ImperialIO$Hyperthyroidism[x], ImperialIO$Hypophysitis[x], ImperialIO$Hyperglycaemia[x], 	ImperialIO$Hotflush[x])
  ImperialIO$Neurological[x] <- max(ImperialIO$Insomnia[x],	ImperialIO$Seizure[x],	ImperialIO$Syncope[x],	ImperialIO$Vertigo[x],	ImperialIO$Cordcompression[x],	ImperialIO$Dizziness[x],	ImperialIO$Neuropathy[x], ImperialIO$Neuropathy[x], ImperialIO$Handfootsyndrome[x], ImperialIO$Headache[x], ImperialIO$Neckpain[x])
  ImperialIO$SkinJoints[x] <- max(ImperialIO$Ulweakness[x], ImperialIO$Chestwallpain[x],	ImperialIO$Pelvispain[x],	ImperialIO$Backpain[x], ImperialIO$Rash[x], ImperialIO$Arthralgia[x], ImperialIO$Myositis[x], ImperialIO$Injectionsite[x], ImperialIO$Skinulceration[x])
  ImperialIO$Renal[x] <- max(ImperialIO$Haematuria[x],	ImperialIO$Proteinuria[x],	ImperialIO$UrinaryRetention[x],	ImperialIO$UTI[x],	ImperialIO$Electrolytedisturbance[x],	ImperialIO$Creatininerise[x])
  ImperialIO$Haem[x] <- max(ImperialIO$Lymphocytonpenia[x],	ImperialIO$Neutropenia[x], ImperialIO$Thrombocytopenia[x],	ImperialIO$Anaemia[x],	ImperialIO$INRrise[x],	ImperialIO$Thromboticevent[x])
  ImperialIO$Others[x] <- max(ImperialIO$Fever[x],	ImperialIO$Sepsis[x],	ImperialIO$Cytokinerelease[x],ImperialIO$Fatigue[x], ImperialIO$Insomnia[x], ImperialIO$Chelitis[x], ImperialIO$Stomatitis[x], ImperialIO$Mucositis[x], ImperialIO$Proteinuria[x], ImperialIO$Vaginalbleed[x], ImperialIO$Anaphylaxis[x],	ImperialIO$Coryza[x], ImperialIO$Pruritis[x], ImperialIO$Nospecificpain[x],	ImperialIO$Lymphnodepain[x])
}

tbl <- table(ImperialIO$GI)
tbl <- table(ImperialIO$HPB)
tbl <- table(ImperialIO$Pulmonary)
tbl <- table(ImperialIO$CVS)
tbl <- table(ImperialIO$Renal)
tbl <- table(ImperialIO$Endocrine)
tbl <- table(ImperialIO$Neurological)
tbl <- table(ImperialIO$SkinJoints)
tbl <- table(ImperialIO$Haem)
tbl <- table(ImperialIO$Others)
cbind(tbl,prop.table(tbl))

#########Correlation: AE by System with Body Metrics and PS############ 

install.packages("corrplot")
library(corrplot)

install.packages("GGally")
library(GGally)

#Milan 
aSMA <- as.numeric(Milan$aSMA)
MHU <- as.numeric(Milan$MHU)
BMI <- as.numeric(Milan$BMI)
BSA <- as.numeric(Milan$BSA)
Composite <- as.numeric(Milan$Composite)
PS <- as.numeric(Milan$`Performance Status`)
Survival <- as.numeric(Milan$NewOS)
`Highest Grade` <- as.numeric(Milan$HighestGrade)
NAE <- Milan$NAE

GI <- Milan$GI
Pulmonary <- Milan$Pulmonary
CVS <- Milan$CVS
Renal <- Milan$Renal
Endocrine <- Milan$Endocrine
Neurological <- Milan$Neurological
SkinJoints <- Milan$SkinJoints
Others <- Milan$Others

x_vars <- data.frame(aSMA, MHU, BMI, BSA, Composite, PS)
y_vars <- data.frame(NAE, `Highest Grade`, GI, Pulmonary, CVS, Renal, Endocrine, Neurological, SkinJoints, Others)

#Imperial
aSMA <- as.numeric(Imperial$aSMA)
MHU <- as.numeric(Imperial$MHU)
BMI <- as.numeric(Imperial$BMI)
BSA <- as.numeric(Imperial$BSA)
Composite <- as.numeric(Imperial$Composite)
Survival <- as.numeric(Imperial$OS)
`Highest Grade` <- as.numeric(Imperial$HighestGrade)
NAE <- Imperial$NAE
`Hospital Stay` <- as.numeric(Imperial$HospitalStay)

GI <- Imperial$GI
HPB <- Imperial$HPB
Pulmonary <- Imperial$Pulmonary
CVS <- Imperial$CVS
Renal <- Imperial$Renal
Endocrine <- Imperial$Endocrine
Neurological <- Imperial$Neurological
SkinJoints <- Imperial$SkinJoints
Haem <- Imperial$Haem
Others <- Imperial$Others

x_vars <- data.frame(aSMA, MHU, BMI, BSA, Composite)
y_vars <- data.frame(`Hospital Stay`, NAE, `Highest Grade`, GI, HPB, Pulmonary, CVS, Renal, Endocrine, Neurological, SkinJoints, Haem, Others)

#ImperialIO
aSMA <- as.numeric(ImperialIO$aSMA)
MHU <- as.numeric(ImperialIO$MHU)
BMI <- as.numeric(ImperialIO$BMI)
BSA <- as.numeric(ImperialIO$BSA)
Composite <- as.numeric(ImperialIO$Composite)
#Survival <- as.numeric(ImperialIO$OS)
`Highest Grade` <- as.numeric(ImperialIO$HighestGrade)
NAE <- ImperialIO$NAE
`Hospital Stay` <- as.numeric(ImperialIO$HospitalStay)

GI <- ImperialIO$GI
HPB <- ImperialIO$HPB
Pulmonary <- ImperialIO$Pulmonary
CVS <- ImperialIO$CVS
Renal <- ImperialIO$Renal
Endocrine <- ImperialIO$Endocrine
Neurological <- ImperialIO$Neurological
SkinJoints <- ImperialIO$SkinJoints
Haem <- ImperialIO$Haem
Others <- ImperialIO$Others

x_vars <- data.frame(aSMA, MHU, BMI, BSA, Composite)
y_vars <- data.frame(`Hospital Stay`, NAE, `Highest Grade`, GI, HPB, Pulmonary, CVS, Renal, Endocrine, Neurological, SkinJoints, Haem, Others)

# Compute correlation matrix
corr_matrix <- cor(x_vars, y_vars, use="pairwise.complete.obs")

# Plot heatmap
corrplot(corr_matrix, method="color", col=colorRampPalette(c("blue", "white", "red"))(200), 
         tl.col="black", tl.cex=0.8, addCoef.col="black", number.cex=0.7)

############Prediction targets: 1. Highest severity 2. Incidence of adverse events 3. Stratify by cohort################

#Find the most severe

for (x in 1:nrow(Milan)) {
  Milan$Highest[x] <- max(Milan[x,19:41])
  
  if (Milan$Highest[x] > 0){
  Milan$Incidence[x] <- 1
  }
  else {
    Milan$Incidence[x] <- 0
  }
}
for (x in 1:nrow(Imperial)){
  Imperial$Highest[x] <- max(Imperial[x,14:48])
  
  if (Imperial$Highest[x] > 0){
    Imperial$Incidence[x] <- 1
  }
  else {
    Imperial$Incidence[x] <- 0
  }
}


Age <- Milan$Age
hist(Age)
summary(Age)

mean(Milan$BMI, na.rm = TRUE)
sd(Milan$BMI, na.rm = TRUE)

mean(Milan$aSMA, na.rm = TRUE)
sd(Milan$aSMA, na.rm = TRUE)

table (Milan$Sex)


mean(Imperial$Age, na.rm= TRUE)
sd(Imperial$Age, na.rm= TRUE)

table (Imperial$Sex)



