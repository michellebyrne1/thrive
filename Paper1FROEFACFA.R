library(readr)
DBHSfro2022 <- read_csv("Documents/2022MergedAlldataLongitudinalFirstResponseOnlyNoMissingData.csv")

View(DBHSfro2022)

library(lavaan)
library("lmerTest")
library("lme4")

#random allocated 10% subsample for EFA and 90% for CFA
#DBHSfro2022$subsample <- sample(factor(rep(1:10, length.out=nrow(DBHSfro2022)), 
                                       labels=paste0("subsample", 1:10)))
#HSEFAfro_subsample <- DBHSfro2022[which(DBHSfro2022$subsample=="subsample1"), ]
#HSCFAfro_subsample <- DBHSfro2022[which(DBHSfro2022$subsample=="subsample2"|DBHSfro2022$subsample=="subsample3"|DBHSfro2022$subsample=="subsample4"|DBHSfro2022$subsample=="subsample5"|DBHSfro2022$subsample=="subsample6"|DBHSfro2022$subsample=="subsample7"|DBHSfro2022$subsample=="subsample8"|DBHSfro2022$subsample=="subsample9"|DBHSfro2022$subsample=="subsample10"), ]

HSEFAfro_subsample <- read_csv("Documents/HSEFAfro_subsample.csv")
HSCFAfro_subsample <- read_csv("Documents/HSCFAfro_subsample.csv")

#EFA analysis on 10% subsample 
ef5fro <- '
efa("efa")*f1 +
efa("efa")*f2 +
efa("efa")*f3 +
efa("efa")*f4 +
efa("efa")*f5 =~
HSIntimatePartner + HSFriend + HSParent + HSOtherRelative + HSMonashMHP + HSCommunityMHP + HSPhoneorOnlineEmergencyService + HSMonashGP + HSCommunityGP + HSMinisterOrReligiousLeader + HSMedicalProfThroughTelehealth + HSDigitalApps + HSHealthOnlineForums
'
efa_f5fro <-
  cfa(model = ef5fro,
      data = HSEFAfro_subsample,
      rotation = "oblimin",
      estimator = "ML",
      missing = "pairwise")
summary(efa_f5fro, fit.measures = T, standardized = T)
fitmea_efa_f5ALL <- fitmeasures(efa_f5fro,
                                   fit.measures = c("cfi","ecvi","mfi", "chisq", "pvalue", "rmsea", "srmr"))
View(fitmea_efa_f5ALL)

#write_excel_csv(HSEFAfro_subsample, file = "HSEFAfro_subsample.csv")
#write_excel_csv(HSCFAfro_subsample, file = "HSCFAfro_subsample.csv")

#EFA analysis on 10% subsample excluding intimate partner
ef5froNIP <- '
efa("efa")*f1 +
efa("efa")*f2 +
efa("efa")*f3 +
efa("efa")*f4 +
efa("efa")*f5 =~
HSFriend + HSParent + HSOtherRelative + HSMonashMHP + HSCommunityMHP + HSPhoneorOnlineEmergencyService + HSMonashGP + HSCommunityGP + HSMedicalProfThroughTelehealth + HSDigitalApps + HSHealthOnlineForums + HSMinisterOrReligiousLeader 
'
efa_f5froNIP <-
  cfa(model = ef5froNIP,
      data = HSEFAfro_subsample,
      rotation = "oblimin",
      estimator = "ML",
      missing = "pairwise")
summary(efa_f5froNIP, fit.measures = T, standardized = T)

fitmea_efa_f5NIP <- fitmeasures(efa_f5froNIP,
                                fit.measures = c("cfi","ecvi","mfi", "chisq", "pvalue", "rmsea", "srmr"))
View(fitmea_efa_f5NIP)

#EFA analysis on 10% subsample excluding intimate partner, four factors only
ef4froNIP <- '
efa("efa")*f1 +
efa("efa")*f2 +
efa("efa")*f3 +
efa("efa")*f4 =~
HSFriend + HSParent + HSOtherRelative + HSMonashMHP + HSCommunityMHP + HSPhoneorOnlineEmergencyService + HSMonashGP + HSCommunityGP + HSMedicalProfThroughTelehealth + HSDigitalApps + HSHealthOnlineForums + HSMinisterOrReligiousLeader
'
efa_f4froNIP <-
  cfa(model = ef4froNIP,
      data = HSEFAfro_subsample,
      rotation = "oblimin",
      estimator = "ML",
      missing = "pairwise")
summary(efa_f4froNIP, fit.measures = T, standardized = T)

#fitmea_efa_f4froNIP <- fitmeasures(efa_f4froNIP, fit.measures = c("cfi","ecvi","mfi", "chisq", "pvalue", "rmsea", "srmr"))
#summary(efa_f4froNIP, fit.measures = c("cfi","ecvi","mfi", "chisq", "pvalue", "rmsea", "srmr"), standardized = T)
fitmea_efa_f4froNIP <- fitmeasures(efa_f4froNIP,
                                fit.measures = c("cfi","ecvi","mfi", "chisq", "pvalue", "rmsea", "srmr"))
#summary(efa_f4froNIP, fitmea_efa_f4froNIP)
View(fitmea_efa_f4froNIP)
#Running CFA from successful EFA 3 analysis excl. IP
FourFactor_modelfro <- "
Factor1 =~ HSFriend + HSParent + HSOtherRelative
Factor2 =~ HSCommunityMHP + HSCommunityGP
Factor3 =~ HSMonashGP + HSMonashMHP
Factor4 =~ HSMedicalProfThroughTelehealth + HSDigitalApps + HSHealthOnlineForums + HSMinisterOrReligiousLeader + HSPhoneorOnlineEmergencyService
"
fit_4f_model <- cfa(FourFactor_modelfro, data = HSCFAfro_subsample, estimator = "ML", 
                    missing = "available.cases", verbose = T)
summary(fit_4f_model, fit.measures = T, standardized = T)
Subsetdata_factor<-predict(fit_4f_model)

fitmea_cfa_f4froNIP <- fitmeasures(fit_4f_model,
                                   fit.measures = c("cfi","ecvi","mfi", "chisq", "pvalue", "rmsea", "srmr"))
View(fitmea_cfa_f4froNIP)
#attempting to fit the latent variables for the 5 factors for subset
HSCFAfro_subsampleBIND<-cbind(HSCFAfro_subsample, Subsetdata_factor)
View(HSCFAfro_subsampleBIND)

#Printing final data
write_excel_csv(HSEFAfro_subsample, file = "HSEFAfro_subsample.csv")
write_excel_csv(HSCFAfro_subsample, file = "HSCFAfro_subsample.csv")
write_excel_csv(HSCFAfro_subsampleBIND, file = "HSCFAfro_subsampleBIND.csv")

library(dplyr)
HelpSeeking <- select(DBHSfro2022, HSIntimatePartner, HSFriend, HSParent, HSOtherRelative, HSMonashMHP, HSCommunityMHP, HSPhoneorOnlineEmergencyService, HSMonashGP, HSCommunityGP, HSMinisterOrReligiousLeader, HSMedicalProfThroughTelehealth, HSDigitalApps, HSHealthOnlineForums)
PerceivedStress <- select(DBHSfro2022, PSS1, PSS2, PSS3, PSS4R, PSS5R, PSS6, PSS7R, PSS8, PSS9, PSS10)
Anxiety <- select(DBHSfro2022, ANX1, ANX2, ANX3, ANX4)
install.packages("psych")
library(psych)
alpha(HelpSeeking)
alpha(PerceivedStress)
alpha(Anxiety)
install.packages("summarytools")
table(DBHSfro2022$Gender)

View(HSCFAfro_subsampleBIND)

#Running ANCOVA analysis with age, gender, YOD, and international vs domestic status as covariates, help seeking factors as outcomes and Wellbeing grouping as grouping variable



#MRA analysis for stress
stressMRA <- lm(formula = PSStotal ~ Factor1 + Factor2 + Factor3 + Factor4 + HSIntimatePartner, data=HSCFAfro_subsampleBIND)
summary(stressMRA)
coefficients(stressMRA) # model coefficients
confint(stressMRA, level=0.95) # CIs for model parameters
fitted(stressMRA) # predicted values
residuals(stressMRA) # residuals
anova(stressMRA) # anova table
vcov(stressMRA) # covariance matrix for model parameters
#influence(stressMRA) # regression diagnostics

AnxietyMRA <- lm(formula = ANXtotal ~ Factor1 + Factor2 + Factor3 + Factor4 + HSIntimatePartner, data=HSCFAfro_subsampleBIND)
summary(AnxietyMRA)
coefficients(AnxietyMRA) # model coefficients
confint(AnxietyMRA, level=0.95) # CIs for model parameters
fitted(AnxietyMRA) # predicted values
residuals(AnxietyMRA) # residuals
anova(AnxietyMRA) # anova table
vcov(AnxietyMRA) # covariance matrix for model parameters

#assumption checking
#Normality
HSIntimatePartner<-(HSCFAfro_subsampleBIND$HSIntimatePartner)
Factor1<-(HSCFAfro_subsampleBIND$Factor1)
Factor2<-(HSCFAfro_subsampleBIND$Factor2)
Factor3<-(HSCFAfro_subsampleBIND$Factor3)
Factor4<-(HSCFAfro_subsampleBIND$Factor4)
Stress<-(HSCFAfro_subsampleBIND$PSStotal)
Anxiety<-(HSCFAfro_subsampleBIND$ANXtotal)
hist(HSIntimatePartner, col='steelblue', main='Intimate Partner')#Violated twin peaks
hist(Factor1, col='steelblue', main='Factor 1')#Assumed
hist(Factor2, col='steelblue', main='Factor 2')#Assumed
hist(Factor3, col='steelblue', main='Factor 3')#positively skewed
hist(Factor4, col='steelblue', main='Factor 4')#positively skewed
hist(Stress, col='steelblue', main='Total stress')#Assumed
hist(Anxiety, col='steelblue', main='Anxiety')#Assumed

qqnorm(HSIntimatePartner)
qqline(HSIntimatePartner) #Violated and weird 
qqnorm(Factor3)
qqline(Factor3)#seems okay, some deviation at ends
qqnorm(Factor4)
qqline(Factor4)#seems okay, some deviation at ends

Testingsample <- data.frame(HSIP = c(HSIntimatePartner), F1 = c(Factor1), F2 = c(Factor2), F3 = c(Factor3), F4 = c(Factor4), Stress = c(Stress), Anxiety = c(Anxiety))
#ead(Testingsample)
#Testingsample$mahalanobis<-mahalanobis(Testingsample, colMeans(Testingsample), cov = Testingsample)
#Mahalanobis distance isnt working, cov requires a square 
write_excel_csv(HSCFAfro_subsampleBIND, file = "HSCFAfro_subsampleBIND.csv")
#Visualising correlations with stress
library(ggplot2)
r <- round(cor(HSCFAfro_subsampleBIND$HSIntimatePartner, HSCFAfro_subsampleBIND$PSStotal), 2)
p <- cor.test(HSCFAfro_subsampleBIND$HSIntimatePartner, HSCFAfro_subsampleBIND$PSStotal)$p.value
ggplot(HSCFAfro_subsampleBIND, aes(y=HSIntimatePartner, x=PSStotal)) + 
  geom_point() + 
  geom_smooth(method="lm", col="black") + 
  annotate("text", x=5.5, y=5.5, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=5, y=5, label=paste0("p = ", round(p, 3)), hjust=0) +
  theme_classic() 

r <- round(cor(HSCFAfro_subsampleBIND$Factor1, HSCFAfro_subsampleBIND$PSStotal), 2)
p <- cor.test(HSCFAfro_subsampleBIND$Factor1, HSCFAfro_subsampleBIND$PSStotal)$p.value
ggplot(HSCFAfro_subsampleBIND, aes(y=Factor1, x=PSStotal)) + 
  geom_point() + 
  geom_smooth(method="lm", col="black") + 
  annotate("text", x=40, y=7, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=40, y=7, label=paste0("p = ", round(p, 3)), hjust=0) +
  theme_classic() 

r <- round(cor(HSCFAfro_subsampleBIND$Factor2, HSCFAfro_subsampleBIND$PSStotal), 2)
p <- cor.test(HSCFAfro_subsampleBIND$Factor2, HSCFAfro_subsampleBIND$PSStotal)$p.value
ggplot(HSCFAfro_subsampleBIND, aes(y=Factor2, x=PSStotal)) + 
  geom_point() + 
  geom_smooth(method="lm", col="black") + 
  annotate("text", x=40, y=7, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=40, y=7, label=paste0("p = ", round(p, 3)), hjust=0) +
  theme_classic()

r <- round(cor(HSCFAfro_subsampleBIND$Factor3, HSCFAfro_subsampleBIND$PSStotal), 2)
p <- cor.test(HSCFAfro_subsampleBIND$Factor3, HSCFAfro_subsampleBIND$PSStotal)$p.value
ggplot(HSCFAfro_subsampleBIND, aes(y=Factor3, x=PSStotal)) + 
  geom_point() + 
  geom_smooth(method="lm", col="black") + 
  annotate("text", x=40, y=7, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=40, y=7, label=paste0("p = ", round(p, 3)), hjust=0) +
  theme_classic()

r <- round(cor(HSCFAfro_subsampleBIND$Factor4, HSCFAfro_subsampleBIND$PSStotal), 2)
p <- cor.test(HSCFAfro_subsampleBIND$Factor4, HSCFAfro_subsampleBIND$PSStotal)$p.value
ggplot(HSCFAfro_subsampleBIND, aes(y=Factor4, x=PSStotal)) + 
  geom_point() + 
  geom_smooth(method="lm", col="black") + 
  annotate("text", x=40, y=7, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=40, y=7, label=paste0("p = ", round(p, 3)), hjust=0) +
  theme_classic()
#Visualising correlations with anxiety
r <- round(cor(HSCFAfro_subsampleBIND$HSIntimatePartner, HSCFAfro_subsampleBIND$ANXtotal), 2)
p <- cor.test(HSCFAfro_subsampleBIND$HSIntimatePartner, HSCFAfro_subsampleBIND$ANXtotal)$p.value
ggplot(HSCFAfro_subsampleBIND, aes(y=HSIntimatePartner, x=ANXtotal)) + 
  geom_point() + 
  geom_smooth(method="lm", col="black") + 
  annotate("text", x=40, y=7, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=40, y=7, label=paste0("p = ", round(p, 3)), hjust=0) +
  theme_classic() 

r <- round(cor(HSCFAfro_subsampleBIND$Factor1, HSCFAfro_subsampleBIND$ANXtotal), 2)
p <- cor.test(HSCFAfro_subsampleBIND$Factor1, HSCFAfro_subsampleBIND$ANXtotal)$p.value
ggplot(HSCFAfro_subsampleBIND, aes(y=Factor1, x=ANXtotal)) + 
  geom_point() + 
  geom_smooth(method="lm", col="black") + 
  annotate("text", x=40.5, y=7.5, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=40, y=7, label=paste0("p = ", round(p, 3)), hjust=0) +
  theme_classic() 

r <- round(cor(HSCFAfro_subsampleBIND$Factor2, HSCFAfro_subsampleBIND$ANXtotal), 2)
p <- cor.test(HSCFAfro_subsampleBIND$Factor2, HSCFAfro_subsampleBIND$ANXtotal)$p.value
ggplot(HSCFAfro_subsampleBIND, aes(y=Factor2, x=ANXtotal)) + 
  geom_point() + 
  geom_smooth(method="lm", col="black") + 
  annotate("text", x=40.5, y=7.5, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=40, y=7, label=paste0("p = ", round(p, 3)), hjust=0) +
  theme_classic()

r <- round(cor(HSCFAfro_subsampleBIND$Factor3, HSCFAfro_subsampleBIND$ANXtotal), 2)
p <- cor.test(HSCFAfro_subsampleBIND$Factor3, HSCFAfro_subsampleBIND$ANXtotal)$p.value
ggplot(HSCFAfro_subsampleBIND, aes(y=Factor3, x=ANXtotal)) + 
  geom_point() + 
  geom_smooth(method="lm", col="black") + 
  annotate("text", x=40.5, y=7.5, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=40, y=7, label=paste0("p = ", round(p, 3)), hjust=0) +
  theme_classic()

r <- round(cor(HSCFAfro_subsampleBIND$Factor4, HSCFAfro_subsampleBIND$ANXtotal), 2)
p <- cor.test(HSCFAfro_subsampleBIND$Factor4, HSCFAfro_subsampleBIND$ANXtotal)$p.value
ggplot(HSCFAfro_subsampleBIND, aes(y=Factor4, x=ANXtotal)) + 
  geom_point() + 
  geom_smooth(method="lm", col="black") + 
  annotate("text", x=40.5, y=7.5, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=40, y=7, label=paste0("p = ", round(p, 3)), hjust=0) +
  theme_classic()

r <- round(cor(HSCFAfro_subsampleBIND$HSIntimatePartner, HSCFAfro_subsampleBIND$ANXtotal), 2)
p <- cor.test(HSCFAfro_subsampleBIND$HSIntimatePartner, HSCFAfro_subsampleBIND$ANXtotal)$p.value
ggplot(HSCFAfro_subsampleBIND, aes(y=HSIntimatePartner, x=ANXtotal)) + 
  geom_point() + 
  geom_smooth(method="lm", col="black") + 
  annotate("text", x=40.5, y=7.5, label=paste0("r = ", r), hjust=0) +
  annotate("text", x=40, y=7, label=paste0("p = ", round(p, 3)), hjust=0) +
  theme_classic() 

#New analysis focusing on ANCOVA 
write_excel_csv(HSCFAfro_subsampleBIND2, file = "HSCFAfro_subsampleBIND2.csv")
CFAfro_ANCOVA <- read_csv("Documents/HSCFAfro_subsampleBIND2.csv")
case_when(CFAfro_ANCOVA$WHOgrouping %in% 0:28 ~ 1, CFAfro_ANCOVA$WHOgrouping %in% 29:49 ~ 2, CFAfro_ANCOVA$WHOgrouping %in% 50:100 ~ 3)
