library(ggplot2)
library(patternplot)
library(lme4)
library(visreg)
library(tidyr)
library(magrittr) 
library(dplyr)
library(ggpubr)
library(emmeans)
library(lsmeans)
library(MuMIn)
library(multcomp)
library(zoo)
library(ggpmisc)
library(rstatix)
library(lmerTest)


setwd('/Users/alex/AlexBadea_MyCodes/wind-main/')
masterFile<-'APOE_MWM_AB281221.csv'
outpath='./Figures/learning/'

info<-read.csv(masterFile, header=TRUE)
df<-data.frame(info)
df<-subset(df, !is.nan(df$NormSWDist))

df2 <-subset(df, (APOE=='E22') & (Age==12) & (Diet=='Control')) #Keeps only Genotypes 2/2, 3/3, and 4/4
df2m<-subset(df, (APOE=='E22') & (Age==12) & (Diet=='Control') & (Sex=='M')) #Keeps only Genotypes 2/2, 3/3, and 4/4
df2f<-subset(df, (APOE=='E22') & (Age==12) & (Diet=='Control') & (Sex=='F')) #Keeps only Genotypes 2/2, 3/3, and 4/4
length(unique(df2$Animal))
length(unique(df2m$Animal))
length(unique(df2f$Animal))
mean(df2$Age.months, na.rm = TRUE)
sd(as.numeric(df2$Age.months), na.rm = TRUE)

df3<-subset(df, (APOE=='E33') & (Age==12) & (Diet=='Control'))
df3m<-subset(df, (APOE=='E33') & (Age==12) & (Diet=='Control') & (Sex=='M')) #Keeps only Genotypes 2/2, 3/3, and 4/4
df3f<-subset(df, (APOE=='E33') & (Age==12) & (Diet=='Control') & (Sex=='F')) #Keeps only Genotypes 2/2, 3/3, and 4/4
length(unique(df3$Animal))
length(unique(df3m$Animal))
length(unique(df3f$Animal))
mean(df3$Age.months, na.rm = TRUE)
sd(as.numeric(df3$Age.months), na.rm = TRUE)

df4<-subset(df, (APOE=='E44') & (Age==12) & (Diet=='Control'))
df4m<-subset(df, (APOE=='E44') & (Age==12)& (Diet=='Control') & (Sex=='M')) #Keeps only Genotypes 2/2, 3/3, and 4/4
df4f<-subset(df, (APOE=='E44') & (Age==12) & (Diet=='Control') & (Sex=='F')) #Keeps only Genotypes 2/2, 3/3, and 4/4
length(unique(df4$Animal))
length(unique(df4m$Animal))
length(unique(df4f$Animal))
mean(df4$Age.months, na.rm = TRUE)
sd(as.numeric(df4$Age.months), na.rm = TRUE)

AnimalTable<-matrix(nrow=5, ncol=6)
AnimalTable[1,]=c('', '', '', '', '', '')
AnimalTable[2,]=c('Genotype', 'Total', 'Males', 'Females', 'Mean_Age', 'SD_Age')
AnimalTable[3,]=c('APOE2', length(unique(df2$Animal)), length(unique(df2m$Animal)), length(unique(df2f$Animal)), mean(df2$Age.months, na.rm = TRUE), sd(as.numeric(df2$Age.months), na.rm = TRUE))
AnimalTable[4,]=c('APOE3', length(unique(df3$Animal)), length(unique(df3m$Animal)), length(unique(df3f$Animal)), mean(df3$Age.months, na.rm = TRUE), sd(as.numeric(df3$Age.months), na.rm = TRUE))
AnimalTable[5,]=c('APOE4', length(unique(df4$Animal)), length(unique(df4m$Animal)), length(unique(df4f$Animal)), mean(df4$Age.months, na.rm = TRUE), sd(as.numeric(df4$Age.months), na.rm = TRUE))

myanimalfile<-paste(outpath,'AnimalGroups_APOE234.csv', sep='')
#write.table(mytTable, file=myfile, col.names = mycolnames , sep = "," , row.names = F,append=TRUE)
write.table(AnimalTable, file=myanimalfile, sep=",", row.names=F, append=TRUE, col.names=F)

#-----------
dfAll<-rbind(df2, df3, df4)

#distance calibration in AnyMaze
dfAll$Distance=dfAll$Distance/10
dfAll$NE.distance=dfAll$NE.distance/10
dfAll$NW.distance=dfAll$NW.distance/10
dfAll$SE.distance=dfAll$SE.distance/10
dfAll$SW.distance=dfAll$SW.distance/10


#learning trials only in dfFin
df1<-subset(dfAll, (Stage=='Day1'))
df2<-subset(dfAll, (Stage=='Day2'))
df3<-subset(dfAll, (Stage=='Day3'))
df4<-subset(dfAll, (Stage=='Day4'))
df5<-subset(dfAll, (Stage=='Day5'))
dfFin<-rbind(df1,df2,df3,df4,df5) #dfFin contains info from all regular trials, no probes

#Normalize time and distance in target region
#dfFin$NormSWTime<-dfFin$SW.time/dfFin$Duration
#dfFin$NormSWDist<-dfFin$SW.distance/dfFin$Distance
#dfFin<-subset(dfFin, (NormSWTime <= 1)) #is this necessary???


dfFin$Animal <- as.factor(dfFin$Animal)
dfFin$APOE <- as.factor(dfFin$APOE)
dfFin$Sex <- as.factor(dfFin$Sex)
dfFin$Stage <- as.factor(dfFin$Stage)
dfFin$Treatment <- as.factor(dfFin$Treatment)
dfFin$Code <- as.factor(dfFin$Code)
dfFin$Trial <- as.factor(dfFin$Trial)
dfFin$NormSWDist<-as.numeric(dfFin$NormSWDist)
dfFin$NormSWDTime<-as.numeric(dfFin$NormSWTime)


#Averages 4 trials per day for each mouse
dfAveraged<-aggregate(.~Animal+APOE+Sex+Test+Treatment+Code+Stage, dfFin, mean, na.action=na.pass)
dfAveraged<-aggregate(x = dfFin, by = list(dfFin$Animal,dfFin$APOE,dfFin$Code,dfFin$Sex,dfFin$Stage,dfFin$Treatment), FUN = "mean", na.action=na.pass)

dfAveraged$APOE <- NULL
dfAveraged$Sex <- NULL
dfAveraged$Stage <- NULL
dfAveraged$Treatment <- NULL
dfAveraged$Code <- NULL
dfAveraged$Trial <- NULL
dfAveraged$Animal <- NULL

colnames(dfAveraged[1:6])<-colnames(dfFin[1:6])

dfAveraged$Animal <- dfAveraged$Group.1
dfAveraged$APOE <- dfAveraged$Group.2
dfAveraged$Code <- dfAveraged$Group.3
dfAveraged$Sex <- dfAveraged$Group.4
dfAveraged$Stage <- dfAveraged$Group.5
dfAveraged$Treatment <- dfAveraged$Group.6

sink(paste(outpath, "LearningDist_Tukeymodels.txt"))

#https://stats.stackexchange.com/questions/405661/repeated-measures-regression-in-r
#https://stats.stackexchange.com/questions/136495/addressing-note-results-may-be-misleading-due-to-involvement-in-interactions

lm_Distance_trials <- lmer(Distance ~ APOE*Sex*Stage+(1|Animal), dfAveraged, REML = TRUE)

library(effsize)
dfAveraged_tmps <-
  dfAveraged %>% 
  mutate(cut = factor(APOE, levels = c('E33', 'E44')))
levels(dfAveraged_tmps$cut)
cohen.d(Distance ~ cut, dfAveraged_tmps)

anova(lm_Distance_trials)
eta_squared(lm_Distance_trials, alternative='two.sided')
effectsize::cohens_f(lm_Distance_trials, alternative='two.sided')
r.squaredGLMM(lm_Distance_trials)
qqnorm(resid(lm_Distance_trials))

emmeans(lm_Distance_trials , ~ APOE*Sex, contr="tukey")
summary(glht(lm_Distance_trials, emm(pairwise ~ APOE)))
summary(glht(lm_Distance_trials, emm(pairwise ~ APOE*Sex, adjust="sidak"))) #consider save to file
summary(glht(lm_Distance_trials, emm(pairwise ~ APOE| Sex)))
summary(glht(lm_Distance_trials, emm(pairwise ~ Sex | APOE)))
sink()

sink(paste(outpath, "LearningNormSWDist_Tukeymodels.txt"))
lm_NormSWDist_trials <- lmer(NormSWDist ~ APOE*Sex*Stage+(1|Animal), dfAveraged, REML=TRUE)


anova(lm_NormSWDist_trials)
eta_squared(lm_NormSWDist_trials,alternative='two.sided')
effectsize::cohens_f(lm_NormSWDist_trials,alternative='two.sided')
r.squaredGLMM(lm_NormSWDist_trials)
qqnorm(resid(lm_NormSWDist_trials))

#emmeans(lm_Distance_trials , ~ APOE*Sex, contr="tukey")
summary(glht(lm_NormSWDist_trials, emm(pairwise ~ APOE)))
summary(glht(lm_NormSWDist_trials, emm(pairwise ~ APOE*Sex, adjust="sidak"))) #consider save to file
summary(glht(lm_NormSWDist_trials, emm(pairwise ~ APOE| Sex)))
summary(glht(lm_NormSWDist_trials, emm(pairwise ~ Sex | APOE)))
sink()

sink(paste(outpath, "Winding_LearningDay1Trial4_Tukeymodels.txt"))
dfFinS1T4<-subset(dfFin, dfFin$Stage=='Day1' & dfFin$Trial==4)
lm_Day1Trial4_trials <- lm(Winding ~ APOE*Sex, dfFinS1T4)
anova(lm_Day1Trial4_trials)
eta_squared(lm_Day1Trial4_trials,alternative='two.sided')
effectsize::cohens_f(lm_Day1Trial4_trials, alternative='two.sided')
r.squaredGLMM(lm_Day1Trial4_trials)
qqnorm(resid(lm_Day1Trial4_trials))


#emmeans(lm_Day1Trial4_trials , ~ APOE*Sex, contr="tukey")
summary(glht(lm_Day1Trial4_trials, emm(pairwise ~ APOE)))
summary(glht(lm_Day1Trial4_trials, emm(pairwise ~ APOE*Sex, adjust="sidak"))) #consider save to file
summary(glht(lm_Day1Trial4_trials, emm(pairwise ~ APOE| Sex)))
summary(glht(lm_Day1Trial4_trials, emm(pairwise ~ Sex | APOE)))
sink()


sink(paste(outpath, "Winding_Learning_Tukeymodels.txt"))
lm_winding_trials <- lmer(Winding ~ APOE*Sex*Stage+(1|Animal), dfAveraged, REML = TRUE)
anova(lm_winding_trials)

eta_squared(lm_winding_trials,alternative='two.sided')
effectsize::cohens_f(lm_winding_trials,alternative='two.sided')
r.squaredGLMM(lm_winding_trials)
qqnorm(resid(lm_winding_trials))


#emmeans(lm_Day1Trial4_trials , ~ APOE*Sex, contr="tukey")
summary(glht(lm_winding_trials, emm(pairwise ~ APOE)))
summary(glht(lm_winding_trials, emm(pairwise ~ APOE*Sex, adjust="sidak"))) #consider save to file
summary(glht(lm_winding_trials, emm(pairwise ~ APOE| Sex)))
summary(glht(lm_winding_trials, emm(pairwise ~ Sex | APOE)))
sink()

#Extract probe trials only
dfp1<-subset(dfAll, (Stage=='Probe_D5'))
dfp2<-subset(dfAll, (Stage=='Probe_D8'))



sink(paste(outpath, "Dist_D5_Tukeymodels.txt"))
lm_Distance_d5 <- lm(Distance ~ APOE*Sex, dfp1)
anova(lm_Distance_d5)
eta_squared(lm_Distance_d5,alternative='two.sided')
effectsize::cohens_f(lm_Distance_d5,alternative='two.sided')
r.squaredGLMM(lm_Distance_d5)
qqnorm(resid(lm_Distance_d5))
#emmeans(lm_Distance_d5 , ~ APOE*Sex, contr="tukey")
summary(glht(lm_Distance_d5, emm(pairwise ~ APOE)))
summary(glht(lm_Distance_d5, emm(pairwise ~ APOE*Sex, adjust="sidak"))) #consider save to file
summary(glht(lm_Distance_d5, emm(pairwise ~ APOE| Sex)))
summary(glht(lm_Distance_d5, emm(pairwise ~ Sex | APOE)))
sink()


sink(paste(outpath, "NormSWDist_D5_Tukeymodels.txt"))
lm_NormSWDist_d5 <- lm(NormSWDist ~ APOE*Sex, dfp1)
anova(lm_NormSWDist_d5)
eta_squared(lm_NormSWDist_d5,alternative='two.sided')
effectsize::cohens_f(lm_NormSWDist_d5,alternative='two.sided')
r.squaredGLMM(lm_NormSWDist_d5)
qqnorm(resid(lm_NormSWDist_d5))
emmeans(lm_NormSWDist_d5 , ~ APOE*Sex, contr="tukey")
summary(glht(lm_NormSWDist_d5, emm(pairwise ~ APOE)))
summary(glht(lm_NormSWDist_d5, emm(pairwise ~ APOE*Sex, adjust="sidak"))) #consider save to file
summary(glht(lm_NormSWDist_d5, emm(pairwise ~ APOE| Sex)))
summary(glht(lm_NormSWDist_d5, emm(pairwise ~ Sex | APOE)))
sink()


sink(paste(outpath, "Winding_D5_Tukeymodels.txt"))
lm_Winding_d5 <- lm(Winding ~ APOE*Sex, dfp1)
anova(lm_Winding_d5)
eta_squared(lm_Winding_d5,alternative='two.sided')
effectsize::cohens_f(lm_Winding_d5,alternative='two.sided')
r.squaredGLMM(lm_Winding_d5)
qqnorm(resid(lm_Winding_d5))
summary(glht(lm_Winding_d5, emm(pairwise ~ APOE)))
emmeans(lm_Winding_d5 , ~ APOE*Sex, contr="tukey")
#33333
summary(glht(lm_Winding_d5, emm(pairwise ~ APOE*Sex, adjust="sidak"))) #consider save to file
summary(glht(lm_Winding_d5, emm(pairwise ~ APOE| Sex)))
summary(glht(lm_Winding_d5, emm(pairwise ~ Sex | APOE)))
#3333
sink()


sink(paste(outpath, "Dist_D8_Tukeymodels.txt"))
lm_Distance_d8 <- lm(Distance ~ APOE*Sex, dfp2)
anova(lm_Distance_d8)
eta_squared(lm_Distance_d8,alternative='two.sided')
effectsize::cohens_f(lm_Distance_d8,alternative='two.sided')
r.squaredGLMM(lm_Distance_d8)
qqnorm(resid(lm_Distance_d8))
summary(glht(lm_Distance_d8, emm(pairwise ~ APOE*Sex, adjust="sidak"))) #consider save to file
summary(glht(lm_Distance_d8, emm(pairwise ~ APOE| Sex)))
summary(glht(lm_Distance_d8, emm(pairwise ~ Sex | APOE)))
sink()

#NormSWDist
sink(paste(outpath, "NormSWDist_D8_Tukeymodels.txt"))
lm_NormSWDist_d8 <- lm(NormSWDist ~ APOE*Sex, dfp2)
anova(lm_NormSWDist_d8)
eta_squared(lm_NormSWDist_d8,alternative='two.sided')
effectsize::cohens_f(lm_NormSWDist_d8,alternative='two.sided')
r.squaredGLMM(lm_NormSWDist_d8)
qqnorm(resid(lm_NormSWDist_d8))
summary(glht(lm_NormSWDist_d8, emm(pairwise ~ APOE*Sex, adjust="sidak"))) #consider save to file
summary(glht(lm_NormSWDist_d8, emm(pairwise ~ APOE| Sex)))
summary(glht(lm_NormSWDist_d8, emm(pairwise ~ Sex | APOE)))
sink()


sink(paste(outpath, "Winding_D8_Tukeymodels.txt"))
'genotype only, no sex'
anova(lm(Winding ~ APOE, dfp2))

lm_Winding_d8 <- lm(Winding ~ APOE*Sex, dfp2)
anova(lm_Winding_d8)
eta_squared(lm_Winding_d8,alternative='two.sided')
effectsize::cohens_f(lm_Winding_d8,alternative='two.sided')
r.squaredGLMM(lm_Winding_d8)
qqnorm(resid(lm_Winding_d8))
summary(glht(lm_NormSWDist_d8, emm(pairwise ~ APOE)))
summary(glht(lm_Winding_d8, emm(pairwise ~ APOE*Sex, adjust="sidak"))) #consider save to file
summary(glht(lm_Winding_d8, emm(pairwise ~ APOE| Sex)))
summary(glht(lm_Winding_d8, emm(pairwise ~ Sex | APOE)))

sink()
