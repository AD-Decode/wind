library(ggplot2)
library(patternplot)
library(lme4)
library(visreg)
library(tidyr)
library(magrittr) 
library(dplyr)
library(reshape2)
library(ggpubr)

dirin_trial<-'/Users/alex/AlexBadea_MyPapers/DavidDunson/dropbox/WindingNumber/TrackingData/Trials/'


dirin<-'/Users/alex/AlexBadea_MyPapers/DavidDunson/dropbox/WindingNumber/TrackingData/Probes/'
dirout<-'/Users/alex/AlexBadea_MyPapers/DavidDunson/Figures/'
  
animal<- '190715_8' # APOE2 red
animal<- '190715_6' # APOE2 red
#animal<- '190715_9' # APOE2 red
Genotype2<-'APOE22'
Genotype<-Genotype2

filein214<-paste(dirin_trial,animal,'_Day1_T4_positions.csv',sep='') #genotype 2, day1, trial 4
filein25<-paste(dirin,animal,'_Probe_D5_T1_positions.csv',sep='') # genotype 2, day 5
filein28<-paste(dirin,animal,'_Probe_D8_T1_positions.csv',sep='') # genotype 2, day 8

mytraj_data214 <- read.csv(filein214)
mytraj_data214 <- na.omit(mytraj_data214)
mytraj_data214['Genotype']=rep(Genotype2, length(mytraj_data214$Centre.position.X))
mytraj_data214['Animal']=rep(animal, length(mytraj_data214$Centre.position.X))

mytraj_data25 <- read.csv(filein25)
mytraj_data25 <- na.omit(mytraj_data25)
mytraj_data25['Genotype']=rep(Genotype2, length(mytraj_data25$Centre.position.X))
mytraj_data25['Animal']=rep(animal, length(mytraj_data25$Centre.position.X))

mytraj_data28 <- read.csv(filein28)
mytraj_data28 <- na.omit(mytraj_data28)
mytraj_data28['Genotype']=rep(Genotype2, length(mytraj_data28$Centre.position.X))
mytraj_data28['Animal']=rep(animal, length(mytraj_data28$Centre.position.X))


animal<-'190909_9' #APOE3 green
#animal<-'190909_1' #APOE3 green
#animal<-'190715_1' #APOE3 green

Genotype3<-'APOE33'

filein314<-paste(dirin_trial,animal,'_Day1_T4_positions.csv',sep='') #genotype 3, day1, trial 4
filein35<-paste(dirin,animal,'_Probe_D5_T1_positions.csv',sep='') # genotype 3, day 5
filein38<-paste(dirin,animal,'_Probe_D8_T1_positions.csv',sep='') # genotype 3, day 8

mytraj_data314 <- read.csv(filein314)
mytraj_data314 <- na.omit(mytraj_data314)
mytraj_data314['Genotype']=rep(Genotype3, length(mytraj_data314$Centre.position.X))
mytraj_data314['Animal']=rep(animal, length(mytraj_data314$Centre.position.X))

mytraj_data35 <- read.csv(filein35)
mytraj_data35 <- na.omit(mytraj_data35)
mytraj_data35['Genotype']=rep(Genotype3, length(mytraj_data35$Centre.position.X))
mytraj_data35['Animal']=rep(animal, length(mytraj_data35$Centre.position.X))

mytraj_data38 <- read.csv(filein38)
mytraj_data38 <- na.omit(mytraj_data38)
mytraj_data38['Genotype']=rep(Genotype3, length(mytraj_data38$Centre.position.X))
mytraj_data38['Animal']=rep(animal, length(mytraj_data38$Centre.position.X))


animal<-'191006_6' #APOE4  blue
#animal<-'191006_9' #APOE4  blue
Genotype4<-'APOE44'

filein414<-paste(dirin_trial,animal,'_Day1_T4_positions.csv',sep='') #genotype 4, day1, trial 4
filein45<-paste(dirin,animal,'_Probe_D5_T1_positions.csv',sep='') # genotype 3, day 5
filein48<-paste(dirin,animal,'_Probe_D8_T1_positions.csv',sep='') # genotype 3, day 8

mytraj_data414 <- read.csv(filein414)
mytraj_data414 <- na.omit(mytraj_data414)
mytraj_data414['Genotype']=rep(Genotype4, length(mytraj_data414$Centre.position.X))
mytraj_data414['Animal']=rep(animal, length(mytraj_data414$Centre.position.X))

mytraj_data45 <- read.csv(filein45)
mytraj_data45 <- na.omit(mytraj_data45)
mytraj_data45['Genotype']=rep(Genotype4, length(mytraj_data45$Centre.position.X))
mytraj_data45['Animal']=rep(animal, length(mytraj_data45$Centre.position.X))

mytraj_data48 <- read.csv(filein48)
mytraj_data48 <- na.omit(mytraj_data48)
mytraj_data48['Genotype']=rep(Genotype4, length(mytraj_data48$Centre.position.X))
mytraj_data48['Animal']=rep(animal, length(mytraj_data48$Centre.position.X))

mytraj_data14<-rbind(mytraj_data214, mytraj_data314, mytraj_data414)
mytraj_data5<-rbind(mytraj_data25, mytraj_data35, mytraj_data45)
mytraj_data8<-rbind(mytraj_data28, mytraj_data38, mytraj_data48)

plotD14<-ggplot(data=mytraj_data14, aes(x=Centre.position.X, y=Centre.position.Y, group = Genotype, color=Genotype)) +  
  scale_color_manual(values=c('blueviolet', 'chartreuse1', 'red'))+
  geom_path (linetype=1, size=0.5)+
  theme_classic()+
  facet_grid(. ~ Genotype, scales='fixed', space='fixed') + 
  coord_fixed(ratio=1)+
  #stat_smooth(method = "lm") +
  #background_grid(major = 'xy', minor = "none") + # add thin horizontal lines 
  theme_bw()+
  ggtitle(paste('A:Learning Trial 4,  Day 1'))

plotD5<-ggplot(data=mytraj_data5, aes(x=Centre.position.X, y=Centre.position.Y, group = Genotype, color=Genotype)) +  
  scale_color_manual(values=c('blueviolet', 'chartreuse1', 'red'))+
  geom_path (linetype=1, size=0.5)+
  theme_classic()+
  facet_grid(. ~ Genotype, scales='fixed', space='fixed') + 
  coord_fixed(ratio=1)+
  #stat_smooth(method = "lm") +
  #background_grid(major = 'xy', minor = "none") + # add thin horizontal lines 
  theme_bw()+
  ggtitle(paste('B:Probe Day 5'))
  

plotD8<-ggplot(data=mytraj_data8, aes(x=Centre.position.X, y=Centre.position.Y, group = Genotype, color=Genotype)) +  
  scale_color_manual(values=c('blueviolet', 'chartreuse1', 'red'))+
  geom_path (linetype=1, size=0.5)+
  theme_classic()+
  facet_grid(. ~ Genotype, scales='fixed', space='fixed') + 
  coord_fixed(ratio=1)+
  #stat_smooth(method = "lm") +
  #background_grid(major = 'xy', minor = "none") + # add thin horizontal lines 
  theme_bw()+
  ggtitle(paste('C:Probe Day 8'))
 

ggarrange(plotD14, plotD5, plotD8,
          #labels = c("A: Day 5", "B: Day 8"),
          ncol = 1, nrow = 3)

ggsave(paste(dirout,'my3animals', 'Probe_D1T4D5D8_combo_example3.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=7, height=9, unit=c("in"), dpi=300)


#myx<-mytraj_data$Centre.position.X
#myy<-mytraj_data$Centre.position.Y

# #APOE4
# filein <- '/Users/alex/AlexBadea_MyPapers/DavidDunson/dropbox/WindingNumber/TrackingData/Probe_D8/191006_7_Probe_D8_T1_positions.csv'
# filein <- '/Users/alex/AlexBadea_MyPapers/DavidDunson/dropbox/WindingNumber/TrackingData/Probe_D8/191006_1_Probe_D8_T1_positions.csv'
# filein <- '/Users/alex/AlexBadea_MyPapers/DavidDunson/dropbox/WindingNumber/TrackingData/Probe_D8/191006_1_Probe_D8_T1_positions.csv'
# filein <- '/Users/alex/AlexBadea_MyPapers/DavidDunson/dropbox/WindingNumber/TrackingData/Probe_D8/191006_3_Probe_D8_T1_positions.csv'
# filein <- '/Users/alex/AlexBadea_MyPapers/DavidDunson/dropbox/WindingNumber/TrackingData/Probe_D8/191006_9_Probe_D8_T1_positions.csv'
# filein <- '/Users/alex/AlexBadea_MyPapers/DavidDunson/dropbox/WindingNumber/TrackingData/Probe_D8/191006_6_Probe_D8_T1_positions.csv'
# #APOE3
# 
# filein <- '/Users/alex/AlexBadea_MyPapers/DavidDunson/dropbox/WindingNumber/TrackingData/Probe_D8/190909_12_Probe_D8_T1_positions.csv'
# filein <- '/Users/alex/AlexBadea_MyPapers/DavidDunson/dropbox/WindingNumber/TrackingData/Probe_D8/190909_10_Probe_D8_T1_positions.csv'
# #filein <- '/Users/alex/AlexBadea_MyPapers/DavidDunson/dropbox/WindingNumber/TrackingData/Probe_D8/190909_11_Probe_D8_T1_positions.csv'



# ggplot(data=mytraj_data, aes(x=Centre.position.X, y=Centre.position.Y, group = Genotype)) +  
#   #geom_point()+
#   geom_smooth(method=lm, se=FALSE, fullrange=TRUE, aes(fill=Genotype))+
#   geom_path (linetype=1, size=0.5, color=Genotype)+
#   theme_classic()+
#   facet_grid(. ~ Genotype) + 
#   stat_smooth(method = "lm") +
#   #background_grid(major = 'xy', minor = "none") + # add thin horizontal lines 
#   theme_bw() +
#   theme(axis.text.x = element_text(face="bold",  size=14, angle=0),
#         axis.text.y = element_text(face="bold", size=14, angle=0),
#         axis.line.x = element_line(colour = 'black', size=0.5),
#         axis.line.y = element_line(colour = 'black', size=0.5),
#         # panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         # panel.border = element_blank(),
#         panel.background = element_blank())

#ggdraw() + 
#  draw_plot(plot1, 0, .5, 1, .5)

#ggsave(paste(dirout,'Example1', 'Probe_D8.pdf',sep=''), plot = last_plot(), device='pdf',
#       scale=1, width=5, height=5, unit=c("in"), dpi=300)