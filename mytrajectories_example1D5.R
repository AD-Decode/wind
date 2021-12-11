library(ggplot2)
library(patternplot)
library(lme4)
library(visreg)
library(tidyr)
library(magrittr) 
library(dplyr)
library(reshape2)
library(ggpubr)


dirin<-'/Users/alex/AlexBadea_MyPapers/DavidDunson/dropbox/WindingNumber/TrackingData/Probes/'
dirout<-'/Users/alex/AlexBadea_MyPapers/DavidDunson/Figures/'
  
animal<- '190715_8' # APOE2 red
animal<- '190715_6' # APOE2 red
animal<- '190715_9' # APOE2 red
Genotype2<-'APOE22'
Genotype<-Genotype2
filein2<-paste(dirin,animal,'_Probe_D8_T1_positions.csv',sep='')
mytraj_data2 <- read.csv(filein2)
mytraj_data2 <- na.omit(mytraj_data2)
mytraj_data2['Genotype']=rep(Genotype2, length(mytraj_data2$Centre.position.X))
mytraj_data2['Animal']=rep(animal, length(mytraj_data2$Centre.position.X))


animal<-'190909_9' #APOE3 green
animal<-'190909_1' #APOE3 green
animal<-'190715_1' #APOE3 green

Genotype3<-'APOE33'
filein3<-paste(dirin,animal,'_Probe_D8_T1_positions.csv',sep='')
mytraj_data3 <- read.csv(filein3)
mytraj_data3 <- na.omit(mytraj_data3)
mytraj_data3['Genotype']=rep(Genotype3, length(mytraj_data3$Centre.position.X))
mytraj_data3['Animal']=rep(animal, length(mytraj_data3$Centre.position.X))


animal<-'191006_6' #APOE4  blue
animal<-'191006_9' #APOE4  blue
Genotype4<-'APOE44'
filein4<-paste(dirin,animal,'_Probe_D8_T1_positions.csv',sep='')
mytraj_data4 <- read.csv(filein4)
mytraj_data4 <- na.omit(mytraj_data4)
mytraj_data4['Genotype']=rep(Genotype4, length(mytraj_data4$Centre.position.X))
mytraj_data4['Animal']=rep(animal, length(mytraj_data4$Centre.position.X))

mytraj_data<-rbind(mytraj_data2, mytraj_data3,mytraj_data4)


plotD8<-ggplot(data=mytraj_data, aes(x=Centre.position.X, y=Centre.position.Y, group = Genotype, color=Genotype)) +  
  scale_color_manual(values=c('blueviolet', 'chartreuse1', 'red'))+
  geom_path (linetype=1, size=0.5)+
  theme_classic()+
  facet_grid(. ~ Genotype) + 
  coord_fixed(ratio=1)+
  #stat_smooth(method = "lm") +
  #background_grid(major = 'xy', minor = "none") + # add thin horizontal lines 
  theme_bw() +
  ggtitle(paste('Probe Day D8'))




ggsave(paste(dirout,'my3animals', 'Probe_D8.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)


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

ggsave(paste(dirout,'Example1', 'Probe_D8.pdf',sep=''), plot = last_plot(), device='pdf',
       scale=1, width=5, height=5, unit=c("in"), dpi=300)