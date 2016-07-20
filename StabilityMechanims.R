####packages required 
library(RCurl)
library(vegan)
library(ggplot2)
library(plyr)
library(gridExtra)
library(lme4)

####load data from GitHub
ZooCounts.URL<-getURL("https://raw.githubusercontent.com/Monsauce/Stability-mechanisms/master/ZooplanktonTowCounts.csv")
ZooCounts<-read.csv(text=ZooCounts.URL)

ZooDensity.URL<-getURL("https://raw.githubusercontent.com/Monsauce/Stability-mechanisms/master/ZooplanktonTowDensity.csv")
ZooDensity<-read.csv(text=ZooDensity.URL)

Isotopes.URL<-getURL("https://raw.githubusercontent.com/Monsauce/Stability-mechanisms/master/StableIsotopes.csv")
Isotopes<-read.csv(text=Isotopes.URL)

####copepod community composition 
#Run NMDS
Bray<-vegdist(ZooCounts[2:11], method="bray")
NMDS<-metaMDS(Bray, k=2)
Points<-as.data.frame(NMDS$points)

#plot NMDS
Points$Site<- "Null"
Points$Site[1:9] <- "Reference"
Points$Site[9:16] <- "Farm"

Figure1A<-ggplot(Points, aes(x = MDS1, y = MDS2))+ geom_point(aes(colour=Site,size=3))+xlab("Site")+
  ylab("NMDS 1")+ xlab("NMDS 2")+theme_minimal()+
  scale_colour_manual(values=c("black", "grey"))

#run ANOSIM
ZooANOSIM<-anosim(Bray, ZooCounts[,1])
ZooANOSIM

#density analyses 
VarianceZoo<-ddply(.data=ZooDensity, .variables=.(Treatment, Stage, Replicate), .fun= summarise, mean = mean(Density))

AverageZoo<-ddply(.data=ZooDensity, .variables=.(Treatment, Stage), .fun= summarise, mean = mean(Density), se=sd(Density)/sqrt(length(Density)))

Figure1B<-ggplot(AverageZoo, aes(x =Stage , y = mean, fill=Treatment))+geom_bar(stat = "identity",position="dodge")+xlab("Site")+
  ylab("Density (in/L)")+
  theme_minimal()+
  scale_fill_manual(values=c("black", "grey"))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,position=position_dodge(.9))

#make Figure 1 
Figure1<-grid.arrange(Figure1A,Figure1B,ncol=2)

#run ANOVA between farm and reference sites 
DensityANOVA<-aov(mean~Treatment*Stage, VarianceZoo)
Tukey<-TukeyHSD(DensityANOVA)

####stable isotopes 
#subset 1+ and 2+ data 
Mussels1<-Isotopes[Isotopes$Site=="1+" & Isotopes$Tissue=="M",] 
Mussels2<-Isotopes[Isotopes$Site=="2+" & Isotopes$Tissue=="M",] 

#run linear model 
lm.mussels1<-lm(Delta15N~Length, data=Mussels1)
lm.mussels2<-lm(Delta15N~Length, data=Mussels2)

#get confidence intervals 
Seston<-Isotopes[Isotopes$Tissue%in%c("S"),]
SestonMean<-ddply(.data=Seston, .variables=.(Site), .fun= summarise, Mean=mean(Delta15N), sd=sd(Delta15N), CI=1.96*(sd(Delta15N)/sqrt(length(Total_DW))))

Adult<-Isotopes[Isotopes$Tissue%in%c("A"),]
AdultMean<-ddply(.data=Adult, .variables=.(Site), .fun= summarise, Mean=mean(Delta15N), sd=sd(Delta15N), CI=1.96*(sd(Delta15N)/sqrt(length(Total_DW))))

Nauplii<-Isotopes[Isotopes$Tissue%in%c("N"),]
NaupliiMean<-ddply(.data=Nauplii, .variables=.(Site), .fun= summarise, Mean=mean(Delta15N), sd=sd(Delta15N), CI=1.96*(sd(Delta15N)/sqrt(length(Total_DW))))


Figure2A<-ggplot(Mussels1, aes(x = Length, y = Delta15N))+ geom_point()+ylab("∆N")+
  xlab("Length (mm)")+
  theme_minimal()+
  scale_color_manual(values=c("Black"))+
  stat_smooth(method=lm, aes(colour="Grey"))+
  geom_hline(aes(yintercept=3.400000),colour="Grey")+
  geom_hline(aes(yintercept=3.664243), colour="Grey", linetype="dashed")+#seston
  geom_hline(aes(yintercept=3.135757), colour="Grey", linetype="dashed")+
  geom_hline(aes(yintercept=5.854444), colour="Black")+#nauplii
  geom_hline(aes(yintercept=6.068801), colour="Black", linetype="dashed")+
  geom_hline(aes(yintercept=5.640087), colour="Black", linetype="dashed")+
  geom_hline(aes(yintercept=7.406667), colour="Dark Grey")+#adults
  geom_hline(aes(yintercept=7.643855), colour="Dark Grey", linetype="dashed")+
  geom_hline(aes(yintercept=7.169479), colour="Dark Grey", linetype="dashed")+
  theme(legend.position="none")+
  scale_y_continuous(limits=c(2, 8.5))


Figure2B<-ggplot(Mussels2, aes(x = Length, y = Delta15N))+ geom_point()+ylab("∆N")+
  xlab("Length (mm)")+
  theme_minimal()+
  scale_color_manual(values=c("Black"))+
  stat_smooth(method=lm, aes(colour="Grey"))+
  geom_hline(aes(yintercept=3.490667),colour="Grey")+#seston
  geom_hline(aes(yintercept=3.689138), colour="Grey", linetype="dashed")+
  geom_hline(aes(yintercept=3.292196), colour="Grey", linetype="dashed")+
  geom_hline(aes(yintercept=6.138421), colour="Black")+#nauplii
  geom_hline(aes(yintercept=6.28842), colour="Black", linetype="dashed")+
  geom_hline(aes(yintercept=5.988422), colour="Black", linetype="dashed")+
  geom_hline(aes(yintercept=7.518182), colour="Dark Grey")+#adults
  geom_hline(aes(yintercept=7.698437), colour="Dark Grey", linetype="dashed")+
  geom_hline(aes(yintercept=7.337928), colour="Dark Grey", linetype="dashed")+
  theme(legend.position="none")+
  scale_y_continuous(limits=c(2, 8.5))

#plot Figure 2
Figure2<-grid.arrange(Figure2A,Figure2B, ncol=2)

#binomial model
#subset mussels
Mussels<-Isotopes[Isotopes$Tissue%in%c("M"),]

#catergorize into binary varibles if greater than upper CI
Mussels$Binomial.Nauplii<-ifelse(Mussels$Delta15N>6.28842, 1, 0)

Mussels$Binomial.Adults<-ifelse(Mussels$Delta15N>7.337928, 1, 0)

#run model
logit.model.nauplii<-glmer(Binomial.Nauplii~Length + (1|Site), data=Mussels, family=binomial)
coef(logit.model.nauplii)

#plot Figure 3
Figure3<-ggplot(Mussels, aes(x=Length, y=Binomial.Nauplii))+geom_point()+theme_minimal()+
  geom_smooth(method = "glm", family="binomial", colour='black',se=FALSE)+
  ylab("Higher ∆N than nauplii")+
  xlab("Mussel length (mm)")







