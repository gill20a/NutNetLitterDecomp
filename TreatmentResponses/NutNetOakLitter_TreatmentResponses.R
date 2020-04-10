# install.packages("multcomp")
library(multcomp)
library(nlme)
library(Matrix)
library(Hmisc)
library(MuMIn)
library(Polychrome)
library(maps)   
library(mapdata)
library(ggplot2)
library(ggpubr)


setwd("~/Dropbox/UMN/NutNet Litter Decay/")
Nutdata<-read.csv("20200218_Plot_curve_summary_redacted_Nutdata.csv", header=T, sep=",")
newmatNut<-cbind(Nutdata$ppm_K, Nutdata$ppm_Ca, Nutdata$ppm_Mg, Nutdata$ppm_Na)
Nutdata$BaseCat<-rowSums(newmatNut)

#######Remove extreme outlier
Nutdata[, "weibull.mrt"][Nutdata[, "weibull.mrt"] > 2000] <- NA
Nutdata$single.k<-as.numeric(as.character(Nutdata$single.k))
par(mfcol=c(1,2))
dev.off()
hist(Nutdata$weibull.mrt)

max(Nutdata[, "weibull.mrt"], na.rm=T)
# Paragraph one results section
hist(Nutdata$Prop_Init_C_Mass, xlab="Prop. Initial C Remaining", main="", col="dark grey", border=F)
abline(v=mean(Nutdata$Prop_Init_C_Mass, na.rm=T), col="dark red", lty=2, lwd=3)
abline(v=median(Nutdata$Prop_Init_C_Mass, na.rm=T), col="purple", lty=2, lwd=3)
legend("topright", c("mean", "median"), col=c("dark red", "purple"), lwd=c(3,3), bty="n" )
tapply(Nutdata$Prop_Init_C_Mass, Nutdata$Site_Text, mean, na.rm=T)


########################Distribution of Asymp A values collapsing to 0 created 2020/2/18
TbN<-table(Nutdata$N, Nutdata$Asymp_collapse)
TbP<-table(Nutdata$P, Nutdata$Asymp_collapse)
TbK<-table(Nutdata$K, Nutdata$Asymp_collapse)
TbTrt<-table(Nutdata$Trt, Nutdata$Asymp_collapse)
TbSite<-table(Nutdata$Site_Text, Nutdata$Asymp_collapse)

chisq.test(TbN)
chisq.test(TbP)
chisq.test(TbK)
chisq.test(TbTrt)
chisq.test(TbSite)

######################Histograms for predictor distributions
#############################Figure 2-5
####Figure 2
par(mfrow=c(3,3), mar=c(4,4,1,1), mgp=c(2.2,1,0))
hist((Nutdata$MAT_v2), xlab="MAT (deg C)", main="")
hist((Nutdata$MAP_v2), xlab="MAP (mm)", main="")
hist((Nutdata$pct_C), xlab="Soil Carbon (%)", main="")
hist((Nutdata$ppm_P), xlab="Soil Phosphorus (ppm)", main="")
hist((Nutdata$ppm_Mn), xlab="Soil Manganese (ppm)", main="")
hist((Nutdata$pH), xlab="Soil pH", main="")
hist((Nutdata$N_Dep), xlab="N Deposition (kg/m2/yr)", main="")
hist((Nutdata$Total_mass_mean), xlab="Aboveground Biomass (g/m2)", main="")
hist(Nutdata$Ground_PAR_mean, xlab="Ground par", main="")


#Distribution of response variables
#Figure 3
par(mfrow=c(3,3))
hist((Nutdata$single.k), main="", xlab="single k")
hist((Nutdata$weibull.alpha), main="", xlab="weibull.alpha")
hist((Nutdata$weibull.half), main="", xlab="weibull 10% mass loss")
hist((Nutdata$weibull.quarter), main="", xlab="weibull 25% mass loss")
hist((Nutdata$weibull.half.life), main="", xlab="weibull half life")
hist((Nutdata$weibull.mrt), main="", xlab="weibull mrt")
hist((Nutdata$MRT_HL), main="", xlab="MRT_HL")

####Figure 4 - Predictors after transformation
par(mfrow=c(3,3), mar=c(4,4,1,1), mgp=c(2.2,1,0))
hist((Nutdata$MAT_v2), xlab="MAT (deg C)", main="")
hist((Nutdata$MAP_v2), xlab="MAP (mm)", main="")
hist(log(Nutdata$pct_C), xlab="log(Soil Carbon) (%)", main="")
hist(log(Nutdata$ppm_P), xlab="log(Soil Phosphorus) (ppm)", main="")
hist(log(Nutdata$ppm_Mn), xlab="log(Soil Manganese) (ppm)", main="")
hist((Nutdata$pH), xlab="Soil pH", main="")
hist(log(Nutdata$N_Dep), xlab="log(N Deposition) (kg/m2/yr)", main="")
hist(sqrt(Nutdata$Total_mass_mean), xlab="sqrt(Aboveground Biomass) (g/m2)", main="")
hist((Nutdata$Ground_PAR_mean), xlab="Ground par", main="")

####Figure 5 - Responses after transformation
par(mfrow=c(3,3))
hist(log(Nutdata$single.k), main="", xlab="log(single k)")
hist(log(Nutdata$weibull.alpha), main="", xlab="log(weibull alpha)")
hist(sqrt(Nutdata$weibull.tenth), main="", xlab="sqrt(weibull 10% mass loss)")
hist((Nutdata$weibull.quarter), main="", xlab="weibull 25% mass loss")
hist(sqrt(Nutdata$weibull.half.life), main="", xlab="sqrt(weibull half life)")
hist(log(Nutdata$weibull.mrt), main="", xlab="log(weibull mrt)")
hist(log(Nutdata$MRT_HL), main="", xlab="log(MRT_HL)")

colnames(Nutdata)

###################### Figure 1 Map
long<-unique(Nutdata$longitude)
length(long)
lat<-unique(Nutdata$latitude)
length(lat)
sites<-unique(Nutdata$Site_Text)
mapdat<-cbind(long, lat, as.character(sites))
mapdat<-as.data.frame(mapdat)
mapdat2<-mapdat[c(1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,8,20), ]
#
dev.off()
col_vector=c("#4e79a7","#f28e2b", "#e15759","#76b7b2", "black", "#59a14f", "#edc948","#b07aa1","#ff9da7",
             "#9c755f", "#bab0ac","#4e79a7","#f28e2b", "#e15759","#76b7b2","black", "#59a14f", "#edc948","#b07aa1","#ff9da7")
# png(filename="map20200408.png", width = 11, height = 11, units = "cm", res=600)
map("world",boundary=T, interior=F, col="dark grey", ylim=c(-60, 90), mar=rep(0,4))
points(as.numeric(as.character(mapdat2$long)), as.numeric(as.character(mapdat2$lat)), col=col_vector,
       pch=c(rep(c(1,2,3,8),5)))
legend("left", c("Bogong", "Boulder","Bunchgrass", "Burrawan", "CBGB-Iowa", "Cedar Creek","Cowichan",
                 "Elliott", "Hall's Prarie", "Hopland", "Kinypanial", "Lookout", "McLaughlin", "Sagehen",
                 "Sedgwick", "Sheep", "Sierra", "Spindletop", "UNC", "Val Mustair"), 
       col=col_vector, pch=c(rep(c(1,2,3,8),5)),bty="n", cex=0.69)
dev.off()
dev.off(2)

########################### Table 2 Correlations
cor.mat<-cbind(abs(Nutdata$latitude), Nutdata$MAT_v2, Nutdata$MAP_v2, Nutdata$AI, Nutdata$pct_C,
               Nutdata$pct_N,  Nutdata$ppm_P,Nutdata$BaseCat, Nutdata$ppm_Mn,Nutdata$pH, Nutdata$PET, Nutdata$N_Dep, Nutdata$Total_mass_mean)
colnames(cor.mat)<-c("abs(latitude)","MAT", "MAP","AI", "pct_C", "pct_N", "ppm_P", "BaseCat","ppm_Mn", "pH", "PET", "N_dep", 
                     "total mass")
head(cor.mat)
mydata.rcorr = rcorr((cor.mat))
mydata.rcorr


### Mixed model selection
#################################################
###########Dredge Model Selection 
##            Data included in table S2
############################################################
##10% mass loss
#Remove replicates that do not have complete predictor datasets
colnames(data)
data2<-Nutdata[complete.cases(Nutdata[ ,c("weibull.tenth",
                                          "MAT_v2", 
                                          "MAP_v2",  "N_Dep","Total_mass_mean","pct_C", "ppm_Mn", 
                                          "ppm_P", "pH")]),]

model <- lme(sqrt(weibull.tenth) ~ MAT_v2+MAP_v2+log(pct_C)+log(ppm_P)+log(ppm_Mn)+
               pH+log(N_Dep)+sqrt(Total_mass_mean),random=~1|Site_Text/Block,
             data = data2, na.action = na.fail)
summary(model)
options(na.action = "na.fail") 
select_model <- dredge(model)
select_model
summary(model.avg(select_model, subset = delta <= 3))

##Random intercept model parameters for candidate models
fit1<-lme(sqrt(weibull.tenth)~1, random=~1|Site_Text, data =data2, na.action=na.fail)
summary(fit1)
r.squaredGLMM(fit1)

#####################################
#25% mass loss

data2<-Nutdata[complete.cases(Nutdata[ ,c("weibull.quarter",
                                          "MAT_v2", 
                                          "MAP_v2",  "N_Dep","total_mass","pct_C", "ppm_Mn", 
                                          "ppm_P", "pH")]),]

model <- lme((weibull.quarter) ~ MAT_v2+MAP_v2+log(pct_C)+log(ppm_P)+log(ppm_Mn)
             +pH+log(N_Dep)+sqrt(total_mass),random=~1|Site_Text/Block,
             data = data2, na.action = na.fail)
summary(model)
options(na.action = "na.fail") 
select_model <- dredge(model)
select_model
summary(model.avg(select_model, subset = delta <= 3))

fit1<-lme((weibull.quarter)~1, random=~1|Site_Text, data =data2, na.action=na.fail)
summary(fit1)
AIC(fit1)
r.squaredGLMM(fit1)

fit1<-lme((weibull.quarter)~pH, random=~1|Site_Text, data =data2, na.action=na.fail)
summary(fit1)
AIC(fit1)
r.squaredGLMM(fit1)


################~~~~~~~~~~~~~~~~~~~~~~~~~
#50% mass loss
data2<-Nutdata[complete.cases(Nutdata[ ,c("weibull.half.life",
                                          "MAT_v2", 
                                          "MAP_v2",  "N_Dep","total_mass","pct_C", "ppm_Mn", 
                                          "ppm_P", "pH")]),]

model <- lme(sqrt(weibull.half.life) ~ MAT_v2+MAP_v2+log(pct_C)+log(ppm_P)+log(ppm_Mn)
             +pH+log(N_Dep)+sqrt(total_mass),random=~1|Site_Text/Block,
             data = data2, na.action = na.fail)
summary(model)
options(na.action = "na.fail") 
select_model <- dredge(model)
select_model
summary(model.avg(select_model, subset = delta <= 3))

model1<-lme(sqrt(weibull.half.life) ~ log(N_Dep),random=~1|Site_Text,
            data = data2, na.action = na.fail)
summary(model1)
r.squaredGLMM(model1)

model2 <- lme(sqrt(weibull.half.life) ~ log(N_Dep)+MAT_v2,random=~1|Site_Text,
              data = data2, na.action = na.fail)
summary(model2)
r.squaredGLMM(model2)
#AIC of interaction model is higher and interaction not significant
model2a <- lme(sqrt(weibull.half.life) ~ log(N_Dep)*MAT_v2,random=~1|Site_Text,
               data = data2, na.action = na.fail)
summary(model2a)
r.squaredGLMM(model2a)

model3<-lme(sqrt(weibull.half.life) ~ 1,random=~1|Site_Text,
            data = data2, na.action = na.fail)
summary(model3)
r.squaredGLMM(model3)


###############################weibull mrt
data2<-Nutdata[complete.cases(Nutdata[ ,c("weibull.mrt",
                                          "MAT_v2", 
                                          "MAP_v2",  "N_Dep","total_mass","pct_C", "ppm_Mn", 
                                          "ppm_P", "pH")]),]

model <- lme(log(weibull.mrt) ~ MAT_v2+MAP_v2+log(pct_C)+log(ppm_P)+log(ppm_Mn)
             +pH+log(N_Dep)+sqrt(total_mass),random=~1|Site_Text/Block,
             data = data2, na.action = na.fail)
summary(model)
options(na.action = "na.fail") 
select_model <- dredge(model)
select_model
summary(model.avg(select_model, subset = delta <= 3))

model1<-lme(log(weibull.mrt) ~ log(N_Dep),random=~1|Site_Text,
            data = data2, na.action = na.fail)
summary(model1)
r.squaredGLMM(model1)

model2<-lme(log(weibull.mrt) ~ log(N_Dep)+MAT_v2,random=~1|Site_Text,
            data = data2, na.action = na.fail)
summary(model2)
r.squaredGLMM(model2)

model3<-lme(log(weibull.mrt) ~ log(N_Dep)+MAT_v2+log(ppm_P),random=~1|Site_Text,
            data = data2, na.action = na.fail)
summary(model3)
r.squaredGLMM(model3)
######

########### Weibull alpha
data2<-Nutdata[complete.cases(Nutdata[ ,c("weibull.alpha",
                                          "MAT_v2", 
                                          "MAP_v2",  "N_Dep","total_mass","pct_C", "ppm_Mn", 
                                          "ppm_P", "pH")]),]

model <- lme(log(weibull.alpha) ~ MAT_v2+MAP_v2+log(pct_C)+log(ppm_P)+log(ppm_Mn)
             +pH+log(N_Dep)+sqrt(total_mass),random=~1|Site_Text/Block,
             data = data2, na.action = na.fail)
summary(model)
options(na.action = "na.fail") 
select_model <- dredge(model)
select_model
summary(model.avg(select_model, subset = delta <= 3))

model1<-lme(log(weibull.alpha) ~ log(N_Dep),random=~1|Site_Text,
            data = data2, na.action = na.fail)
summary(model1)
r.squaredGLMM(model1)
model2 <- lme(log(weibull.alpha) ~ 1,random=~1|Site_Text,
              data = data2, na.action = na.fail)
summary(model2)
r.squaredGLMM(model2)


#########################single k
Nutdata$single.k<-as.numeric(as.character(Nutdata$single.k))
data2<-Nutdata[complete.cases(Nutdata[ ,c("single.k",
                                          "MAT_v2", 
                                          "MAP_v2",  "N_Dep","total_mass","pct_C", "ppm_Mn", 
                                          "ppm_P", "pH")]),]


model <- lme(log(single.k) ~ MAT_v2+MAP_v2+log(pct_C)+log(ppm_P)+log(ppm_Mn)
             +pH+log(N_Dep)+sqrt(total_mass),random=~1|Site_Text/Block,
             data = data2, na.action = na.fail)
summary(model)
options(na.action = "na.fail") 
select_model <- dredge(model)
select_model
summary(model.avg(select_model, subset = delta <= 3))

model1<-lme(log(single.k) ~ log(N_Dep),random=~1|Site_Text,
            data = data2, na.action = na.fail)
summary(model1)
r.squaredGLMM(model1)
model2 <- lme(log(single.k) ~ 1,random=~1|Site_Text,
              data = data2, na.action = na.fail)
summary(model2)
r.squaredGLMM(model2)
model3 <- lme(log(single.k) ~ log(N_Dep)+MAT_v2,random=~1|Site_Text,
              data = data2, na.action = na.fail)
summary(model3)
r.squaredGLMM(model3)


########### asymp k
dim(data2)
data2<-Nutdata[complete.cases(Nutdata[ ,c("asymp.k",
                                          "MAT_v2", 
                                          "MAP_v2",  "N_Dep","total_mass","pct_C", "ppm_Mn", 
                                          "ppm_P", "pH")]),]

model <- lme(log(asymp.k) ~ MAT_v2+MAP_v2+log(pct_C)+log(ppm_P)+log(ppm_Mn)
             +pH+log(N_Dep)+sqrt(total_mass),random=~1|Site_Text/Block,
             data = data2, na.action = na.fail)
summary(model)
options(na.action = "na.fail") 
select_model <- dredge(model)
select_model
summary(model.avg(select_model, subset = delta <= 3))

model1<-lme(log(asymp.k) ~ 1,random=~1|Site_Text,
            data = data2, na.action = na.fail)
summary(model1)
r.squaredGLMM(model1)




########### asymp A
data2<-Nutdata[complete.cases(Nutdata[ ,c("asymp.A",
                                          "MAT_v2", 
                                          "MAP_v2",  "N_Dep","total_mass","pct_C", "ppm_Mn", 
                                          "ppm_P", "pH")]),]

model <- lme((asymp.A) ~ MAT_v2+MAP_v2+log(pct_C)+log(ppm_P)+log(ppm_Mn)
             +pH+log(N_Dep)+sqrt(total_mass),random=~1|Site_Text/Block,
             data = data2, na.action = na.fail)
summary(model)
options(na.action = "na.fail") 
select_model <- dredge(model)
select_model
summary(model.avg(select_model, subset = delta <= 3))

model1<-lme((asymp.A) ~ 1,random=~1|Site_Text,
            data = data2, na.action = na.fail)
summary(model1)
r.squaredGLMM(model1)
model2 <- lme((asymp.A) ~ log(N_Dep),random=~1|Site_Text,
              data = data2, na.action = na.fail)
summary(model2)
r.squaredGLMM(model2)

########################################################################
############### Soil pH Figure

##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

col_vector=c("#4e79a7","#f28e2b", "#e15759","#76b7b2", "black", "#59a14f", "#edc948","#b07aa1","#ff9da7",
             "#bab0ac","#4e79a7", "#e15759","#76b7b2","black", "#59a14f", "#edc948","#b07aa1","#ff9da7")
lty_vector<-c(rep(1,9),rep(3,9))

data2<-Nutdata[complete.cases(Nutdata[ ,c("weibull.quarter",
                                          "MAT_v2", 
                                          "MAP_v2",  "N_Dep","total_mass","pct_C", "ppm_Mn", 
                                          "ppm_P", "pH")]),]
Form <- (data2$weibull.quarter) ~ data2$pH
M.lm <- gls((weibull.quarter) ~ pH, data = data2)
plot(Form, data = data2, xlab = "pH", ylab = "25% mass loss")
abline(M.lm, lwd = 4)
plot(M.lm) #residuals fairly evenly distributed
# 


fit1<-lme((weibull.quarter)~pH, random=~1|Site_Text, data =data2, na.action=na.fail) #random intercept model
summary(fit1)
F0 <- fitted(fit1, level = 0)
F1 <- fitted(fit1, level = 1) ## the fitted, random values for each site -  this will plot the trend line by site
I <- order(data2$pH); Polygon.Area.100s <- sort(data2$pH) 

unique(data2$Site_Text)
AllSites <- unique(data2$Site_Text) 

png(filename="SoilpH.png", width = 11, height = 11, units = "cm", res=600)
par(mfcol=c(2,2), mar=c(0,4,3,1),mgp=c(1.75, 0.5, 0))
plot(Polygon.Area.100s, F0[I], lwd = 5, type = "l", ylab = expression('Weibull '*italic(t[1/4])*''), xaxt="n",xlab = "pH", ylim = c(0.25, 3)) ## draws average trend line


## following for-loop draws sites with individual trend lines
for(i in 1:length(AllSites)) {
  x1 <- data2$pH[data2$Site_Text == AllSites[i]]
  y1 <- F1[data2$Site_Text == AllSites[i]]
  K <- order(x1)
  lines(sort(x1), y1[K], lty = lty_vector[i],col=col_vector[i], lwd=2)
}
# legend("topleft", legend=as.character(AllSites), col=col_vector[1:15], bty="n", lty=1, cex=0.75, lwd=3)
mtext("(a)",side=3, adj=0.95, line=-1.25, cex=0.75)



form<- (weibull.quarter)~pH
ctrl <- lmeControl(opt='optim');
fit2 <- lme(form, random = ~1 + pH|Site_Text, data = data2, control=ctrl, method="REML", na.action=na.fail) ## the "~1" specifies the random intecerpt and "Polygon.Area.100|fSpecies" allows the slope to vary by hummingbird species
summary(fit2)
AIC(fit1, fit2)
anova(fit1, fit2)
r.squaredGLMM(fit1)
r.squaredGLMM(fit2)

F0 <- fitted (fit2, level = 0)
F1 <- fitted(fit2, level = 1)
I <- order(data2$pH); Polygon.Area.100s <- sort(data2$pH)

par(mar=c(3,4,0,1), mgp=c(1.75, 0.5, 0))
plot(Polygon.Area.100s, F0[I], lwd = 5, type = "l", ylab = expression('Weibull '*italic(t[1/4])*''), xlab = "Soil pH", ylim = c(0.25, 3))

for(i in 1:length(AllSites)) {
  x1 <- data2$pH[data2$Site_Text == AllSites[i]]
  y1 <- F1[data2$Site_Text == AllSites[i]]
  K <- order(x1)
  lines(sort(x1), y1[K], lwd = 2, lty = lty_vector[i], col = col_vector[i])
}
mtext("(b)",side=3, adj=0.95, line=-1.25, cex=0.75)

# plot.new()
par(mar=c(0,0,3,1))
plot.new()
legend("topleft", legend=c("Bogong", "Boulder","Bunchgrass", "Burrawan", "CBGB-Iowa", "Cedar Creek","Cowichan",
                           "Elliott", "Hall's Prarie", "Kinypanial", "Lookout",  "Sagehen",
                           "Sedgwick", "Sheep", "Sierra", "Spindletop", "UNC", "Val Mustair"), col=col_vector, bty="n", lty = lty_vector, cex=0.52, lwd=3)


dev.off()

##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
My_Theme = theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 14))

Nutdata$Sqrt_weibull.half.life<-sqrt(Nutdata$weibull.half.life)
Nutdata$Log_weibull.mrt<-log(Nutdata$weibull.mrt)
Nutdata$Log_weibull.alpha<-log(Nutdata$weibull.alpha)
Nutdata$Log_single.k<-log(as.numeric(as.character(Nutdata$single.k)))
Nutdata$Log_asymp.k<-log(Nutdata$asymp.k)
Nutdata$Log_N_Dep<-log(Nutdata$N_Dep)


####
model1<-lme(Sqrt_weibull.half.life~Log_N_Dep, data=Nutdata, random=~1|Site_Text/Block,na.action = na.omit)
summary(model1)
r.squaredGLMM(model1)

P<-plot_model(model1, type = "pred",  ci.lvl = 0.95,terms = c("Log_N_Dep"), colors=c("dark grey", "red"),title = "") + My_Theme
Q<-P +  xlab("Log (N Deposition Rate)") + ylab(expression('Sqrt (Weibull '*italic(t[1/2])*')')) +
  theme(legend.position = c(0.125, 0.1), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))

HalfLife_NDep<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_blank(), axis.line = element_line(colour = "black"),
                         plot.margin = margin(0.25, 0.5, 0.25, 0.25,  "cm"),)+annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('A')), cex=5)
HalfLife_NDep



###
model1<-lme(Log_weibull.mrt~Log_N_Dep, data=Nutdata, random=~1|Site_Text/Block,na.action = na.omit)
summary(model1)
r.squaredGLMM(model1)

P<-plot_model(model1, type = "pred",  ci.lvl = 0.95,terms = c("Log_N_Dep"), colors=c("dark grey", "red"),title = "") + My_Theme
Q<-P +  xlab("Log (N Deposition Rate)") + ylab(expression('Log (Weibull '*italic(MRT)*')')) +
  theme(legend.position = c(0.125, 0.1), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))

MRT_NDep<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                    plot.margin = margin(0.25, 0.55, 0.25, 0.5,  "cm"),)+annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('B')), cex=5)
MRT_NDep
###
model1<-lme(Log_weibull.alpha~Log_N_Dep, data=Nutdata, random=~1|Site_Text/Block,na.action = na.omit)
summary(model1)
r.squaredGLMM(model1)

P<-plot_model(model1, type = "pred",  ci.lvl = 0.95,terms = c("Log_N_Dep"), colors=c("dark grey", "red"),title = "") + My_Theme
Q<-P +  xlab("Log (N Deposition Rate)") + ylab(expression('Log (Weibull '*italic(alpha)*')')) +
  theme(legend.position = c(0.125, 0.1), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))

Alpha_NDep<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"),
                      plot.margin = margin(0.25, 0.5, 0.25, 0.25,  "cm"),)+annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('C')), cex=5)
Alpha_NDep
###
model1<-lme(Log_single.k~Log_N_Dep, data=Nutdata, random=~1|Site_Text/Block,na.action = na.omit)
summary(model1)
r.squaredGLMM(model1)

P<-plot_model(model1, type = "pred",  ci.lvl = 0.95,terms = c("Log_N_Dep"), colors=c("dark grey", "red"),title = "") + My_Theme
Q<-P +  xlab("Log (N Deposition Rate)") + ylab(expression('Log (Single '*italic(k)*')')) +
  theme(legend.position = c(0.125, 0.1), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))
Single_NDep <-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"),
                        plot.margin = margin(0.25, 0.5, 0.25, 0.25,  "cm"),)+annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('D')), cex=5)
Single_NDep

###
model1<-lme(asymp.A~Log_N_Dep, data=Nutdata, random=~1|Site_Text/Block,na.action = na.omit)
summary(model1)
r.squaredGLMM(model1)

P<-plot_model(model1, type = "pred",  ci.lvl = 0.95,terms = c("Log_N_Dep"), colors=c("dark grey", "red"),title = "") + My_Theme
Q<-P +  xlab("Log (N Deposition Rate)") + ylab(expression('Asymptotic '*italic(A)*'')) +
  theme(legend.position = c(0.125, 0.1), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))

Asymp_NDep<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black"),
                      plot.margin = margin(0.25, 0.5, 0.25, 0.25,  "cm"),)+annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('E')), cex=5)
Asymp_NDep



SixPanel<-ggarrange(HalfLife_NDep, Single_NDep,
                    MRT_NDep,Asymp_NDep,
                    Alpha_NDep,NA,
                    ncol = 2, nrow = 3)



ggsave(filename = "NDepFig.pdf",
       plot = SixPanel,
       dpi=300,
       width = 18, height = 22, units = "cm")


#######################################################
### Table S5
model <- lme(sqrt(weibull.tenth)~N*P*K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(anova(model)[,"F-value"])), as.data.frame(as.matrix(anova(model)[,"p-value"])))

model <- lme((weibull.quarter)~N*P*K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(anova(model)[,"F-value"])), as.data.frame(as.matrix(anova(model)[,"p-value"])))

model <- lme(sqrt(weibull.half.life)~N*P*K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(anova(model)[,"F-value"])), as.data.frame(as.matrix(anova(model)[,"p-value"])))

model <- lme(log(weibull.mrt)~N*P*K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(anova(model)[,"F-value"])), as.data.frame(as.matrix(anova(model)[,"p-value"])))

model <- lme(log(weibull.alpha)~N*P*K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(anova(model)[,"F-value"])), as.data.frame(as.matrix(anova(model)[,"p-value"])))

model <- lme(log(single.k)~N*P*K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(anova(model)[,"F-value"])), as.data.frame(as.matrix(anova(model)[,"p-value"])))

model <- lme(log(asymp.k)~N*P*K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(anova(model)[,"F-value"])), as.data.frame(as.matrix(anova(model)[,"p-value"])))

model <- lme((asymp.A)~N*P*K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(anova(model)[,"F-value"])), as.data.frame(as.matrix(anova(model)[,"p-value"])))


########################### Figure X
# tiff(file = "20200404_EarlyResponse.tiff", width = 11, height = 11, units = "cm", res = 300)
png(filename="20200409_EarlyResponse.png", width = 11, height = 11, units = "cm", res=600)
mean<-tapply((Nutdata$weibull.tenth),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$weibull.tenth),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$weibull.tenth,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]

par(mar=c(0.5,5,0.5,1), mfcol=c(6,2), mgp=c(2,0.5,0))
plot(mat$mean, ylab=expression('Weibull '*italic(t[1/10])*''), ylim=c(min(lower),max(upper)), tck=-0.01,xaxt="n", xlab="", las=1, cex.lab=1, cex.axis=0.8)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
points(mat$mean, pch=16)
# axis(1, las=2,  at=c(1,2,3,4,5,6,7,8,9,10), labels=c("Control", "Fence", "K", "N", "NK", "NP", "NPK", "NPK+Fen.","P", "PK"), cex.axis=0.9)
mtext(side=1, text="(A)",line=-1.25, adj=0.02, cex=0.5)
box()


mean<-tapply((Nutdata$weibull.quarter),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$weibull.quarter),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$weibull.quarter,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]

plot(mat$mean, ylab=expression('Weibull '*italic(t[1/4])*''), ylim=c(min(lower),max(upper)), tck=-0.01, xaxt="n", xlab="", las=1, cex.lab=1, cex.axis=0.8 )
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
points(mat$mean, pch=16)
# axis(1, las=2,  at=c(1,2,3,4,5,6,7,8,9,10), labels=c("Control", "Fence", "K", "N", "NK", "NP", "NPK", "NPK+Fen.","P", "PK"), cex.axis=0.9)
mtext(side=1, text="(B)",line=-1.25, adj=0.02, cex=0.5)
box()

mean<-tapply((Nutdata$weibull.half.life),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$weibull.half.life),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$weibull.half.life,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]

plot(mat$mean, ylab=expression('Weibull '*italic(t[1/2])*''), ylim=c(min(lower),max(upper)), tck=-0.01,xaxt="n", xlab="", las=1, cex.lab=1, cex.axis=0.8)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
points(mat$mean, pch=16)
# axis(1, las=2,  at=c(1,2,3,4,5,6,7,8,9,10), labels=c("Control", "Fence", "K", "N", "NK", "NP", "NPK", "NPK+Fen.","P", "PK"), cex.axis=0.9)
mtext(side=1, text="(C)",line=-1.25, adj=0.02, cex=0.5)
box()


mean<-tapply((Nutdata$asymp.k),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$asymp.k),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$asymp.k,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]

plot(mat$mean, ylab=expression('Asymptotic '*italic(k)*''), ylim=c(min(lower),max(upper)), tck=-0.01, xaxt="n", xlab="", las=1, cex.lab=1, tck=-0.01, cex.axis=0.8)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
points(mat$mean, pch=16)
axis(1, las=2,  at=c(1,2,3,4,5,6,7,8), labels=c("Control", "N", "P", "K", "NP", "NK", "PK", "NPK"), cex.axis=1, tck=-0.01)
mtext(side=3, text="(D)",line=-1, adj=0.02, cex=0.5)
box()
dev.off()
dev.off(4)
##########################


########################### Figure X2
# tiff(file = "20200404_LateResponse.tiff", width = 11, height = 11, units = "cm", res = 300)
png(filename="20200404_LateResponse.png", width = 11, height = 11, units = "cm", res=600)
mean<-tapply((Nutdata$weibull.mrt),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$weibull.mrt),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$weibull.mrt,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]

par(mar=c(0.5,5,0.5,1), mfcol=c(6,2), mgp=c(2,0.5,0))
plot(mat$mean, ylab=expression('Weibull '*MRT*''), ylim=c(min(lower),max(upper)), tck=-0.01,xaxt="n", xlab="", las=1, cex.lab=1, cex.axis=0.8)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
points(mat$mean, pch=16)
# axis(1, las=2,  at=c(1,2,3,4,5,6,7,8,9,10), labels=c("Control", "Fence", "K", "N", "NK", "NP", "NPK", "NPK+Fen.","P", "PK"), cex.axis=0.9)
mtext(side=3, text="(A)",line=-1.25, adj=0.02, cex=0.5)
box()


mean<-tapply((Nutdata$weibull.alpha),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$weibull.alpha),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$weibull.alpha,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]

plot(mat$mean, ylab=expression('Weibull '*italic(alpha)*''), ylim=c(min(lower),max(upper)), tck=-0.01, xaxt="n", xlab="", las=1, cex.lab=1, cex.axis=0.8 )
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
points(mat$mean, pch=16)
# axis(1, las=2,  at=c(1,2,3,4,5,6,7,8,9,10), labels=c("Control", "Fence", "K", "N", "NK", "NP", "NPK", "NPK+Fen.","P", "PK"), cex.axis=0.9)
mtext(side=3, text="(B)",line=-1.25, adj=0.02, cex=0.5)
box()

mean<-tapply((Nutdata$single.k),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$single.k),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$single.k,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]

plot(mat$mean, ylab=expression('Single '*italic(k)*''), ylim=c(min(lower),max(upper)), tck=-0.01,xaxt="n", xlab="", las=1, cex.lab=1, cex.axis=0.8)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
points(mat$mean, pch=16)
# axis(1, las=2,  at=c(1,2,3,4,5,6,7,8,9,10), labels=c("Control", "Fence", "K", "N", "NK", "NP", "NPK", "NPK+Fen.","P", "PK"), cex.axis=0.9)
mtext(side=3, text="(C)",line=-1.25, adj=0.02, cex=0.5)
box()


mean<-tapply((Nutdata$asymp.A),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$asymp.A),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$asymp.A,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]

plot(mat$mean, ylab=expression('Asymptotic '*italic(A)*''), ylim=c(min(lower),max(upper)), tck=-0.01, xaxt="n", xlab="", las=1, cex.lab=1, tck=-0.01, cex.axis=0.8)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
points(mat$mean, pch=16)
axis(1, las=2,  at=c(1,2,3,4,5,6,7,8), labels=c("Control", "N", "P", "K", "NP", "NK", "PK", "NPK"), cex.axis=1, tck=-0.01)
mtext(side=3, text="(D)",line=-1, adj=0.02, cex=0.5)
box()
dev.off()




######Interaction Models revised 2020/1/13
#Revised again 2020/02/18 when ALG ran asymptotic data across all sites rather than just initial 3 for which the asymp best fit most sites/treatments
model <- lme(sqrt(weibull.tenth)~abs(latitude)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.quarter)~abs(latitude)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.half.life)~abs(latitude)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.mrt)~abs(latitude)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.alpha)~abs(latitude)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(single.k)~abs(latitude)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(asymp.k)~abs(latitude)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme((asymp.A)~abs(latitude)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]

##
model <- lme(sqrt(weibull.tenth)~pH+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.quarter)~pH+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.half.life)~pH+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.mrt)~pH+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.alpha)~pH+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(single.k)~pH+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(asymp.k)~pH+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme((asymp.A)~pH+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]

##
model <- lme(sqrt(weibull.tenth)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.quarter)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.half.life)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.mrt)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.alpha)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(single.k)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(asymp.k)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme((asymp.A)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]


##
model <- lme(sqrt(weibull.tenth)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.quarter)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.half.life)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.mrt)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.alpha)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(single.k)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(asymp.k)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme((asymp.A)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K+log(N_Dep):N+log(N_Dep):P++log(N_Dep):K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]

##
model <- lme(sqrt(weibull.tenth)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.quarter)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.half.life)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.mrt)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.alpha)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(single.k)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(asymp.k)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme((asymp.A)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]

##
model <- lme(sqrt(weibull.tenth)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.quarter)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.half.life)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.mrt)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.alpha)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(single.k)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(asymp.k)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme((asymp.A)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]

##
model <- lme(sqrt(weibull.tenth)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.quarter)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.half.life)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.mrt)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.alpha)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(single.k)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(asymp.k)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme((asymp.A)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]

##
model <- lme(sqrt(weibull.tenth)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.quarter)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.half.life)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.mrt)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.alpha)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(single.k)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(asymp.k)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme((asymp.A)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]

##
model <- lme(sqrt(weibull.tenth)~sqrt(Total_mass_mean)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.quarter)~sqrt(total_mass)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(sqrt(weibull.half.life)~sqrt(total_mass)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.mrt)~sqrt(Total_mass_mean)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(weibull.alpha)~sqrt(total_mass)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(single.k)~sqrt(total_mass)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme(log(asymp.k)~sqrt(total_mass)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
model <- lme((asymp.A)~sqrt(total_mass)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]





#######################################################################################
####
model <- lme(Log_single.k~Log_N_Dep*N*P*K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)

P<-plot_model(model, type = "pred",  ci.lvl = 0.95,terms = c("Log_N_Dep", "N"),title = "") + My_Theme
Q<-P +  xlab("Log(N Deposition)") + ylab("Log(Single k)") +
  theme(legend.position = c(0.125, 0.1), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))


R<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             plot.margin = margin(2, 1, 1, 1, "cm"),)+annotate("text", x = 2.5, y = -0.5, label = "*N Dep, *N, *N Dep x N", cex=5)
R
