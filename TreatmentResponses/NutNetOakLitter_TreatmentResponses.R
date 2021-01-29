urlfile = "https://raw.githubusercontent.com/gill20a/NutNetLitterDecomp/master/TreatmentResponses/NutNetNPK_OakLitterDecomp_ModelParameters.csv"
Nutdata<-read.csv(url(urlfile))
newmatNut<-cbind(Nutdata$ppm_K, Nutdata$ppm_Ca, Nutdata$ppm_Mg, Nutdata$ppm_Na)
Nutdata$BaseCat<-rowSums(newmatNut)


#######Remove extreme outlier
Nutdata[, "weibull.mrt"][Nutdata[, "weibull.mrt"] > 2000] <- NA
Nutdata$single.k<-as.numeric(as.character(Nutdata$single.k))

###################### Figure 1 Map
long<-unique(Nutdata$longitude)
length(long)
lat<-unique(Nutdata$latitude)
length(lat)
sites<-unique(Nutdata$Site_Text)
mapdat<-cbind(long, lat, as.character(sites))
mapdat<-as.data.frame(mapdat)
mapdat2<-mapdat[c(18,4,3,15,1,6,16,7,8,10,17,9,11,19,16,20,12,13,5,14), ]
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

######################Histograms for distributions
#############################
#Predictors before transformation
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
par(mfrow=c(3,3))
hist((Nutdata$single.k), main="", xlab="single k")
hist((Nutdata$weibull.alpha), main="", xlab="weibull.alpha")
hist((Nutdata$weibull.tenth), main="", xlab="weibull 10% mass loss")
hist((Nutdata$weibull.quarter), main="", xlab="weibull 25% mass loss")
hist((Nutdata$weibull.half.life), main="", xlab="weibull half life")
hist((Nutdata$weibull.mrt), main="", xlab="weibull mrt")
# hist((Nutdata$MRT_HL), main="", xlab="MRT_HL")

####Predictors after transformation
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

####Responses after transformation
par(mfrow=c(3,3))
hist(sqrt(Nutdata$single.k), main="", xlab="sqrt(single k)")
hist(log(Nutdata$weibull.alpha), main="", xlab="log(weibull alpha)")
hist(sqrt(Nutdata$weibull.tenth), main="", xlab="sqrt(weibull 10% mass loss)")
hist((Nutdata$weibull.quarter), main="", xlab="weibull 25% mass loss")
hist(sqrt(Nutdata$weibull.half.life), main="", xlab="sqrt(weibull half life)")
hist(log(Nutdata$weibull.mrt), main="", xlab="log(weibull mrt)")
# hist(log(Nutdata$MRT_HL), main="", xlab="log(MRT_HL)")

########################### Table S2 Correlations
cor.mat<-cbind(abs(Nutdata$latitude), Nutdata$MAT_v2, Nutdata$MAP_v2, Nutdata$AI, Nutdata$pct_C,
               Nutdata$pct_N,  Nutdata$ppm_P,Nutdata$BaseCat, Nutdata$ppm_Mn,Nutdata$pH, Nutdata$PET, Nutdata$N_Dep, Nutdata$Total_mass_mean)
colnames(cor.mat)<-c("abs(latitude)","MAT", "MAP","AI", "pct_C", "pct_N", "ppm_P", "BaseCat","ppm_Mn", "pH", "PET", "N_dep", 
                     "total mass")
head(cor.mat)
mydata.rcorr = rcorr((cor.mat))
mydata.rcorr

#############################################
# Paragraph two results section
hist(Nutdata$Prop_Init_C_Mass, xlab="Prop. Initial C Remaining", main="", col="dark grey", border=F)
abline(v=mean(Nutdata$Prop_Init_C_Mass, na.rm=T), col="dark red", lty=2, lwd=3)
abline(v=median(Nutdata$Prop_Init_C_Mass, na.rm=T), col="purple", lty=2, lwd=3)
legend("topright", c("mean", "median"), col=c("dark red", "purple"), lwd=c(3,3), bty="n" )
tapply(Nutdata$Prop_Init_C_Mass, Nutdata$Site_Text, mean, na.rm=T)

###############################################3
# Abstract
mean<-tapply(Nutdata$weibull.tenth, Nutdata$N, mean, na.rm=T)
1-(mean[2]/mean[1])
mean<-tapply(Nutdata$weibull.quarter, Nutdata$N, mean, na.rm=T)
1-(mean[2]/mean[1])
mean<-tapply(Nutdata$weibull.alpha, Nutdata$N, mean, na.rm=T)
1-(mean[2]/mean[1])
mean<-tapply(Nutdata$asymp.A, Nutdata$N, mean, na.rm=T)
(mean[2]/mean[1])-1
mean<-tapply((Nutdata$weibull.mrt), Nutdata$N, mean, na.rm=T)
(mean[1]/mean[2])

########################### Table S2 Correlations
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
##            Data included in table S5
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
fit1<-lme(sqrt(weibull.tenth)~1, random=~1|Site_Text/Block, data =data2, na.action=na.fail)
summary(fit1)
AIC(fit1)
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

fit1<-lme((weibull.quarter)~1, random=~1|Site_Text/Block, data =data2, na.action=na.fail)
summary(fit1)
AIC(fit1)
r.squaredGLMM(fit1)

fit1<-lme((weibull.quarter)~pH, random=~1|Site_Text/Block, data =data2, na.action=na.fail)
summary(fit1)
AIC(fit1)
r.squaredGLMM(fit1)

## random slope and intercept
fit1<-lme((weibull.quarter)~pH, random=~1 + pH|Site_Text, data =data2, na.action=na.fail)
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

model2 <- lme(sqrt(weibull.half.life) ~ log(N_Dep)+MAT_v2,random=~1|Site_Text/Block,
              data = data2, na.action = na.fail)
summary(model2)
r.squaredGLMM(model2)
AIC(model)
#AIC of interaction model is higher and interaction not significant
model2a <- lme(sqrt(weibull.half.life) ~ log(N_Dep)*MAT_v2,random=~1|Site_Text/Block,
               data = data2, na.action = na.fail)
summary(model2a)
r.squaredGLMM(model2a)

model1<-lme(sqrt(weibull.half.life) ~ log(N_Dep),random=~1|Site_Text/Block,
            data = data2, na.action = na.fail)
summary(model1)
r.squaredGLMM(model1)


model3<-lme(sqrt(weibull.half.life) ~ 1,random=~1|Site_Text/Block,
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

model1<-lme(log(weibull.mrt) ~ log(N_Dep),random=~1|Site_Text/Block,
            data = data2, na.action = na.fail)
summary(model1)
r.squaredGLMM(model1)

model2<-lme(log(weibull.mrt) ~ log(N_Dep)+MAT_v2,random=~1|Site_Text/Block,
            data = data2, na.action = na.fail)
summary(model2)
r.squaredGLMM(model2)

model3<-lme(log(weibull.mrt) ~ log(N_Dep)+MAT_v2+log(ppm_P),random=~1|Site_Text/Block,
            data = data2, na.action = na.fail)
summary(model3)


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

model1<-lme(log(weibull.alpha) ~ log(N_Dep),random=~1|Site_Text/Block,
            data = data2, na.action = na.fail)
summary(model1)
r.squaredGLMM(model1)

model2 <- lme(log(weibull.alpha) ~ 1,random=~1|Site_Text/Block,
              data = data2, na.action = na.fail)
summary(model2)
r.squaredGLMM(model2)

model3<-lme(log(weibull.alpha) ~ log(N_Dep)+pH,random=~1|Site_Text/Block,
            data = data2, na.action = na.fail)
summary(model3)
r.squaredGLMM(model3)


#########################single k
Nutdata$single.k<-as.numeric(as.character(Nutdata$single.k))
data2<-Nutdata[complete.cases(Nutdata[ ,c("single.k",
                                          "MAT_v2", 
                                          "MAP_v2",  "N_Dep","total_mass","pct_C", "ppm_Mn", 
                                          "ppm_P", "pH")]),]


model <- lme(sqrt(single.k) ~ MAT_v2+MAP_v2+log(pct_C)+log(ppm_P)+log(ppm_Mn)
             +pH+log(N_Dep)+sqrt(total_mass),random=~1|Site_Text/Block,
             data = data2, na.action = na.fail)
summary(model)
options(na.action = "na.fail") 
select_model <- dredge(model)
select_model
summary(model.avg(select_model, subset = delta <= 3))

model1<-lme(sqrt(single.k) ~ log(N_Dep),random=~1|Site_Text/Block,
            data = data2, na.action = na.fail)
summary(model1)
r.squaredGLMM(model1)
model2 <- lme(sqrt(single.k) ~ 1,random=~1|Site_Text/Block,
              data = data2, na.action = na.fail)
summary(model2)
r.squaredGLMM(model2)



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

model1<-lme(log(asymp.k) ~ 1,random=~1|Site_Text/Block,
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

model1<-lme((asymp.A) ~ 1,random=~1|Site_Text/Block,
            data = data2, na.action = na.fail)
summary(model1)
r.squaredGLMM(model1)
model2 <- lme((asymp.A) ~ log(N_Dep),random=~1|Site_Text/Block,
              data = data2, na.action = na.fail)
summary(model2)
r.squaredGLMM(model2)


########### Nmax
data2<-Nutdata[complete.cases(Nutdata[ ,c("Nmax",
                                          "MAT_v2", 
                                          "MAP_v2",  "N_Dep","total_mass","pct_C", "ppm_Mn", 
                                          "ppm_P", "pH")]),]
hist(log(Nutdata$Nmax))
model <- lme(log(Nmax) ~ MAT_v2+MAP_v2+log(pct_C)+log(ppm_P)+log(ppm_Mn)
             +pH+log(N_Dep)+sqrt(total_mass),random=~1|Site_Text/Block,
             data = data2, na.action = na.fail)
summary(model)
options(na.action = "na.fail") 
select_model <- dredge(model)
select_model
summary(model.avg(select_model, subset = delta <= 3))

model1<-lme((Nmax) ~ 1,random=~1|Site_Text/Block,
            data = data2, na.action = na.fail)
summary(model1)
r.squaredGLMM(model1)
model2 <- lme(log(Nmax) ~ log(ppm_Mn),random=~1|Site_Text/Block,
              data = data2, na.action = na.fail)
summary(model2)
r.squaredGLMM(model2)

########### T_NNmax
hist(sqrt(Nutdata$T_nmax))
data2<-Nutdata[complete.cases(Nutdata[ ,c("T_nmax",
                                          "MAT_v2", 
                                          "MAP_v2",  "N_Dep","total_mass","pct_C", "ppm_Mn", 
                                          "ppm_P", "pH")]),]
hist((Nutdata$Nmax))
model <- lme(sqrt(T_nmax) ~ MAT_v2+MAP_v2+log(pct_C)+log(ppm_P)+log(ppm_Mn)
             +pH+log(N_Dep)+sqrt(total_mass),random=~1|Site_Text/Block,
             data = data2, na.action = na.fail)
summary(model)
options(na.action = "na.fail") 
select_model <- dredge(model)
select_model
summary(model.avg(select_model, subset = delta <= 3))



### Figure 2
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
My_Theme = theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 14))

Nutdata$Sqrt_weibull.half.life<-sqrt(Nutdata$weibull.half.life)
Nutdata$Log_weibull.mrt<-log(Nutdata$weibull.mrt)
Nutdata$Log_weibull.alpha<-log(Nutdata$weibull.alpha)
Nutdata$Sqrt_single.k<-sqrt(as.numeric(as.character(Nutdata$single.k)))
Nutdata$Log_asymp.k<-log(Nutdata$asymp.k)
Nutdata$Log_N_Dep<-log(Nutdata$N_Dep)
###########################################
My_Theme = theme(
  axis.title.x = element_text(size = 12),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 14))
model1<-lme(Sqrt_weibull.half.life~Log_N_Dep, data=Nutdata, random=~1|Site_Text/Block,na.action = na.omit)
summary(model1)
r.squaredGLMM(model1)
P<-plot_model(model1, type = "pred",  ci.lvl = 0.95,terms = c("Log_N_Dep"), colors=c("dark grey", "red"),title = "") + My_Theme
Q<-P + xlab(bquote('Log (N Dep.)'*~'(kg N '~ha^-1 ~year^-1*')')) + ylab(expression('Sqrt (Weibull '*italic(t[1/2])*') (years)')) +
  theme(legend.position = c(0.125, 0.1), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid")) 
 

HalfLife_NDep<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                         plot.margin = margin(0.25, 0.5, 0.25, 0.25,  "cm"),)+annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('a')), cex=5)
HalfLife_NDep



###
model1<-lme(Log_weibull.mrt~Log_N_Dep, data=Nutdata, random=~1|Site_Text/Block,na.action = na.omit)
summary(model1)
r.squaredGLMM(model1)

P<-plot_model(model1, type = "pred",  ci.lvl = 0.95,terms = c("Log_N_Dep"), colors=c("dark grey", "red"),title = "") + My_Theme
Q<-P +  xlab(bquote('Log (N Dep.)'*~'(kg N '~ha^-1 ~year^-1*')'))  + ylab(expression('Log (Weibull '*italic(MRT)*') (years)')) +
  theme(legend.position = c(0.125, 0.1), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))

MRT_NDep<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                    plot.margin = margin(0.25, 0.55, 0.25, 0.5,  "cm"),)+annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('b')), cex=5)
MRT_NDep
###
model1<-lme(Log_weibull.alpha~Log_N_Dep, data=Nutdata, random=~1|Site_Text/Block,na.action = na.omit)
summary(model1)
r.squaredGLMM(model1)

P<-plot_model(model1, type = "pred",  ci.lvl = 0.95,terms = c("Log_N_Dep"), colors=c("dark grey", "red"),title = "") + My_Theme
Q<-P +  xlab(bquote('Log (N Dep.)'*~'(kg N '~ha^-1 ~year^-1*')')) + ylab(expression('Log (Weibull '*italic(alpha)*')')) +
  theme(legend.position = c(0.125, 0.1), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))

Alpha_NDep<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                      plot.margin = margin(0.25, 0.5, 0.25, 0.25,  "cm"),)+annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('c')), cex=5)
Alpha_NDep
###
model1<-lme(Sqrt_single.k~Log_N_Dep, data=Nutdata, random=~1|Site_Text/Block,na.action = na.omit)
summary(model1)
r.squaredGLMM(model1)

P<-plot_model(model1, type = "pred",  ci.lvl = 0.95,terms = c("Log_N_Dep"), colors=c("dark grey", "red"),title = "") + My_Theme
Q<-P +   xlab(bquote('Log (N Dep.)'*~'(kg N '~ha^-1 ~year^-1*')')) + ylab(expression('Sqrt ('*italic(k[s])*')')) +
  theme(legend.position = c(0.125, 0.1), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))
Single_NDep <-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                        plot.margin = margin(0.25, 0.5, 0.25, 0.25,  "cm"),)+annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('d')), cex=5)
Single_NDep

###
model1<-lme(asymp.A~Log_N_Dep, data=Nutdata, random=~1|Site_Text/Block,na.action = na.omit)
summary(model1)
r.squaredGLMM(model1)

P<-plot_model(model1, type = "pred",  ci.lvl = 0.95,terms = c("Log_N_Dep"), colors=c("dark grey", "red"),title = "") + My_Theme
Q<-P +   xlab(bquote('Log (N Dep.)'*~'(kg N '~ha^-1 ~year^-1*')')) + ylab(expression('Asymptotic '*italic(A)*'')) +
  theme(legend.position = c(0.125, 0.1), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))

Asymp_NDep<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                      plot.margin = margin(0.25, 0.5, 0.25, 0.25,  "cm"),)+annotate("text",colour="black",x = Inf, y = Inf, hjust=2, vjust=1.75, label = expression(bold('e')), cex=5)
Asymp_NDep



SixPanel<-ggarrange(HalfLife_NDep, MRT_NDep, Alpha_NDep,
                    Single_NDep,Asymp_NDep,NA,
                    ncol = 3, nrow = 2)



ggsave(filename = "20210126_2_NDepFig.pdf",
       plot = SixPanel,
       dpi=300,
       width = 18, height = 18, units = "cm")


#######################################################
### Table S6
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

model <- lme(sqrt(single.k)~N*P*K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
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


########################### Figure 3
# tiff(file = "20200404_EarlyResponse.tiff", width = 11, height = 11, units = "cm", res = 300)
# png(filename="20200409_EarlyResponse.png", width = 11, height = 11, units = "cm", res=600)
mean<-tapply((Nutdata$weibull.tenth),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$weibull.tenth),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$weibull.tenth,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]

par(mar=c(0.5,6,2.5,1), mfcol=c(3,2), mgp=c(3,0.5,0))
plot(mat$mean, ylab=expression(''*italic(t[1/10])*'  (years)'), ylim=c(min(lower),max(upper)), tck=-0.01,xaxt="n", xlab="", las=1, cex.lab=2, cex.axis=1.5)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
circle<-rep(.3,8)
symbols(x=position, mat$mean, circles=circle,inches=1/14,ann=F, bg=c("black"), fg=NULL, add=T)
# axis(1, las=2,  at=c(1,2,3,4,5,6,7,8,9,10), labels=c("Control", "Fence", "K", "N", "NK", "NP", "NPK", "NPK+Fen.","P", "PK"), cex.axis=0.9)
mtext(side=3, text="a",line=-1.5, adj=0.02, cex=1, font=2)
box()


mean<-tapply((Nutdata$weibull.quarter),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$weibull.quarter),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$weibull.quarter,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]

par(mar=c(0.5,6,0.5,1))
plot(mat$mean, ylab=expression(''*italic(t[1/4])*'  (years)'), ylim=c(min(lower),max(upper)), tck=-0.01, xaxt="n", xlab="", las=1, cex.lab=2, cex.axis=1.5 )
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
circle<-rep(.3,8)
symbols(x=position, mat$mean, circles=circle,inches=1/14,ann=F, bg=c("black"), fg=NULL, add=T)
# axis(1, las=2,  at=c(1,2,3,4,5,6,7,8,9,10), labels=c("Control", "Fence", "K", "N", "NK", "NP", "NPK", "NPK+Fen.","P", "PK"), cex.axis=0.9)
mtext(side=3, text="b",line=-1.5, adj=0.02, cex=1, font=2)
box()


mean<-tapply((Nutdata$asymp.k),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$asymp.k),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$asymp.k,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]
par(mar=c(4.5,6,0.5,1))
plot(mat$mean, ylab=expression(''*italic(k[a])~' ('~years^-1*')'),xlab="Treatment", ylim=c(min(lower),max(upper)), tck=-0.01, xaxt="n",las=1, cex.lab=2, tck=-0.01, cex.axis=1.5)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
circle<-rep(.3,8)
symbols(x=position, mat$mean, circles=circle,inches=1/14,ann=F, bg=c("black"), fg=NULL, add=T)
axis(1, las=1,  at=c(1,2,3,4,5,6,7,8), labels=c("Control", "N", "P", "K", "NP", "NK", "PK", "NPK"), cex.axis=1.5, tck=-0.01)
mtext(side=3, text="c",line=-1.5, adj=0.02, cex=1, font=2)
box()

##########################


########################### Figure 4
# tiff(file = "20200404_LateResponse.tiff", width = 11, height = 11, units = "cm", res = 300)
# png(filename="20200404_LateResponse.png", width = 11, height = 11, units = "cm", res=600)


mean<-tapply((Nutdata$single.k),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$single.k),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$single.k,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]

par(mar=c(0.5,5.5,5.5,1), mfcol=c(2,2), mgp=c(3,0.5,0))
plot(mat$mean, ylab=expression(''*italic(k[s])~' ('~years^-1*')'), ylim=c(min(lower),max(upper)), tck=-0.01,xaxt="n", xlab="", las=1, cex.lab=2, cex.axis=1)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
circle<-rep(.3,8)
symbols(x=position, mat$mean, circles=circle,inches=1/14,ann=F, bg=c("black"), fg=NULL, add=T)
# axis(1, las=2,  at=c(1,2,3,4,5,6,7,8,9,10), labels=c("Control", "Fence", "K", "N", "NK", "NP", "NPK", "NPK+Fen.","P", "PK"), cex.axis=0.9)
mtext(side=3, text="a",line=-1.5, adj=0.02, cex=1, font=2)
box()

mean<-tapply((Nutdata$asymp.A),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$asymp.A),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$asymp.A,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]

par(mar=c(5.5,5.5,0.5,1))
plot(mat$mean, ylab=expression('Asymptoitic '*italic(A)*''), ylim=c(min(lower),max(upper)), tck=-0.01, xaxt="n", xlab="Treatment", las=1, cex.lab=2, tck=-0.01, cex.axis=1.5)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
circle<-rep(.3,8)
symbols(x=position, mat$mean, circles=circle,inches=1/14,ann=F, bg=c("black"), fg=NULL, add=T)
axis(1, las=1,  at=c(1,2,3,4,5,6,7,8), labels=c("Control", "N", "P", "K", "NP", "NK", "PK", "NPK"), cex.axis=1.5, tck=-0.01)
mtext(side=3, text="b",line=-1.5, adj=0.02, cex=1, font=2)
box()


mean<-tapply((Nutdata$weibull.mrt),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$weibull.mrt),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$weibull.mrt,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]

par(mar=c(0.5,5.5,5.5,1)) 
plot(mat$mean, ylab=expression(''*italic(MRT)~' ('~years*')'), ylim=c(min(lower),max(upper)), tck=-0.01,xaxt="n", xlab="", las=1, cex.lab=2, cex.axis=1.5)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
circle<-rep(.3,8)
symbols(x=position, mat$mean, circles=circle,inches=1/14,ann=F, bg=c("black"), fg=NULL, add=T)
# axis(1, las=2,  at=c(1,2,3,4,5,6,7,8,9,10), labels=c("Control", "Fence", "K", "N", "NK", "NP", "NPK", "NPK+Fen.","P", "PK"), cex.axis=0.9)
mtext(side=3, text="c",line=-1.5, adj=0.02, cex=1, font=2)
box()


mean<-tapply((Nutdata$weibull.alpha),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$weibull.alpha),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$weibull.alpha,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]

par(mar=c(5.5,5.5,0.5,1)) 
plot(mat$mean, ylab=expression(alpha), ylim=c(min(lower),max(upper)), tck=-0.01, xaxt="n", xlab="Treatment", las=1, cex.lab=2, cex.axis=1.5 )
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
circle<-rep(.3,8)
symbols(x=position, mat$mean, circles=circle,inches=1/14,ann=F, bg=c("black"), fg=NULL, add=T)
# axis(1, las=2,  at=c(1,2,3,4,5,6,7,8,9,10), labels=c("Control", "Fence", "K", "N", "NK", "NP", "NPK", "NPK+Fen.","P", "PK"), cex.axis=0.9)
axis(1, las=1,  at=c(1,2,3,4,5,6,7,8), labels=c("Control", "N", "P", "K", "NP", "NK", "PK", "NPK"), cex.axis=1.5, tck=-0.01)
mtext(side=3, text="d",line=-1.5, adj=0.02, cex=1, font=2)
box()
dev.off()

##################### Figure 5
######################################################################
mean<-tapply(Nutdata$single.k, Nutdata$N, mean, na.rm=T )
t1<-seq(0,10,length.out=10000)
par(mfrow=c(2,3), oma=c(0,0,0,0), mar=c(5,5.5,2,1), mgp=c(3,0.8,0), cex.lab=2,  cex.main=2, cex.axis=1.5)
MLdata<-read.csv("5.23.2019 Nutnet Masterlist.csv",  header=T, sep=",")
head(MLdata)
plot(MLdata$Yrs_Since_Deployment, MLdata$Prop_Init_C_Mass, pch=1, col="grey87",xlab="Time (years)",ylab="Prop. C Remaining", main="Single Exponential", xlim=c(0,10), ylim=c(0,1))
lines(t1,exp(-mean[1]*t1), col="grey50", lty=2, type="l", lwd=2)
lines(t1,exp(-mean[2]*t1), col="forest green", lty=1, lwd=2)
legend("topright", c("Control","+Nitrogen"), col=c("grey50","forest green"), lty=c(1,1), bty="n", cex=1.5, lwd=c(2,2))
mtext(side=1, text="a",line=-1.4, adj=0.02, cex=1, font=2)

mean_ka<-tapply(Nutdata$asymp.k, Nutdata$N, mean, na.rm=T )
mean_A<-tapply(Nutdata$asymp.A, Nutdata$N, mean, na.rm=T )
plot(MLdata$Yrs_Since_Deployment, MLdata$Prop_Init_C_Mass, pch=16, col="grey87", xlab="Time (years)",ylab="Prop. C Remaining", main="Asymptotic Exponential", xlim=c(0,10), ylim=c(0,1))
lines(t1,mean_A[1]+((1-mean_A[1])*exp(-mean_ka[1]*t1)), col="grey50", lty=1)
lines(t1,mean_A[2]+((1-mean_A[2])*exp(-mean_ka[2]*t1)), col="forest green", lty=2)
legend("topright", c("Control","+Nitrogen"), col=c("grey50","forest green"), lty=c(1,1), bty="n", cex=1.5, lwd=c(2,2))
mtext(side=1, text="b",line=-1.4, adj=0.02, cex=1, font=2)

mean_alpha<-tapply(Nutdata$weibull.alpha, Nutdata$N, mean, na.rm=T )
mean_beta<-tapply(Nutdata$weibull.beta, Nutdata$N, mean, na.rm=T )
weibull.fit<-exp(-(t1/mean_beta[1])^mean_alpha[1])
plot(MLdata$Yrs_Since_Deployment, MLdata$Prop_Init_C_Mass, pch=16, col="grey87",xlab="Time (years)",ylab="Prop. C Remaining", main="Weibull", xlim=c(0,10), ylim=c(0,1))
lines(t1,exp(-(t1/mean_beta[1])^mean_alpha[1]), col="grey50", lty=1)
lines(t1,exp(-(t1/mean_beta[2])^mean_alpha[2]), col="forest green", lty=2)
legend("topright", c("Control","+Nitrogen"), col=c("grey50","forest green"), lty=c(1,1), bty="n", cex=1.25, lwd=c(2,2))
mtext(side=1, text="c",line=-1.4, adj=0.02, cex=1, font=2)


###### ANCOVAS Table  S7
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
model <- lme(sqrt(single.k)~abs(latitude)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
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
model <- lme(sqrt(single.k)~pH+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
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
model <- lme(sqrt(single.k)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
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
model <- lme(sqrt(single.k)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
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
model <- lme(sqrt(single.k)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
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
model <- lme(sqrt(single.k)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
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
model <- lme(sqrt(single.k)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
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
model <- lme(sqrt(single.k)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
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
model <- lme(sqrt(single.k)~sqrt(total_mass)+N+P+K+N:P+N:K+P:K+N:P:K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
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
Nutdata$Sqrt_single.k<-sqrt(Nutdata$single.k)
model <- lme(Sqrt_single.k~Log_N_Dep*N*P*K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)

My_Theme = theme(
  axis.title.x = element_text(size = 18),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 18))

P<-plot_model(model, type = "pred",  ci.lvl = 0.95,terms = c("Log_N_Dep", "N"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P +  xlab(bquote('Sqrt (N Deposition)'*~'(kg N '~ha^-1 ~year^-1*')')) +ylab(bquote('Log (Single '*italic(k)~') (' ~years^-1*')'))+
  theme(legend.position = c(0.1, 0.9), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))


R<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
             plot.margin = margin(2, 1, 1, 1, "cm"),)

R

############################### Nitrogen content figures
######## N results figure

hist(log(Nutdata$Nmax))
Nutdata$Log_Nmax<-log(Nutdata$Nmax)
model <- lme((Log_Nmax)~N*P*K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(anova(model)[,"F-value"])), as.data.frame(as.matrix(anova(model)[,"p-value"])))


model <- lme((T_nmax)~N*P*K,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(anova(model)[,"F-value"])), as.data.frame(as.matrix(anova(model)[,"p-value"])))


par(mfcol=c(2,2),mar=c(0.5,5,3.5,1))
boxplot(log(Nutdata$Nmax)~Nutdata$N, na.action = na.omit, xaxt="n", xlab="", ylim=c(-3,-1),col=c("grey50", "forest green"), ylab="Log(Max litter N pool) (g)")
  mtext(side=3, text="a",line=-1.5, adj=0.05, cex=1.5, font=2)
mtext(side=3, text="p < 0.0001",line=-1.5, adj=0.98, cex=1)

par(mar=c(3.5,5,0.5,1))
boxplot(Nutdata$T_nmax~Nutdata$N, na.action = na.omit, ylim=c(0,8),names=c("Control", "+Nitrogen"), xlab="Treatment", col=c( "grey50", "forest green"), ylab="Time to max litter N (years)")
mtext(side=3, text="b",line=-1.5, adj=0.05, cex=1.5, font=2)
mtext(side=3, text="p = 0.031",line=-1.5, adj=0.98, cex=1)

par(mar=c(0.5,1,3.5,5))
boxplot(log(Nutdata$Nmax)~Nutdata$P, na.action = na.omit, xaxt="n",xlab="",  ylim=c(-3,-1),col=c("grey50", "#8b228b"), ylab="", yaxt="n")
mtext(side=3, text="c",line=-1.5, adj=0.05, cex=1.5, font=2)
mtext(side=3, text="p = 0.0643",line=-1.5, adj=0.98, cex=1)

par(mar=c(3.5,1,0.5,5))
boxplot(Nutdata$T_nmax~Nutdata$P, na.action = na.omit, ylim=c(0,8),names=c("Control", "+Phosphorus"), xlab="Treatment",  col=c("grey50", "#8b228b"), ylab="", yaxt="n")
mtext(side=3, text="d",line=-1.5, adj=0.05, cex=1.5, font=2)
mtext(side=3, text="p = 0.12",line=-1.5, adj=0.98, cex=1)

