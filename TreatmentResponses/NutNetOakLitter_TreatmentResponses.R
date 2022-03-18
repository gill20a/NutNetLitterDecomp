urlfile = "https://raw.githubusercontent.com/gill20a/NutNetLitterDecomp/master/TreatmentResponses/NutNetNPK_OakLitterDecomp_ModelParameters.csv"
data<-read.csv(url(urlfile))

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
library(sjPlot)

newmatNut<-cbind(Nutdata$ppm_K, Nutdata$ppm_Ca, Nutdata$ppm_Mg, Nutdata$ppm_Na)
Nutdata$BaseCat<-rowSums(newmatNut)
Nutdata$MAPDist<-Nutdata$MAP_WET_M_v2/Nutdata$MAP_v2

#######Remove extreme outlier
Nutdata[, "weibull.mrt"][Nutdata[, "weibull.mrt"] > 2000] <- NA
Nutdata$single.k<-as.numeric(as.character(Nutdata$single.k))

##########################
TbN<-table(Nutdata$C_Length)
TbN/length((Nutdata$C_Length))
TbN<-table(Nutdata$N_Length)

######################### Table S1
soilC<-tapply(Nutdata$pct_C, Nutdata$site_name, mean, na.rm=T)
soilN<-tapply(Nutdata$pct_N, Nutdata$site_name, mean, na.rm=T)
ppmMn<-tapply(Nutdata$ppm_Mn, Nutdata$site_name, mean, na.rm=T)
MAPDist<-tapply(Nutdata$MAPDist, Nutdata$site_name, mean, na.rm=T)
BaseCat<-tapply(Nutdata$BaseCat, Nutdata$site_name, mean, na.rm=T)
TotalMass<-tapply(Nutdata$Total_mass_mean, Nutdata$site_name, mean, na.rm=T)
MAPDist<-tapply(Nutdata$MAPDist, Nutdata$site_name, mean, na.rm=T)
GroundPAR<-tapply(Nutdata$Ground_PAR_mean, Nutdata$site_name, mean, na.rm=T)

###################### Figure 1 Map
long<-unique(Nutdata$longitude)
length(long)
lat<-unique(Nutdata$latitude)
length(lat)
sites<-unique(Nutdata$Site_Text)
mapdat<-cbind(long, lat, as.character(sites))
mapdat<-as.data.frame(mapdat)
mapdat2<-mapdat[c(9,11,16,3,5,1,15,12,7,17,8,20,19,4,18,13,10,2,6,14), ]

col_vector=c("#4e79a7","#f28e2b", "#e15759","#76b7b2", "black", "#59a14f", "#edc948","#b07aa1","#ff9da7",
             "#9c755f", "#bab0ac","#4e79a7","#f28e2b", "#e15759","#76b7b2","black", "#59a14f", "#edc948","#b07aa1","#ff9da7")
# png(filename="map20200408.png", width = 11, height = 11, units = "cm", res=600)


l <- layout(matrix(c(1, 1, 1,2,3, 4), nrow=2, byrow=T))
par(mar=c(0.5,0,0,0), mgp=c(2.5,0.8,0), oma=c(0,1,0,1))
map("world",boundary=T, interior=F, col="dark grey", ylim=c(-60, 90), mar=rep(0,4))
points(as.numeric(as.character(mapdat2$long)), as.numeric(as.character(mapdat2$lat)), col=col_vector, cex=2,lwd=2, 
       pch=c(rep(c(1,2,3,8),5)), font=2)
legend("left", c("Bogong", "Boulder","Bunchgrass", "Burrawan", "CBGB-Iowa", "Cedar Creek","Cowichan",
                 "Elliott", "Hall's Prarie", "Hopland", "Kinypanial", "Lookout", "McLaughlin", "Sagehen",
                 "Sedgwick", "Sheep", "Sierra", "Spindletop", "UNC", "Val Mustair"),
       col=col_vector, pch=c(rep(c(1,2,3,8),5)),bty="n", cex=1)
mtext(c('a'), side=3, line=-2.5,adj=0.98, cex=1.25, font=2)
box()

par(mar=c(4,4,0.5,1), mgp=c(2.5,0.8,0))
plot(tapply(Nutdata$MAT_v2, Nutdata$Site_Text, mean, na.rm=T), 
     tapply(Nutdata$MAP_v2, Nutdata$Site_Text, mean, na.rm=T),
     col=col_vector, pch=c(rep(c(1,2,3,8),5)), xlab=c("MAT (\u00B0C)"),
     ylab=c("MAP (mm)"), cex=2, cex.lab=1.75, cex.axis=1.25, lwd=2)
fit1<-lm(tapply(Nutdata$MAP_v2, Nutdata$Site_Text, mean, na.rm=T)~tapply(Nutdata$MAT_v2, Nutdata$Site_Text, mean, na.rm=T))
summary(fit1)
mtext(c('b'), side=3, line=-2.25,adj=0.05, cex=1.25, font=2)

plot(tapply(Nutdata$MAT_v2, Nutdata$Site_Text, mean, na.rm=T), 
     tapply(Nutdata$MAPDist, Nutdata$Site_Text, mean, na.rm=T),
     col=col_vector, pch=c(rep(c(1,2,3,8),5)), xlab=c("MAT (\u00B0C)"),
     ylab=c("Precipitation Distribution"), cex=2, cex.lab=1.75, cex.axis=1.25, lwd=2)
fit1<-lm(tapply(Nutdata$MAP_v2, Nutdata$Site_Text, mean, na.rm=T)~tapply(Nutdata$MAT_v2, Nutdata$Site_Text, mean, na.rm=T))
summary(fit1)
mtext(c('c'), side=3, line=-2.25,adj=0.05, cex=1.25, font=2)

par(mar=c(4,4,0.5,0))
plot(tapply(Nutdata$MAP_v2, Nutdata$Site_Text, mean, na.rm=T), 
     tapply(Nutdata$MAPDist, Nutdata$Site_Text, mean, na.rm=T),
     col=col_vector, pch=c(rep(c(1,2,3,8),5)), xlab=c("MAP (mm)"),
     ylab=c("Precipitation Distribution"), cex=2, cex.lab=1.75, cex.axis=1.25, lwd=2)
fit1<-lm(tapply(Nutdata$MAP_v2, Nutdata$Site_Text, mean, na.rm=T)~tapply(Nutdata$MAT_v2, Nutdata$Site_Text, mean, na.rm=T))
summary(fit1)
mtext(c('d'), side=3, line=-2.25,adj=0.95, cex=1.25, font=2)

################ Correlation Matrix of parameters
Nutdata$Sqrt_weibull.tenth<-sqrt(Nutdata$weibull.tenth)
Nutdata$Sqrt_weibull.quarter<-sqrt(Nutdata$weibull.quarter)
Nutdata$Sqrt_weibull.half.life<-sqrt(Nutdata$weibull.half.life)
Nutdata$Log_weibull.mrt<-log(Nutdata$weibull.mrt)
Nutdata$Log_weibull.alpha<-log(Nutdata$weibull.alpha)
Nutdata$Sqrt_single.k<-sqrt(Nutdata$single.k)
Nutdata$Log_asymp.k<-log(Nutdata$asymp.k)

parameter_correlation<-cbind(Nutdata$Sqrt_weibull.tenth,Nutdata$Sqrt_weibull.quarter,Nutdata$Sqrt_weibull.half.life,
                             Nutdata$Log_weibull.mrt, Nutdata$Log_weibull.alpha, Nutdata$Sqrt_single.k,
                             Nutdata$Log_asymp.k, Nutdata$asymp.A)
colnames(parameter_correlation)<-c("Sqrt(Weibull t1/10)", "Sqrt(Weibull t1/4)", "Sqrt(Weibull t1/2)",
                                   "Log(Weibull MRT)", "Log(Weibull alpha)", "Sqrt(Single k)",
                                   "Log(Asymp k)", "Asymp A")

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


########################### Figure 2
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
plot(mat$mean, ylab=expression(''*italic(t[1/10])*'  (years)'), ylim=c(min(lower),max(upper)), yaxt="n", tck=-0.01,xaxt="n", xlab="", las=1, cex.lab=2, cex.axis=1.5)
axis(side = 2, at = seq(0.4, 0.8, by = 0.05), las=2, cex.axis=1.4, tck=-0.01)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
circle<-rep(.3,8)
symbols(x=position, mat$mean, circles=circle,inches=1/14,ann=F, bg=c("black"), fg=NULL, add=T)
# axis(1, las=2,  at=c(1,2,3,4,5,6,7,8,9,10), labels=c("Control", "Fence", "K", "N", "NK", "NP", "NPK", "NPK+Fen.","P", "PK"), cex.axis=0.9)
mtext(side=3, text="a",line=-1.5, adj=0.02, cex=1.25, font=2)
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
plot(mat$mean, ylab=expression(''*italic(t[1/4])*'  (years)'), ylim=c(min(lower),max(upper)),yaxt="n", tck=-0.01, xaxt="n", xlab="", las=1, cex.lab=2, cex.axis=1.5 )
axis(side = 2, at = seq(1.3, 1.8, by = 0.05), las=2, cex.axis=1.4, tck=-0.01)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
circle<-rep(.3,8)
symbols(x=position, mat$mean, circles=circle,inches=1/14,ann=F, bg=c("black"), fg=NULL, add=T)
# axis(1, las=2,  at=c(1,2,3,4,5,6,7,8,9,10), labels=c("Control", "Fence", "K", "N", "NK", "NP", "NPK", "NPK+Fen.","P", "PK"), cex.axis=0.9)
mtext(side=3, text="b",line=-1.5, adj=0.02, cex=1.25, font=2)
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
plot(mat$mean, ylab=expression(''*italic(k[a])~' ('~years^-1*')'),xlab="Treatment", yaxt="n", ylim=c(min(lower),max(upper)), tck=-0.01, xaxt="n",las=1, cex.lab=2, tck=-0.01, cex.axis=1.5)
axis(side = 2, at = seq(0.2, 0.5, by = 0.05), las=2, cex.axis=1.4, tck=-0.01)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA, yaxt="n")
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
circle<-rep(.3,8)
symbols(x=position, mat$mean, circles=circle,inches=1/14,ann=F, bg=c("black"), fg=NULL, add=T)
axis(1, las=1,  at=c(1,2,3,4,5,6,7,8), labels=c("Control", "N", "P", "K", "NP", "NK", "PK", "NPK"), cex.axis=1.5, tck=-0.01)
mtext(side=3, text="c",line=-1.5, adj=0.02, cex=1.25, font=2)
box()
##########################

########################### Figure 3
mean<-tapply((Nutdata$single.k),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$single.k),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$single.k,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]

par(mar=c(0.5,5.5,5.5,1), mfrow=c(2,2), mgp=c(3,0.5,0))
plot(mat$mean, ylab=expression(''*italic(k[s])~' ('~years^-1*')'), ylim=c(min(lower),max(upper)), tck=-0.01,xaxt="n", xlab="", las=1, cex.lab=2, cex.axis=1.5)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
circle<-rep(.3,8)
symbols(x=position, mat$mean, circles=circle,inches=1/14,ann=F, bg=c("black"), fg=NULL, add=T)
# axis(1, las=2,  at=c(1,2,3,4,5,6,7,8,9,10), labels=c("Control", "Fence", "K", "N", "NK", "NP", "NPK", "NPK+Fen.","P", "PK"), cex.axis=0.9)
mtext(side=3, text="a",line=-1.5, adj=0.02, cex=1.25, font=2)
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
mtext(side=3, text="b",line=-1.5, adj=0.02, cex=1, font=2)
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
mtext(side=3, text="c",line=-1.5, adj=0.02, cex=1, font=2)
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
plot(mat$mean, ylab=expression('Asymptotic '*italic(A)*''), ylim=c(min(lower),max(upper)), tck=-0.01, xaxt="n", xlab="Treatment", las=1, cex.lab=2, tck=-0.01, cex.axis=1.5)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
circle<-rep(.3,8)
symbols(x=position, mat$mean, circles=circle,inches=1/14,ann=F, bg=c("black"), fg=NULL, add=T)
axis(1, las=1,  at=c(1,2,3,4,5,6,7,8), labels=c("Control", "N", "P", "K", "NP", "NK", "PK", "NPK"), cex.axis=0.85, tck=-0.01)
mtext(side=3, text="d",line=-1.5, adj=0.02, cex=1.25, font=2)
box()



################## Graphical abstract
mean<-tapply((Nutdata$weibull.tenth),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$weibull.tenth),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$weibull.tenth,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]

par(mar=c(0.5,6,5.5,1), mfcol=c(2,2), mgp=c(3,0.5,0))
plot(mat$mean, ylab=expression(''*italic(t[1/10])*'  (years)'), ylim=c(min(lower),max(upper)), yaxt="n", tck=-0.01,xaxt="n", xlab="", las=1, cex.lab=2, cex.axis=1.5, cex.main=1.75, main='Early-stage \n decomposition')
axis(side = 2, at = seq(0.4, 0.8, by = 0.05), las=2, cex.axis=1.4, tck=-0.01)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
circle<-rep(.3,8)
symbols(x=position, mat$mean, circles=circle,inches=1/14,ann=F, bg=c("black"), fg=NULL, add=T)
# axis(1, las=2,  at=c(1,2,3,4,5,6,7,8,9,10), labels=c("Control", "Fence", "K", "N", "NK", "NP", "NPK", "NPK+Fen.","P", "PK"), cex.axis=0.9)
# mtext(side=3, text="a",line=-1.5, adj=0.02, cex=1.25, font=2)
box()

mean<-tapply((Nutdata$asymp.k),Nutdata$Trt,  mean, na.rm=T)
se<-tapply((Nutdata$asymp.k),Nutdata$Trt,  sd, na.rm=T)/sqrt(tapply(Nutdata$asymp.k,Nutdata$Trt,  nnzero, na.counted=F))
upper<-mean+se
lower<-mean-se
upper
lower
mat<-as.data.frame(cbind(mean, se, upper, lower))
mat <- mat[c(1,3,7,2,5,4,8,6),]
par(mar=c(5.5,6,0.5,1))
plot(mat$mean, ylab=expression(''*italic(k[a])~' ('~years^-1*')'),xlab="Treatment", yaxt="n", ylim=c(min(lower),max(upper)), tck=-0.01, xaxt="n",las=1, cex.lab=2, tck=-0.01, cex.axis=1.5)
axis(side = 2, at = seq(0.2, 0.5, by = 0.05), las=2, cex.axis=1.4, tck=-0.01)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA, yaxt="n")
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
circle<-rep(.3,8)
symbols(x=position, mat$mean, circles=circle,inches=1/14,ann=F, bg=c("black"), fg=NULL, add=T)
axis(1, las=1,  at=c(1,2,3,4,5,6,7,8), labels=c("Control", "N", "P", "K", "NP", "NK", "PK", "NPK"), cex.axis=0.85, tck=-0.01)
# mtext(side=3, text="c",line=-1.5, adj=0.02, cex=1.25, font=2)
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
plot(mat$mean, ylab=expression(''*italic(MRT)~' ('~years*')'), ylim=c(min(lower),max(upper)), tck=-0.01,xaxt="n", xlab="", las=1, cex.lab=2, cex.axis=1.5, cex.main=1.75, main='Late-stage \n decomposition')
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
circle<-rep(.3,8)
symbols(x=position, mat$mean, circles=circle,inches=1/14,ann=F, bg=c("black"), fg=NULL, add=T)
# axis(1, las=2,  at=c(1,2,3,4,5,6,7,8,9,10), labels=c("Control", "Fence", "K", "N", "NK", "NP", "NPK", "NPK+Fen.","P", "PK"), cex.axis=0.9)
# mtext(side=3, text="b",line=-1.5, adj=0.02, cex=1, font=2)
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
plot(mat$mean, ylab=expression('Asymptotic '*italic(A)*''), ylim=c(min(lower),max(upper)), tck=-0.01, xaxt="n", xlab="Treatment", las=1, cex.lab=2, tck=-0.01, cex.axis=1.5)
rect(0, mat$lower[1], 11, mat$upper[1], col="light grey", border=NA)
position<-c(1,2,3,4,5,6,7,8)
errbar(position, mat$mean, mat$upper, mat$lower, cap=0.015, col="#4e79a7", xaxt="n", add=T, lwd=2)
circle<-rep(.3,8)
symbols(x=position, mat$mean, circles=circle,inches=1/14,ann=F, bg=c("black"), fg=NULL, add=T)
axis(1, las=1,  at=c(1,2,3,4,5,6,7,8), labels=c("Control", "N", "P", "K", "NP", "NK", "PK", "NPK"), cex.axis=0.85, tck=-0.01)
# mtext(side=3, text="d",line=-1.5, adj=0.02, cex=1.25, font=2)
box()

##########################################################
##################### Figure 4
######################################################################
mean<-tapply(Nutdata$single.k, Nutdata$N, mean, na.rm=T )
t1<-seq(0,10,length.out=10000)
par(mfrow=c(2,3), oma=c(0,0,0,0), mar=c(5,5.5,2,1), mgp=c(3,0.8,0), cex.lab=2,  cex.main=1.75, cex.axis=1.5)
MLdata<-read.csv("5.23.2019 Nutnet Masterlist.csv",  header=T, sep=",")
head(MLdata)
plot(NULL, pch=1, col="grey87",xlab="Time (years)",ylab="Prop. C Remaining", main="Single Exponential", xlim=c(0,10), ylim=c(0,1))
lines(t1,exp(-mean[1]*t1), col="grey50", lty=2, type="l", lwd=3)
lines(t1,exp(-mean[2]*t1), col="forest green", lty=1, lwd=3)
legend("topright", c("Control","+Nitrogen"), col=c("grey50","forest green"), lty=c(1,1), bty="n", cex=1.5, lwd=c(2,2))
mtext(side=1, text="a",line=-1.4, adj=0.02, cex=1.5, font=2)

mean_ka<-tapply(Nutdata$asymp.k, Nutdata$N, mean, na.rm=T )
mean_A<-tapply(Nutdata$asymp.A, Nutdata$N, mean, na.rm=T )
plot(NULL, pch=16, col="grey87", xlab="Time (years)",ylab="Prop. C Remaining", main="Asymptotic Exponential", xlim=c(0,10), ylim=c(0,1))
lines(t1,mean_A[1]+((1-mean_A[1])*exp(-mean_ka[1]*t1)), col="grey50", lwd=3)
lines(t1,mean_A[2]+((1-mean_A[2])*exp(-mean_ka[2]*t1)), col="forest green", lwd=3)
legend("topright", c("Control","+Nitrogen"), col=c("grey50","forest green"), lty=c(1,1), bty="n", cex=1.5, lwd=c(2,2))
mtext(side=1, text="b",line=-1.4, adj=0.02, cex=1.5, font=2)

mean_alpha<-tapply(Nutdata$weibull.alpha, Nutdata$N, mean, na.rm=T )
mean_beta<-tapply(Nutdata$weibull.beta, Nutdata$N, mean, na.rm=T )
weibull.fit<-exp(-(t1/mean_beta[1])^mean_alpha[1])
plot(NULL, pch=16, col="grey87",xlab="Time (years)",ylab="Prop. C Remaining", main="Weibull", xlim=c(0,10), ylim=c(0,1))
lines(t1,exp(-(t1/mean_beta[1])^mean_alpha[1]), col="grey50", lwd=3)
lines(t1,exp(-(t1/mean_beta[2])^mean_alpha[2]), col="forest green", lwd=3)
legend("topright", c("Control","+Nitrogen"), col=c("grey50","forest green"), lty=c(1,1), bty="n", cex=1.25, lwd=c(2,2))
mtext(side=1, text="c",line=-1.4, adj=0.02, cex=1.5, font=2)


##########################################################################

###### ANCOVAS Table  S7

##
model <- lme(sqrt(weibull.tenth)~pH+N+P+K+N:P+N:K+P:K+N:P:K+N:pH+P:pH+K:pH,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),
      as.data.frame(as.matrix(summary(model)$tTable[,"p-value"])))[c(2:12),]
model <- lme(sqrt(weibull.quarter)~pH+N+P+K+N:P+N:K+P:K+N:P:K+N:pH+P:pH+K:pH,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),
      as.data.frame(as.matrix(summary(model)$tTable[,"p-value"])))[c(2:12),]
model <- lme(sqrt(weibull.half.life)~pH+N+P+K+N:P+N:K+P:K+N:P:K+N:pH+P:pH+K:pH,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),
      as.data.frame(as.matrix(summary(model)$tTable[,"p-value"])))[c(2:12),]
model <- lme(log(weibull.mrt)~pH+N+P+K+N:P+N:K+P:K+N:P:K+N:pH+P:pH+K:pH,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),
      as.data.frame(as.matrix(summary(model)$tTable[,"p-value"])))[c(2:12),]
model <- lme(log(weibull.alpha)~pH+N+P+K+N:P+N:K+P:K+N:P:K+N:pH+P:pH+K:pH,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),
      as.data.frame(as.matrix(summary(model)$tTable[,"p-value"])))[c(2:12),]
model <- lme(sqrt(single.k)~pH+N+P+K+N:P+N:K+P:K+N:P:K+N:pH+P:pH+K:pH,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),
      as.data.frame(as.matrix(summary(model)$tTable[,"p-value"])))[c(2:12),]
model <- lme(log(asymp.k)~pH+N+P+K+N:P+N:K+P:K+N:P:K+N:pH+P:pH+K:pH,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),
      as.data.frame(as.matrix(summary(model)$tTable[,"p-value"])))[c(2:12),]
model <- lme((asymp.A)~pH+N+P+K+N:P+N:K+P:K+N:P:K+N:pH+P:pH+K:pH,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),
      as.data.frame(as.matrix(summary(model)$tTable[,"p-value"])))[c(2:12),]

##

##

##
model <- lme(sqrt(weibull.tenth)~MAPDist+N+P+K+N:P+N:K+P:K+N:P:K+N:MAPDist+P:MAPDist+K:MAPDist,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(weibull.quarter)~MAPDist+N+P+K+N:P+N:K+P:K+N:P:K+N:MAPDist+P:MAPDist+K:MAPDist,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(weibull.half.life)~MAPDist+N+P+K+N:P+N:K+P:K+N:P:K+N:MAPDist+P:MAPDist+K:MAPDist,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(weibull.mrt)~MAPDist+N+P+K+N:P+N:K+P:K+N:P:K+N:MAPDist+P:MAPDist+K:MAPDist,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(weibull.alpha)~MAPDist+N+P+K+N:P+N:K+P:K+N:P:K+N:MAPDist+P:MAPDist+K:MAPDist,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(single.k)~MAPDist+N+P+K+N:P+N:K+P:K+N:P:K+N:MAPDist+P:MAPDist+K:MAPDist,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(asymp.k)~MAPDist+N+P+K+N:P+N:K+P:K+N:P:K+N:MAPDist+P:MAPDist+K:MAPDist,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme((asymp.A)~MAPDist+N+P+K+N:P+N:K+P:K+N:P:K+N:MAPDist+P:MAPDist+K:MAPDist,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)

##
model <- lme(sqrt(weibull.tenth)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(pct_C)+P:log(pct_C)+K:log(pct_C),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(weibull.quarter)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(pct_C)+P:log(pct_C)+K:log(pct_C),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(weibull.half.life)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(pct_C)+P:log(pct_C)+K:log(pct_C),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(weibull.mrt)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(pct_C)+P:log(pct_C)+K:log(pct_C),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(weibull.alpha)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(pct_C)+P:log(pct_C)+K:log(pct_C),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(single.k)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(pct_C)+P:log(pct_C)+K:log(pct_C),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(asymp.k)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(pct_C)+P:log(pct_C)+K:log(pct_C),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme((asymp.A)~log(pct_C)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(pct_C)+P:log(pct_C)+K:log(pct_C),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)



##
model <- lme(sqrt(weibull.tenth)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(N_Dep)+P:log(N_Dep)+K:log(N_Dep),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:12),]
anova(model)
model <- lme(sqrt(weibull.quarter)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(N_Dep)+P:log(N_Dep)+K:log(N_Dep),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:12),]
anova(model)
model <- lme(sqrt(weibull.half.life)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(N_Dep)+P:log(N_Dep)+K:log(N_Dep),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:12),]
anova(model)
model <- lme(log(weibull.mrt)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(N_Dep)+P:log(N_Dep)+K:log(N_Dep),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:12),]
anova(model)
model <- lme(log(weibull.alpha)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(N_Dep)+P:log(N_Dep)+K:log(N_Dep),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:12),]
anova(model)
model <- lme(sqrt(single.k)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(N_Dep)+P:log(N_Dep)+K:log(N_Dep),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:12),]
anova(model)
model <- lme(log(asymp.k)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(N_Dep)+P:log(N_Dep)+K:log(N_Dep),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:12),]
anova(model)
model <- lme((asymp.A)~log(N_Dep)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(N_Dep)+P:log(N_Dep)+K:log(N_Dep),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:12),]
anova(model)

##
model <- lme(sqrt(weibull.tenth)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(ppm_P)+P:log(ppm_P)+K:log(ppm_P),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(weibull.quarter)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(ppm_P)+P:log(ppm_P)+K:log(ppm_P),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(weibull.half.life)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(ppm_P)+P:log(ppm_P)+K:log(ppm_P),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(weibull.mrt)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(ppm_P)+P:log(ppm_P)+K:log(ppm_P),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(weibull.alpha)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(ppm_P)+P:log(ppm_P)+K:log(ppm_P),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(single.k)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(ppm_P)+P:log(ppm_P)+K:log(ppm_P),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(asymp.k)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(ppm_P)+P:log(ppm_P)+K:log(ppm_P),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme((asymp.A)~log(ppm_P)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(ppm_P)+P:log(ppm_P)+K:log(ppm_P),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)

##
model <- lme(sqrt(weibull.tenth)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(ppm_Mn)+P:log(ppm_Mn)+K:log(ppm_Mn),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:17),]
anova(model)
model <- lme(sqrt(weibull.quarter)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(ppm_Mn)+P:log(ppm_Mn)+K:log(ppm_Mn),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:17),]
anova(model)
model <- lme(sqrt(weibull.half.life)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(ppm_Mn)+P:log(ppm_Mn)+K:log(ppm_Mn),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:17),]
anova(model)
model <- lme(log(weibull.mrt)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(ppm_Mn)+P:log(ppm_Mn)+K:log(ppm_Mn),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:17),]
anova(model)
model <- lme(log(weibull.alpha)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(ppm_Mn)+P:log(ppm_Mn)+K:log(ppm_Mn),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:17),]
anova(model)
model <- lme(sqrt(single.k)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(ppm_Mn)+P:log(ppm_Mn)+K:log(ppm_Mn),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:17),]
anova(model)
model <- lme(log(asymp.k)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(ppm_Mn)+P:log(ppm_Mn)+K:log(ppm_Mn),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:17),]
anova(model)
summary(model)
model <- lme((asymp.A)~log(ppm_Mn)+N+P+K+N:P+N:K+P:K+N:P:K+N:log(ppm_Mn)+P:log(ppm_Mn)+K:log(ppm_Mn),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:17),]
anova(model)


##
model <- lme(sqrt(weibull.tenth)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K+N:MAT_v2+P:MAT_v2+K:MAT_v2,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(weibull.quarter)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K+N:MAT_v2+P:MAT_v2+K:MAT_v2,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(weibull.half.life)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K+N:MAT_v2+P:MAT_v2+K:MAT_v2,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(weibull.mrt)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K+N:MAT_v2+P:MAT_v2+K:MAT_v2,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(weibull.alpha)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K+N:MAT_v2+P:MAT_v2+K:MAT_v2,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(single.k)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K+N:MAT_v2+P:MAT_v2+K:MAT_v2,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(asymp.k)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K+N:MAT_v2+P:MAT_v2+K:MAT_v2,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme((asymp.A)~MAT_v2+N+P+K+N:P+N:K+P:K+N:P:K+N:MAT_v2+P:MAT_v2+K:MAT_v2,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)

##
model <- lme(sqrt(weibull.tenth)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K+N:MAP_v2+P:MAP_v2+K:MAP_v2,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(weibull.quarter)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K+N:MAP_v2+P:MAP_v2+K:MAP_v2,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(weibull.half.life)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K+N:MAP_v2+P:MAP_v2+K:MAP_v2,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(weibull.mrt)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K+N:MAP_v2+P:MAP_v2+K:MAP_v2,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(weibull.alpha)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K+N:MAP_v2+P:MAP_v2+K:MAP_v2,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(single.k)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K+N:MAP_v2+P:MAP_v2+K:MAP_v2,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(asymp.k)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K+N:MAP_v2+P:MAP_v2+K:MAP_v2,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme((asymp.A)~MAP_v2+N+P+K+N:P+N:K+P:K+N:P:K+N:MAP_v2+P:MAP_v2+K:MAP_v2,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)

##

##
model <- lme(sqrt(weibull.tenth)~sqrt(Total_mass_mean)+N+P+K+N:P+N:K+P:K+N:P:K+N:sqrt(Total_mass_mean)+P:sqrt(Total_mass_mean)+K:sqrt(Total_mass_mean),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(weibull.quarter)~sqrt(Total_mass_mean)+N+P+K+N:P+N:K+P:K+N:P:K+N:sqrt(Total_mass_mean)+P:sqrt(Total_mass_mean)+K:sqrt(Total_mass_mean),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(weibull.half.life)~sqrt(Total_mass_mean)+N+P+K+N:P+N:K+P:K+N:P:K+N:sqrt(Total_mass_mean)+P:sqrt(Total_mass_mean)+K:sqrt(Total_mass_mean),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(weibull.mrt)~sqrt(Total_mass_mean)+N+P+K+N:P+N:K+P:K+N:P:K+N:sqrt(Total_mass_mean)+P:sqrt(Total_mass_mean)+K:sqrt(Total_mass_mean),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(weibull.alpha)~sqrt(Total_mass_mean)+N+P+K+N:P+N:K+P:K+N:P:K+N:sqrt(Total_mass_mean)+P:sqrt(Total_mass_mean)+K:sqrt(Total_mass_mean),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(sqrt(single.k)~sqrt(Total_mass_mean)+N+P+K+N:P+N:K+P:K+N:P:K+N:sqrt(Total_mass_mean)+P:sqrt(Total_mass_mean)+K:sqrt(Total_mass_mean),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme(log(asymp.k)~sqrt(Total_mass_mean)+N+P+K+N:P+N:K+P:K+N:P:K+N:sqrt(Total_mass_mean)+P:sqrt(Total_mass_mean)+K:sqrt(Total_mass_mean),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)
model <- lme((asymp.A)~sqrt(Total_mass_mean)+N+P+K+N:P+N:K+P:K+N:P:K+N:sqrt(Total_mass_mean)+P:sqrt(Total_mass_mean)+K:sqrt(Total_mass_mean),random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
cbind(as.data.frame(as.matrix(as.data.frame(summary(model)$tTable[,"Value"]))),as.data.frame(as.matrix(anova(model)[,"F-value"])),
      as.data.frame(as.matrix(anova(model)[,"p-value"])))[c(2:9),]
anova(model)



##########################################################
# Covariate interactions 

My_Theme = theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 14))

Nutdata$Sqrt_single.k<-sqrt(Nutdata$single.k)
Nutdata$Sqrt_weibull.half.life<-sqrt(Nutdata$weibull.half.life)
Nutdata$Sqrt_single.k<-sqrt(Nutdata$single.k)
Nutdata$Sqrt_weibull.tenth.life<-sqrt(Nutdata$weibull.tenth)
Nutdata$Sqrt_weibull.quarter.life<-sqrt(Nutdata$weibull.quarter)
Nutdata$Log_weibull.alpha<-log(Nutdata$weibull.alpha)
Nutdata$Log_weibull.mrt<-log(Nutdata$weibull.mrt)
Nutdata$Sqrt_Total_mass_mean<-sqrt(Nutdata$Total_mass_mean)
Nutdata$Log_NDep<-log(Nutdata$N_Dep_Ackerman)
Nutdata$Log_ppm_Mn<-log(Nutdata$ppm_Mn)



####################### Mn Figure S1
model <- lme(Sqrt_single.k~Log_ppm_Mn+N+P+K+N:P+N:K+P:K+N:P:K+N:Log_ppm_Mn+P:Log_ppm_Mn+K:Log_ppm_Mn,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("Log_ppm_Mn",  "N"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P +  xlab(bquote('Log (Soil Mn) (ppm)')) +ylab(expression('Sqrt ('*italic(k[s])~') ('~years^-1*')'))+
  theme(legend.position = c(0.15, 0.8), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid")) + guides(col=guide_legend("Nitrogen"))


Mn_R1<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                 plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),) + geom_text(
                   label="a", 
                   x=1.2,
                   y=.3,
                   color = "black", size=6)
Mn_R1

model <- lme(Log_asymp.k~Log_ppm_Mn+N+P+K+N:P+N:K+P:K+N:P:K+N:Log_ppm_Mn+P:Log_ppm_Mn+K:Log_ppm_Mn,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("Log_ppm_Mn",  "N"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P +  xlab(bquote('Log (Soil Mn) (ppm)')) +ylab(expression('Log ('*italic(k[a])~') ('~years^-1*')'))+
  theme(legend.position = c(0.15, 0.8), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))+ guides(col=guide_legend("Nitrogen"))
Mn_R2<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                 plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                   label="b", 
                   x=1.2,
                   y=-2.2,
                   color = "black", size=6)
Mn_R2

SixPanel<-ggarrange(Mn_R1, Mn_R2,
                    NA,NA,
                    ncol = 2, nrow = 2)

ggsave(filename = "20210803_FigS1.pdf",
       plot = SixPanel,
       dpi=300,
       width = 18, height = 22, units = "cm")

##########################################################


#### Soil pH
model <- lme(Sqrt_weibull.quarter~pH+N+P+K+N:P+N:K+P:K+N:P:K+N:pH+P:pH+K:pH,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("pH", "K"),title = "", group.terms = c(1, 2), 
              colors = c("black", "red")) + My_Theme 
Q<-P +  xlab(bquote('soil pH')) +ylab(expression('Sqrt (Weibull '*italic(t[1/4])*') (years)'))+
  theme(legend.position = c(0.15, 0.8), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))+ guides(col=guide_legend("Potassium"))
pH_R1<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                 plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                   label="a", 
                   x=4.25,
                   y=0.2,
                   color = "black", size=6)
pH_R1


model <- lme(Sqrt_weibull.half.life~pH+N+P+K+N:P+N:K+P:K+N:P:K+N:pH+P:pH+K:pH,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("pH", "P"),title = "", group.terms = c(1, 2), 
              colors = c("black", "purple")) + My_Theme 
Q<-P +  xlab(bquote('soil pH')) +ylab(expression('Sqrt (Weibull '*italic(t[1/2])*') (years)'))+
  theme(legend.position = c(0.15, 0.8), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))+ guides(col=guide_legend("Phosphorus"))
pH_R2<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                 plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                   label="b", 
                   x=4.25,
                   y=0.5,
                   color = "black", size=6)
pH_R2

model <- lme(Log_asymp.k~pH+N+P+K+N:P+N:K+P:K+N:P:K+N:pH+P:pH+K:pH,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("pH", "N"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P +  xlab(bquote('soil pH')) +ylab(expression('Log ('*italic(k[a])~') ('~years^-1*')'))+
  theme(legend.position = c(0.15, 0.8), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))+ guides(col=guide_legend("Nitrogen"))
pH_R3<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                 plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                   label="c", 
                   x=4.25,
                   y=-2.2,
                   color = "black", size=6)
pH_R3

SixPanel<-ggarrange(pH_R1, pH_R2,
                    pH_R3,NA,
                    ncol = 2, nrow = 2)

ggsave(filename = "20210803_FigS2.pdf",
       plot = SixPanel,
       dpi=300,
       width = 18, height = 22, units = "cm")

##########################################################


########### N Dep

model <- lme(Sqrt_weibull.half.life~Log_NDep+N+P+K+N:P+N:K+P:K+N:P:K+N:Log_NDep+P:Log_NDep+K:Log_NDep,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("Log_NDep",  "N"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P +  xlab(bquote('Log (N Dep.)'*~'(kg N '~km^-1 ~year^-1*')')) +ylab(expression('Sqrt (Weibull '*italic(t[1/2])*') (years)'))+
  theme(legend.position = c(0.15, 0.8), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))+ guides(col=guide_legend("Nitrogen"))
NDep_R3<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                   plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                     label="a", 
                     x=.4,
                     y=0.5,
                     color = "black", size=6)
NDep_R3

model <- lme(Log_weibull.alpha~Log_NDep+N+P+K+N:P+N:K+P:K+N:P:K+N:Log_NDep+P:Log_NDep+K:Log_NDep,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("Log_NDep",  "N"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P +  xlab(bquote('Log (N Dep.)'*~'(kg N '~km^-1 ~year^-1*')')) +ylab(expression('Log (Weibull '*italic(alpha)*')'))+
  theme(legend.position =  c(0.15, 0.8), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))+ guides(col=guide_legend("Nitrogen"))
NDep_R4<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                   plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                     label="b", 
                     x=.4,
                     y=-1.75,
                     color = "black", size=6)
NDep_R4

model <- lme(Sqrt_single.k~Log_NDep+N+P+K+N:P+N:K+P:K+N:P:K+N:Log_NDep+P:Log_NDep+K:Log_NDep,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("Log_NDep",  "N"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P +  xlab(bquote('Log (N Dep.)'*~'(kg N '~km^-1 ~year^-1*')')) +ylab(expression('Sqrt ('*italic(k[s])~') ('~years^-1*')'))+
  theme(legend.position = c(0.15, 0.8), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))+ guides(col=guide_legend("Nitrogen"))
NDep_R6<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                   plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                     label="c", 
                     x=.4,
                     y=0.25,
                     color = "black", size=6)
NDep_R6

model <- lme(asymp.A~Log_NDep+N+P+K+N:P+N:K+P:K+N:P:K+N:Log_NDep+P:Log_NDep+K:Log_NDep,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("Log_NDep",  "N"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P +  xlab(bquote('Log (N Dep.)'*~'(kg N '~km^-1 ~year^-1*')')) +ylab(expression('Asymptotic '*italic(A)~' '))+
  theme(legend.position = c(0.15, 0.8), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))+ guides(col=guide_legend("Nitrogen"))
NDep_R7<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                   plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                     label="d", 
                     x=.4,
                     y=-.1,
                     color = "black", size=6)
NDep_R7

SixPanel<-ggarrange(NDep_R3,NDep_R4,
                    NDep_R6,NDep_R7,
                    ncol = 2, nrow = 2)
ggsave(filename = "20210803_FigS3.pdf",
       plot = SixPanel,
       dpi=300,
       width = 18, height = 22, units = "cm")


#### MAPDist


model <- lme(Sqrt_weibull.quarter.life~MAPDist+N+P+K+N:P+N:K+P:K+N:P:K+N:MAPDist+P:MAPDist+K:MAPDist,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("MAPDist", "K"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P +  xlab(bquote('Precipitation Distribution')) +ylab(expression('Sqrt (Weibull '*italic(t[1/2])*') (years)'))+
  theme(legend.position = c(0.1, 0.9), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))
MAPDist_R1<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                      plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                        label="a", 
                        x=.075,
                        y=0.75,
                        color = "black", size=6)
MAPDist_R1

model <- lme(Sqrt_weibull.half.life~MAPDist+N+P+K+N:P+N:K+P:K+N:P:K+N:MAPDist+P:MAPDist+K:MAPDist,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("MAPDist", "N"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P +  xlab(bquote('Precipitation Distribution')) +ylab(expression('Sqrt (Weibull '*italic(t[1/2])*') (years)'))+
  theme(legend.position = c(0.1, 0.9), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))
MAPDist_R1<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                      plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                        label="a", 
                        x=.075,
                        y=0.75,
                        color = "black", size=6)
MAPDist_R1

model <- lme(Log_weibull.alpha~MAPDist+N+P+K+N:P+N:K+P:K+N:P:K+N:MAPDist+P:MAPDist+K:MAPDist,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("MAPDist", "N"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P +  xlab(bquote('Precipitation Distribution')) +ylab(bquote('Log (Weibull '*italic(alpha)~')'))+
  theme(legend.position = c(0.1, 0.9), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))
MAPDist_R2<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                      plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                        label="b", 
                        x=.075,
                        y=-1.5,
                        color = "black", size=6)
MAPDist_R2

model <- lme(Log_weibull.mrt~MAPDist+N+P+K+N:P+N:K+P:K+N:P:K+N:MAPDist+P:MAPDist+K:MAPDist,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("MAPDist", "N"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P +  xlab(bquote('Precipitation Distribution')) +ylab(bquote('Log (Weibull '*italic(MRT)~') (years)'))+
  theme(legend.position = c(0.1, 0.9), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))
MAPDist_R3<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                      plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                        label="c", 
                        x=.075,
                        y=0,
                        color = "black", size=6)
MAPDist_R3

model <- lme(Sqrt_single.k~MAPDist+N+P+K+N:P+N:K+P:K+N:P:K+N:MAPDist+P:MAPDist+K:MAPDist,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("MAPDist", "N"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P +  xlab(bquote('Precipitation Distribution')) +ylab(expression('Sqrt ('*italic(k[s])~') ('~years^-1*')'))+
  theme(legend.position = c(0.1, 0.9), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))
MAPDist_R4<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                      plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                        label="d", 
                        x=.075,
                        y=0.3,
                        color = "black", size=6)
MAPDist_R4



model <- lme(asymp.A~MAPDist+N+P+K+N:P+N:K+P:K+N:P:K+N:MAPDist+P:MAPDist+K:MAPDist,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("MAPDist", "N"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P +  xlab(bquote('Precipitation Distribution')) + ylab(expression('Asymptotic '*italic(A)*''))+
  theme(legend.position = c(0.1, 0.9), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))
MAPDist_R5<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                      plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                        label="e", 
                        x=.075,
                        y=-.15,
                        color = "black", size=6)
MAPDist_R5

SixPanel<-ggarrange(MAPDist_R1, MAPDist_R2, MAPDist_R3,
                    MAPDist_R4,MAPDist_R5,NA,
                    ncol = 3, nrow = 2)

ggsave(filename = "20210803_FigS4.pdf",
       plot = SixPanel,
       dpi=300,
       width = 22, height = 22, units = "cm")

#### Total Mass
model <- lme(Sqrt_single.k~Sqrt_Total_mass_mean+N+P+K+N:P+N:K+P:K+N:P:K+N:Sqrt_Total_mass_mean+P:Sqrt_Total_mass_mean+K:Sqrt_Total_mass_mean,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)
P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("Sqrt_Total_mass_mean", "N"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P +  xlab(bquote('Aboveground Biomass (g '~m^-2*')')) +ylab(expression('Sqrt ('*italic(k[s])~') ('~years^-1*')'))+
  theme(legend.position = c(0.1, 0.9), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))
TotalMass_R1<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                        plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                          label="a", 
                          x=6,
                          y=.3,
                          color = "black", size=6)
TotalMass_R1


model <- lme(Log_weibull.alpha~Sqrt_Total_mass_mean+N+P+K+N:P+N:K+P:K+N:P:K+N:Sqrt_Total_mass_mean+P:Sqrt_Total_mass_mean+K:Sqrt_Total_mass_mean,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)

P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("Sqrt_Total_mass_mean", "N"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P + xlab(bquote('Aboveground Biomass (g '~m^-2*')')) +ylab(bquote('Log (Weibull '*italic(alpha)~')'))+
  theme(legend.position = c(0.1, 0.9), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))
TotalMass_R2<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                        plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                          label="b", 
                          x=6,
                          y=-1.5,
                          color = "black", size=6)
TotalMass_R2

model <- lme(Log_weibull.mrt~Sqrt_Total_mass_mean+N+P+K+N:P+N:K+P:K+N:P:K+N:Sqrt_Total_mass_mean+P:Sqrt_Total_mass_mean+K:Sqrt_Total_mass_mean,random=~1|Site_Text/Block,data=Nutdata, na.action=na.omit)
summary(model)
anova(model)

P<-plot_model(model, type = "pred",  show.data=TRUE, ci.lvl = 0.95,terms = c("Sqrt_Total_mass_mean", "N"),title = "", group.terms = c(1, 2), 
              colors = c("black", "forest green")) + My_Theme 
Q<-P + xlab(bquote('Aboveground Biomass (g '~m^-2*')')) +ylab(bquote('Log (Weibull '*italic(MRT)~') (years)'))+
  theme(legend.position = c(0.1, 0.9), legend.background = element_rect(color = F,fill = F, size = 2, linetype = "solid"))
TotalMass_R3<-Q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_rect(colour = "black", fill=NA, size=1), axis.line = element_line(colour = "black"),
                        plot.margin = margin(1, 0.5, 0.5, 0.5, "cm"),)+ geom_text(
                          label="c", 
                          x=6,
                          y=0.5,
                          color = "black", size=6)
TotalMass_R3

SixPanel<-ggarrange(TotalMass_R1, TotalMass_R2,
                    TotalMass_R3,NA,
                    ncol = 2, nrow = 2)

ggsave(filename = "20210803_FigS5.pdf",
       plot = SixPanel,
       dpi=300,
       width = 22, height = 22, units = "cm")


