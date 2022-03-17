#Template to fit single, double, asymptotic, and weibull litter decomposition models to indiviudal plots

#Open necessary libraries
library(readr)
library(stats4)
library(bbmle)
library(mltools)
library(Metrics)

#Download data file. Users will need to update link to reflect their own file path. 
urlfile = "https://raw.githubusercontent.com/gill20a/NutNetLitterDecomp/master/ModelFitting/NutNet_OakLitterDecay_HarvestData.csv"
data<-read.csv(url(urlfile))

#Calculate litter CN ratio
data$CN_Ratio<-data$Prop_Init_C_Mass/data$Fin_N_Mass_g
hist(data$CN_Ratio, breaks=150, main="Litter CN Ratio")
mean(data$CN_Ratio, na.rm=T)-(2*sd(data$CN_Ratio, na.r=T))
abline(v=mean(data$CN_Ratio, na.rm=T), lwd=3, lty=2, col="red")
abline(v=(mean(data$CN_Ratio, na.rm=T)+(2*sd(data$CN_Ratio, na.rm=T))), col="dark green")
mean(data$CN_Ratio, na.rm=T)+(2*sd(data$CN_Ratio, na.rm=T))
quantile(data$CN_Ratio, 0.975, na.rm=T)


#Calculate litter CN ratio as SD
par(mfrow=c(2,2))
pick1<-data$Harvest_No==1
hist(data[pick1, "CN_Ratio"], main="Harvest 1", xlab=c("Litter CN"))
abline(v=mean(data[pick1,"CN_Ratio"], na.rm=T), lwd=3, lty=2, col="red")
abline(v=(mean(data[pick1,"CN_Ratio"], na.rm=T)+(2*sd(data[pick1,"CN_Ratio"], na.rm=T))), col="dark green", lwd=3, lty=3)
abline(v=(mean(data[pick1,"CN_Ratio"], na.rm=T)-(2*sd(data[pick1,"CN_Ratio"], na.rm=T))), col="dark green",lwd=3, lty=3)

#Calculate litter CN ratio as 95% confidence interval
####### Create plot showing liter CN mean and confidnece intervals
par(mfcol=c(3,3), mgp=c(2,1,0), mar=c(4,4,2,1))
pick1<-data$Harvest_No==1
hist(data[pick1, "CN_Ratio"], main="Harvest 1", xlab=c("Litter CN"), border=NULL)
abline(v=mean(data[pick1,"CN_Ratio"], na.rm=T), lwd=3, lty=2, col="steelblue1")
abline(v=(quantile(data[pick1,"CN_Ratio"], 0.975,na.rm=T)), col="seagreen3", lwd=3, lty=3)
abline(v=(quantile(data[pick1,"CN_Ratio"], 0.025,na.rm=T)), col="seagreen3", lwd=3, lty=3)
box()
legend("topright", "(a)", bty="n")

pick1<-data$Harvest_No==2
hist(data[pick1, "CN_Ratio"], main="Harvest 2", xlab=c("Litter CN"), border=NULL)
abline(v=mean(data[pick1,"CN_Ratio"], na.rm=T), lwd=3, lty=2, col="steelblue1")
abline(v=(quantile(data[pick1,"CN_Ratio"], 0.975,na.rm=T)), col="seagreen3", lwd=3, lty=3)
abline(v=(quantile(data[pick1,"CN_Ratio"], 0.025,na.rm=T)), col="seagreen3", lwd=3, lty=3)
box()
legend("topright", "(b)", bty="n")

pick1<-data$Harvest_No==3
hist(data[pick1, "CN_Ratio"], main="Harvest 3", xlab=c("Litter CN"), border=NULL)
abline(v=mean(data[pick1,"CN_Ratio"], na.rm=T), lwd=3, lty=2, col="steelblue1")
abline(v=(quantile(data[pick1,"CN_Ratio"], 0.975,na.rm=T)), col="seagreen3", lwd=3, lty=3)
abline(v=(quantile(data[pick1,"CN_Ratio"], 0.025,na.rm=T)), col="seagreen3", lwd=3, lty=3)
box()
legend("topright", "(c)", bty="n")

pick1<-data$Harvest_No==4
hist(data[pick1, "CN_Ratio"], main="Harvest 4", xlab=c("Litter CN"), border=NULL)
abline(v=mean(data[pick1,"CN_Ratio"], na.rm=T), lwd=3, lty=2, col="steelblue1")
abline(v=(quantile(data[pick1,"CN_Ratio"], 0.975,na.rm=T)), col="seagreen3", lwd=3, lty=3)
abline(v=(quantile(data[pick1,"CN_Ratio"], 0.025,na.rm=T)), col="seagreen3", lwd=3, lty=3)
box()
legend("topright", "(d)", bty="n")

pick1<-data$Harvest_No==5
hist(data[pick1, "CN_Ratio"], main="Harvest 5", xlab=c("Litter CN"), border=NULL)
abline(v=mean(data[pick1,"CN_Ratio"], na.rm=T), lwd=3, lty=2, col="steelblue1")
abline(v=(quantile(data[pick1,"CN_Ratio"], 0.975,na.rm=T)), col="seagreen3", lwd=3, lty=3)
abline(v=(quantile(data[pick1,"CN_Ratio"], 0.025,na.rm=T)), col="seagreen3", lwd=3, lty=3)
box()
legend("topright", "(e)", bty="n")

pick1<-data$Harvest_No==6
hist(data[pick1, "CN_Ratio"], main="Harvest 6", xlab=c("Litter CN"), border=NULL)
abline(v=mean(data[pick1,"CN_Ratio"], na.rm=T), lwd=3, lty=2, col="steelblue1")
abline(v=(quantile(data[pick1,"CN_Ratio"], 0.975,na.rm=T)), col="seagreen3", lwd=3, lty=3)
abline(v=(quantile(data[pick1,"CN_Ratio"], 0.025,na.rm=T)), col="seagreen3", lwd=3, lty=3)
box()
legend("topright", "(f)", bty="n")

pick1<-data$Harvest_No==7
hist(data[pick1, "CN_Ratio"], main="Harvest 7", xlab=c("Litter CN"), border=NULL)
abline(v=mean(data[pick1,"CN_Ratio"], na.rm=T), lwd=3, lty=2, col="steelblue1")
abline(v=(quantile(data[pick1,"CN_Ratio"], 0.975,na.rm=T)), col="seagreen3", lwd=3, lty=3)
abline(v=(quantile(data[pick1,"CN_Ratio"], 0.025,na.rm=T)), col="seagreen3", lwd=3, lty=3)
box()
################################################
### Replace litter N values that fall outside 95% CI CN ratio with NA by harvest
Harvest1<-subset(data,data[, "Harvest_No"] == 1)
Harvest1[, "Fin_N_Mass_g"][Harvest1[, "CN_Ratio"] > (quantile(Harvest1[,"CN_Ratio"], 0.975,na.rm=T))] <- NA
Harvest1[, "Fin_N_Mass_g"][Harvest1[, "CN_Ratio"] < (quantile(Harvest1[,"CN_Ratio"], 0.025,na.rm=T))] <- NA
Harvest2<-subset(data,data[, "Harvest_No"] == 2)
Harvest2[, "Fin_N_Mass_g"][Harvest2[, "CN_Ratio"] > (quantile(Harvest2[,"CN_Ratio"], 0.975,na.rm=T))] <- NA
Harvest2[, "Fin_N_Mass_g"][Harvest2[, "CN_Ratio"] < (quantile(Harvest2[,"CN_Ratio"], 0.025,na.rm=T))] <- NA
Harvest3<-subset(data,data[, "Harvest_No"] == 3)
Harvest3[, "Fin_N_Mass_g"][Harvest3[, "CN_Ratio"] > (quantile(Harvest3[,"CN_Ratio"], 0.975,na.rm=T))] <- NA
Harvest3[, "Fin_N_Mass_g"][Harvest3[, "CN_Ratio"] < (quantile(Harvest3[,"CN_Ratio"], 0.025,na.rm=T))] <- NA
Harvest4<-subset(data,data[, "Harvest_No"] == 4)
Harvest4[, "Fin_N_Mass_g"][Harvest4[, "CN_Ratio"] > (quantile(Harvest4[,"CN_Ratio"], 0.975,na.rm=T))] <- NA
Harvest4[, "Fin_N_Mass_g"][Harvest4[, "CN_Ratio"] < (quantile(Harvest4[,"CN_Ratio"], 0.025,na.rm=T))] <- NA
Harvest5<-subset(data,data[, "Harvest_No"] == 5)
Harvest5[, "Fin_N_Mass_g"][Harvest5[, "CN_Ratio"] > (quantile(Harvest5[,"CN_Ratio"], 0.975,na.rm=T))] <- NA
Harvest5[, "Fin_N_Mass_g"][Harvest5[, "CN_Ratio"] < (quantile(Harvest5[,"CN_Ratio"], 0.025,na.rm=T))] <- NA
Harvest6<-subset(data,data[, "Harvest_No"] == 6)
Harvest6[, "Fin_N_Mass_g"][Harvest6[, "CN_Ratio"] > (quantile(Harvest6[,"CN_Ratio"], 0.975,na.rm=T))] <- NA
Harvest6[, "Fin_N_Mass_g"][Harvest6[, "CN_Ratio"] < (quantile(Harvest6[,"CN_Ratio"], 0.025,na.rm=T))] <- NA
Harvest7<-subset(data,data[, "Harvest_No"] == 7)
Harvest7[, "Fin_N_Mass_g"][Harvest7[, "CN_Ratio"] > (quantile(Harvest7[,"CN_Ratio"], 0.975,na.rm=T))] <- NA
Harvest7[, "Fin_N_Mass_g"][Harvest7[, "CN_Ratio"] < (quantile(Harvest7[,"CN_Ratio"], 0.025,na.rm=T))] <- NA

data<-rbind(Harvest1, Harvest2, Harvest3,Harvest4, Harvest5, Harvest6, Harvest7)

data$Site_Text<-as.character(data$Site_Text)

CBGB<-data[data$Site_Text=="CBGB-Iowa",] #1-10
Sierra<-data[data$Site_Text=="Sierra Foothill REC ",] #11-20
Bunch<-data[data$Site_Text=="Bunchgrass - Andrews LTER ",] #21-30
Boulder<-data[data$Site_Text=="Boulder ",] #31-40
UNC<-data[data$Site_Text=="UNC ",] #41-50
CDR<-data[data$Site_Text=="Cedar Creek LTER ",] #51-60
Elliott<-data[data$Site_Text=="Elliott Chaparral Reserve ",] #61-70
Halls<-data[data$Site_Text=="Halls Prairie ",] #71-80
Lookout<-data[data$Site_Text=="Lookout - Andrews LTER ",] #81-90
Hopland<-data[data$Site_Text=="Hopland REC ",] #91-100
McL<-data[data$Site_Text=="Mclaughlin UCNRS ",] #101-110
Sagehen<-data[data$Site_Text=="Sagehen UCNRS ",] #111-120
Sedgwick<-data[data$Site_Text=="Sedgwick Reserve UCNRS ",] #121-130
Sheep<-data[data$Site_Text=="Sheep Experiment Station ",] #131-140
Spindle<-data[data$Site_Text=="Spindletop Agricultural Farm ",]#141-150
Val<-data[data$Site_Text=="Val Mustair ",]#151-160
Bogong<-data[data$Site_Text=="Bogong High Plains ",]#161-170
Burrawan<-data[data$Site_Text=="Burrawan ",]#171-180
Cowichan<-data[data$Site_Text=="Cowichan ",]#181-190
Fru<-data[data$Site_Text=="Früebüel ",]#191-200
Kiny<-data[data$Site_Text=="Kinypanial ",]#201-210
Caroline<-data[data$Site_Text=="Mt. Caroline ",]#211-220


#vectors for k values (col 1 or 1&2) and sse (col 2)
single.k=array(0,dim=c(1000,2))
double.k=array(0,dim=c(1000,3))
asymptotic.k=array(0,dim=c(1000,3))
weibull=array(0,dim=c(1000,8))

#arrays for obs,pred,res
single.pred=array(0,dim=c(1000,4))
double.pred=array(0,dim=c(1000,4))		
asymptotic.pred=array(0,dim=c(1000,6))	
weibull.pred=array(0,dim=c(1000,7))	

#arrays to store AICc values
a.aicc=array(0,dim=c(1000,10)) #3aicc; 3weights; single, asym, double

#Array to store N model output
Nmodel=array(0,dim=c(1000,2))

#array to store RMSE output
rmse=array(0,dim=c(1000,4)) 

length.array = array(0,dim=c(1000,1)) 
N.length.array = array(0,dim=c(1000,1)) 


#Here, models are fit one at a time. You could make this a loop if you 
#are convinced all your models will converge using the same starting values.
#I never have such luck.

#choose rep or set of CBGB points to fit(here identified by the CBGB column ID)

half.life.calc=function(nls.mod){
  pars= coef(nls.mod)
  hl=pars[1] * (log(2))^(1/pars[2])
  names(hl)="half.life"
  return(hl)
}

quarter.life.calc=function(nls.mod){
  pars= coef(nls.mod)
  hl=pars[1] * (log(1/(3/4)))^(1/pars[2])
  names(hl)="half.life"
  return(hl)
}

tenth.life.calc=function(nls.mod){
  pars= coef(nls.mod)
  hl=pars[1] * (log(1/(9/10)))^(1/pars[2])
  names(hl)="half.life"
  return(hl)
}

mrt.calc=function(nls.mod){
  pars= coef(nls.mod)
  mrt=pars[1] * gamma(1+(1/pars[2]))
  names(mrt)="mrt"
  return(mrt)
}

LLS = function(y,k){
  Mhat=1*exp(-k*xNA$t) #creates a vector of=length to obs of preds
  ssq = sum((Mhat - xNA$Mt)^2)
  sigma = sqrt(ssq/length(xNA$Mt))
  return(-sum(dnorm(xNA$Mt,mean=Mhat,sd=sigma,log = TRUE)))
}

LLA = function(y,k,A){
  Mhat=A+((1-A)*exp(-k*xNA$t)) #creates a vector of=length to obs of preds
  ssq = sum((Mhat - xNA$Mt)^2)
  sigma = sqrt(ssq/length(xNA$Mt))
  return(-sum(dnorm(xNA$Mt,mean=Mhat,sd=sigma,log = TRUE)))
}

LLD = function(y,ks,k,A){
  Mhat=A*exp(-ks*xNA$t)+(1-A)*exp(-k*xNA$t) #creates a vector of=length to obs of preds
  ssq = sum((Mhat - xNA$Mt)^2)
  #browser()
  sigma = sqrt(ssq/length(xNA$Mt))
  return(-sum(dnorm(xNA$Mt,mean=Mhat,sd=sigma,log = TRUE)))
}


#Here, models are fit one at a time. You could make this a loop if you 
#are convinced all your models will converge using the same starting values.
#I never have such luck.

#choose rep or set of CBGB points to fit (here identified by the CBGB column ID)

##CBGB######################################################################################################
unique(CBGB$Plot)
i<-1
s1=subset(CBGB,CBGB$Plot==1)
t=(as.numeric(s1$Yrs_Since_Deployment))
t1<-seq(0,7,length.out=10000)
Mt=s1$Prop_Init_C_Mass
x <- data.frame(t, Mt)
x[3,] <- NA
xNA<- na.exclude(x)
Mt<-xNA$Mt
t<-xNA$t
length.array[i,1] = length(t)

#Weibull function
fit<- nls((Mt) ~ exp(-(t/beta)^alpha), start =list(beta=1, alpha = .81), algorithm="port", lower=c(0.0001,0.0001))
weibull[i,1]=coef(fit)[2]
weibull[i,2]=coef(fit)[1]
weibull[i,3]<-mrt.calc(fit)
weibull[i,4]<-half.life.calc(fit)
weibull[i,5]<-sum(resid(fit)^2)
weibull[i,6]<-quarter.life.calc(fit)
weibull[i,7]<-tenth.life.calc(fit)
weibull[i,8]<-sum(((1*exp((-(xNA$t/weibull[i,2])^weibull[i,1]))) - xNA$Mt)^2)/length(xNA$t[!is.na(xNA$t)])
summary(fit)
weibull.fit<-exp(-(t1/weibull[i,2])^weibull[i,1])
weibull.AIC<-AIC(fit)

#Single pool 
#Method 1: Liklihood function	
singleLL = mle2(minuslogl = LLS, start = list(k = 0.5), data = list(y=Mt),method="L-BFGS-B", lower=c(0),upper=c(1000))
#add 1 to K (in AIC caculation) for estimating sigma
attr(singleLL ,"df") = attr(singleLL,"df")+1
summary(singleLL)
ssq = sum(((1*exp(-single.k[i,1]*xNA$t)) - xNA$Mt)^2)
#Extract k value
single.k[i,1]=coef(singleLL)[1]
single.k[i,2]<-sum(((1*exp(-single.k[i,1]*xNA$t)) - xNA$Mt)^2)/length(xNA$t[!is.na(xNA$t)])
single.exp.fit<-exp(-single.k[i,1]*t1)

#Asymptotic (method 1 only): M(t)=A+(1-A)*exp(-k*t)
#Liklihood function
asymLL = mle2(minuslogl = LLA, start = list(k = 0.1,A=0.2), data = list(y=Mt),method="L-BFGS-B", lower=c(0,0))#,control=list(maxint=10e6))
attr(asymLL ,"df") = attr(asymLL,"df")+1
summary(asymLL)
#Extract k & A value
asymptotic.k[i,1]=coef(asymLL)[1]
asymptotic.k[i,2]=coef(asymLL)[2]
asymptotic.k.fit<-asymptotic.k[i,2]+((1-asymptotic.k[i,2])*exp(-asymptotic.k[i,1]*t1))


#SummaryCBGB
Sample1xx<-AICctab(singleLL, asymLL, nobs=nrow(x), sort=FALSE,base=TRUE, weights=TRUE)
Sample1xx
a.aicc[i,1]=Sample1xx$dAICc[1]
a.aicc[i,2]=Sample1xx$dAICc[2]
a.aicc[i,4]=Sample1xx$weight[1]
a.aicc[i,5]=Sample1xx$weight[2]
a.aicc[i,7]=Sample1xx$AICc[1]
a.aicc[i,8]=Sample1xx$AICc[2]
a.aicc[i,10]=weibull.AIC
a.aicc[i,]

par(mfrow=c(3,6), oma=c(0,0,0,0), mar=c(2,1,1,1), mgp=c(2.4,0.8,0), cex.lab=1.25, cex.main=1.25, cex.axis=1.15)
plot(t,Mt, pch=16, col="dark grey", xlab="Yrs_Since_Deployment",ylab="Prop. C Remaining", main="CBGB Plot 1", xlim=c(0,8), ylim=c(0,1))
lines(t1,exp(-single.k[i,1]*t1), col="grey")
lines(t1, weibull.fit, col="dodgerblue1", lty=1)
lines(t1, asymptotic.k.fit, col="mediumpurple", lty=2)
legend("topright", c("single","weibull"), col=c("grey","dodgerblue1"), lty=c(1,1), bty="n")

t=(as.numeric(s1$Yrs_Since_Deployment))
Nt=s1$Fin_N_Mass_g
if (length(which(!is.na(Nt)))< 4) {
  Nt<-rep(NA, length(Nt))
}
N.length.array[i,1] = sum(!is.na(Nt))

plot(t, Nt, xlim=c(0,7))
Time <- (as.numeric(s1$Yrs_Since_Deployment))
Time2 <- (as.numeric(s1$Yrs_Since_Deployment))^2
quadratic.model<-lm(Nt~0+Time+Time2, offset=rep(Nt[1], length(Time)))
timevalues <- seq(0, 7, 0.01)
predictedcounts <- predict(quadratic.model,list(Time=timevalues, Time2=timevalues^2))
lines(timevalues, predictedcounts, col = "darkgreen", lwd = 3)
Nmat<-as.data.frame(cbind(timevalues, predictedcounts))
Nmodel[i,1]<-Nmat[Nmat$predictedcounts == max(Nmat$predictedcounts), "timevalues"]
Nmodel[i,2]<-max(Nmat$predictedcounts)

xNA$single.predict<-exp(-single.k[i,1]*xNA$t) #For figure
rmse[i,1]<-rmse(xNA$single.predict, Mt)
xNA$asymp.predict<-asymptotic.k[i,2]+((1-asymptotic.k[i,2])*exp(-asymptotic.k[i,1]*xNA$t)) #For figure
rmse[i,2]<-rmse(xNA$asymp.predict, Mt)
xNA$weibull.predict<-exp(-(xNA$t/weibull[i,2])^weibull[i,1])
rmse[i,4]<-rmse(xNA$weibull.predict, Mt)

#########
colnames(a.aicc)<-c("Single_dAICc", "Asym_dAICc", "Double_dAICc", "Single_weight", "Asym_weight", "Double_Weight", "Single_AICc", 
                    "Asym_AICc", "Double_AICc", "weibull.AIC")

k_vals<-cbind(
  single.k,
  double.k,
  asymptotic.k,
  weibull)
head(k_vals)

colnames(k_vals)<-c("single.k", "single.k.div","double.ks", "double.k","double.A", "asymp.k", "asymp.A", "double.empty", "weibull.beta", "weibull.alpha", "weibull.mrt", 
                    "weibull.half.life", "weibull.RSS", "weibull.quarter", "weibull.tenth", "weibull.div")


write.csv(a.aicc, "20210114_Plot_a.aicc_table.csv")
write.csv(k_vals, "20210114_Plot_k_vals_table.csv")
write.csv(rmse, "20210114_RMSE.csv")
write.csv(Nmodel, "20210123_Nmodel.csv")
