## Model Fit Analysis for table S2
urlfile = "https://raw.githubusercontent.com/gill20a/NutNetLitterDecomp/master/ModelFitting/SiteTrt_ModelOutput"
data<-read.csv(url(urlfile))
head(data)

##Model with lowest AIC (single, double, asymp)
tib2<-as.data.frame(data %>%  count(Best_Model_3))
tib2$prob<-tib2$n/sum(tib2$n)*100
tib2
##Table of best models within 3 AIC units (single, double, asymp)
tib2<-as.data.frame(data %>%  count(Best_Model_3all))
tib2$prob<-tib2$n/sum(tib2$n)*100
tib2
##Model with lowest AIC (single, double, asymp, weibull)
tib2<-as.data.frame(data %>%  count(Best_Model_4))
tib2$prob<-tib2$n/sum(tib2$n)*100
tib2
##Table of best models within 3 AIC units (single, double, asymp)
tib2<-as.data.frame(data %>%  count(Best_Model_4all))
tib2$prob<-tib2$n/sum(tib2$n)*100
tib2

#Table S3
tib2<-as.data.frame(data %>% group_by(Treatment) %>%  count(Best_Model_4))
tib2$prop<-(tib2$n/20)*100
tib2

tib2<-as.data.frame(data %>% group_by(Site) %>%  count(Best_Model_4))
tib2$prop<-(tib2$n/8)*100
tib2
