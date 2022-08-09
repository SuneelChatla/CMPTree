### required libraries
rm(list = ls())

# loading source cmp
require(mgcv)
loadNamespace("mgcv")
require(MASS)
source("cmp.R")
#########################

## reading bike-sharing data
hr=read.csv("Bike-Sharing-Dataset/hour.csv",header = T)
head(hr)

########### variable coding
hr$season=with(hr,factor(season,levels = c(1,2,3,4),
                         labels=c("springer", "summer", "fall", "winter")))
hr$yr=with(hr,factor(yr,levels = c(0,1),labels=c( "2011", "2012")))
hr$mnth=with(hr,factor(mnth,levels=c(1:12),
                       labels=c("jan","feb","mar","apr","may","jun",
                                "jul","aug","sep","oct","nov","dec")))
hr$holiday=with(hr,factor(holiday,levels=c(0,1),
                          labels=c("not holiday","holiday")))
hr$weekday=with(hr,factor(weekday,levels = c(0:6),
                          labels = c("sun","mon","tue","wed","thu","fri","sat")))
hr$workingday=with(hr, factor(workingday,levels = c(0,1),
                              labels = c("not working day","working day")))
hr$weathersit=with(hr,factor(weathersit,levels = c(1,2,3,4),
                             labels=c("clear","cloudy","light rain/fog","heavy rain/fog")))
# hr$hr=with(hr,factor(hr,levels = c(0:23),
#                      labels = c("0","1","2","3","4","5","6","7","8","9","10","11",
#                                 "12","13","14","15","16","17","18","19","20","21","22","23"),ordered = TRUE))

hr$day=as.integer(substr(hr$dteday,9,10))

# factors to binary
weathersit=sapply(unique(hr$weathersit), function(x) as.integer(x == hr$weathersit))
colnames(weathersit)[1:4]=c("clear", "cloudy" ,"lightrain","heavyrain")
holiday=sapply(levels(hr$holiday), function(x) as.integer(x == hr$holiday))
colnames(holiday)=c("notholiday","holidayc")
weekday=sapply(levels(hr$weekday), function(x) as.integer(x == hr$weekday))
##
hr=cbind(hr,weathersit,holiday,weekday)
####### Data for a specific month
nov12=subset(hr,yr=="2012"& mnth=="nov")
## For tuning parameters
nov11=subset(hr,yr=="2011"& mnth=="nov")
## For prediction
dec12=subset(hr,yr=="2012"& mnth=="dec")


################################################################################
############ GLM for registered and casual users #############################
################################################################################
regform.cas <- as.formula(casual~day+hr+holiday+weekday+weathersit+
                   atemp+hum+windspeed)

regform.reg <- as.formula(registered~day+hr+holiday+weekday+weathersit+
                    atemp+hum+windspeed)

#
rmse <- function(y1,y2)
{
  if(length(y1) != length(y2)) stop("both vectors are not of same size")
  sqrt(mean((y1-y2)^2))
}
#devtools::install_github("thomas-fung/mpcmp")
#library(mpcmp)

####### Casual users
cmp.reg <- gam.cmp(regform.cas,family = cmp,data = nov12)
cmp.pred <- predict.gam.cmp(cmp.reg,dec12,type="response")$fit
rmse(dec12$casual,cmp.pred)

#
pos.reg <- glm(regform.cas,family = poisson(),data = nov12)
pos.pred <- predict(pos.reg,dec12, type="response")
rmse(dec12$casual, pos.pred)
#
nb.reg <- glm.nb(regform.cas,data = nov12)
nb.pred <- predict(nb.reg, dec12, type="response")
rmse(dec12$casual, nb.pred)
#
mpcmp.reg <-mpcmp::glm.cmp(regform.cas,data = nov12)
mp.pred <- predict(mpcmp.reg, dec12, type="response")
rmse(dec12$casual,mp.pred)

# boxplot
pdf("PEbox-reg.pdf",9,5)
par(mfrow=c(1,2))
boxplot(dec12$casual-cmp.pred,dec12$casual-pos.pred, dec12$casual-nb.pred, dec12$casual-mp.pred,
        main="Prediction errors for Casual users",
        names=c("CMP","Poi", "NB", "MPCMP"),
        las=2,
        horizontal = TRUE,
        notch = TRUE)

######### registered users

cmp.reg2 <- gam.cmp(regform.reg,family = cmp,data = nov12)
cmp.pred2 <- predict.gam.cmp(cmp.reg2,dec12,type="response")$fit
rmse(dec12$registered,cmp.pred2)
#
pos.reg2 <- glm(regform.reg,family = poisson(),data = nov12)
pos.pred2 <- predict(pos.reg2,dec12, type="response")
rmse(dec12$registered, pos.pred2)
#
nb.reg2 <- glm.nb(regform.reg,data = nov12)
nb.pred2 <- predict(nb.reg2, dec12, type="response")
rmse(dec12$registered, nb.pred2)
#
mpcmp.reg2 <-mpcmp::glm.cmp(regform.reg,data = nov12)
mp.pred2 <- predict(mpcmp.reg2, dec12, type="response")
rmse(dec12$registered,mp.pred2)


# boxplot
boxplot(dec12$registered-cmp.pred2,dec12$registered-pos.pred2, dec12$registered-nb.pred2, dec12$registered-mp.pred2,
        main="Prediction errors for Registered users",
        names=c("CMP","Poi", "NB", "MPCMP"),
        las=2,
        horizontal = TRUE,
        notch = TRUE)

dev.off()



######################################################################################
############################### GAM #################################################
######################################################################################

gamform.cas <-as.formula(casual~day+holiday+weekday+weathersit
                 +hum+windspeed +s(hr)+ s(atemp))

##### Casual users
cmp.gam <-gam.cmp(gamform.cas,family = cmp,data = nov12)
cmpgam.pred <- predict.gam.cmp(cmp.gam,dec12,type="response")$fit
rmse(dec12$casual,cmpgam.pred)
#
pos.gam <- gam(gamform.cas,family = poisson(),data = nov12)
posgam.pred <- predict(pos.gam,dec12,type="response")
rmse(dec12$casual,posgam.pred)
#
nb.gam <- gam(gamform.cas,data = nov12,family = nb())
nbgam.pred <- predict(nb.gam,dec12,type="response")
rmse(dec12$casual,nbgam.pred)
#
pdf("gam-cas.pdf",7,4)
plot(cmp.gam,pages=1)
dev.off()

## predictions
pdf("PEbox-gam.pdf",9,4)
par(mfrow=c(1,2))
boxplot(dec12$casual-cmpgam.pred,dec12$casual-posgam.pred, dec12$casual-nbgam.pred,
        main="Prediction errors for Casual users",
        names=c("CMPGAM","PoiGAM", "NBGAM"),
        las=2,
        horizontal = TRUE,
        notch = TRUE)




##### Registered users
gamform.reg <-as.formula(registered~day+holiday+weekday+weathersit
                 +hum+windspeed +s(hr)+ s(atemp))
#
cmp.gam2 <-gam.cmp(gamform.reg,family = cmp,data = nov12)
cmpgam.pred2 <- predict.gam.cmp(cmp.gam2,dec12,type="response")$fit
rmse(dec12$registered,cmpgam.pred2)
#
pos.gam2 <- gam(gamform.reg,family = poisson(),data = nov12)
posgam.pred2 <- predict(pos.gam2, dec12,type="response")
rmse(dec12$registered,posgam.pred2)
#
nb.gam2 <- gam(gamform.reg,data = nov12,family = nb())
nbgam.pred2 <- predict(nb.gam2, dec12,type="response")
rmse(dec12$registered,nbgam.pred2)

#
pdf("gam-reg.pdf",9,4)
plot(cmp.gam2,pages=1)
dev.off()

## prediction

boxplot(dec12$registered-cmpgam.pred2,dec12$registered-posgam.pred2, dec12$registered-nbgam.pred2,
        main="Prediction errors for Registered users",
        names=c("CMPGAM","PoiGAM", "NBGAM"),
        las=2,
        horizontal = TRUE,
        notch = TRUE)

dev.off()




##############################################################################################
#################### Boosting #########################################################
##############################################################################################

formula <- list(as.formula(casual~ hum),as.formula(casual~ windspeed)) #+hr+day
fformula <- list(as.formula(~s(hr) + s(atemp)+weekday),as.formula(~day))
vcformula <- list(as.formula(~   windspeed+ clear + cloudy + lightrain + heavyrain
  + notholiday + holidayc),as.formula(~ clear + cloudy + lightrain + heavyrain +
  notholiday + holidayc +  sun + mon + tue + thu + fri + sat ))


##### Boosting
### Tuning the terminal nodes and the number of iterations
grid.M <- c(15,20,25)
nboost <- seq(1,40,3)
prmse <- matrix(0,length(nboost),length(grid.M))

for(i in 1:length(nboost))
{
  for(j in 1:length(grid.M))
  {
  tres <-   Boost.Tree(formula=formula,fformula=fformula,vcformula=vc1formula,data=nov12,eta=0.01,mini.size=20,M=grid.M[j],nBoost = nboost[i])
  #
  tpred <-predict.boost.cmp(tres,nov11,type = "response")
  prmse[i,j] <- rmse(nov11$casual,tpred)
  }
  cat("Iter:",i ,"\n")
}

### plotting
plot(nboost[1:10], prmse[1:10,1],xlab = "Boosting iterations",ylab = "RMSE", type = "l",lwd=1,ylim = c(106,109))
lines(nboost[1:10], prmse[1:10,2],lty=2)
lines(nboost[1:10], prmse[1:10,3],lty=3)
#lines(nboost[1:10], prmse[1:10,4],lty=4)
legend("top", legend=c("M=15","M=20", "M=25"), lty=c(1,2,3),lwd=2)

###### CMPBoost model fit
### For casual : M=25, nboost=19
### For registered: M=20, nboost=22
cmp.boost <- Boost.Tree(formula = formula,fformula = fformula,vcformula = vcformula,data = nov12, eta = 0.01, M=15,nBoost = 19)
cmpbt1.pred <- predict.boost.cmp(cmp.boost,pred,type="response")
rmse(pred$registered,cmpbt1.pred)

###
## variable importance plot
par(mar = c(8, 3, 0.5, 0.5))
barplot(cmp.bt1$variable.imp.lam/sum(cmp.bt1$variable.imp.lam),names.arg = names(cmp.bt1$varPart1),las=2,cex.names = 1.5)

### partial plot
vnames=c("lightrain")
part.boost <- Partial.Boost(cmp.boost,vnames)

## partial plots for CMPBoost
par(mar = c(5, 5, 0.5, 0.5))
plot(pat.bt1$ulevel.lam[,1][1:2],pat.bt1$coefs.lam[[1]][,2],pch=19,xlab = "lightrain",ylab=expression(beta[1](windspeed)),cex.lab=1.5)
points(pat.bt1$ulevel.lam[,1][1:2],pat.bt1$coefs.lam[[1]][,2],type="l",lty=2)


######################################################################
###################### CMPMOB Trees ##################################
######################################################################

### input model formulation
formula <- list(as.formula(casual~ hum+weekday),as.formula(casual~ hum)) #+hr+day
fformula <- list(as.formula(~s(hr)),as.formula(~atemp+day))
vcformula <- list(as.formula(~   atemp+ hum+ windspeed+ weathersit  + holiday),
                  as.formula(~ atemp+ hum+ windspeed+ weathersit  + holiday ))
### CMP MOB tree
cmp.mob <- cmpmob(formula,fformula,vcformula,data = nov12,family = cmp,sflag = TRUE)
cmpmob.pred <- predict.mob.cmp(cmp.mob, dec12, type="response")
rmse(dec12$casual,cmpmob.pred)

### plot for the tree
plot(cmp.mob$tree)

### partial residual plot for the global effect
plot.gam(cmp.mob$fixed,rug=TRUE)

### Finding the local and global deviance.
node2 <- nodeids(cmp.mob$tree,terminal = TRUE)
nterm2 <- length(node2)
ldev2 <- sum(unlist(nodeapply(cmp.mob$tree, ids = node2, FUN = function(n) info_node(n)$deviance)))
gdev2 <- cmp.mob$glob$deviance



#########################################################################
##################### Plots #############################################
#########################################################################
require(ggplot2)
require(purrr)
require(dplyr)

## hour
# casual
p1 <- hr11 %>%
  group_by(hr) %>%
  summarise(mean_cas=mean(casual), sd_cas=sd(casual), n=n()) %>%
  mutate(se=sd_cas/sqrt(n)) %>%
  ggplot() +
  geom_linerange( aes(x=hr, ymin=mean_cas-se, ymax=mean_cas+se), colour="skyblue", alpha=0.9, size=1.3)+
  geom_bar(aes(x=hr,y=mean_cas), stat="identity", fill="orange", alpha=0.6)+
  xlab("Hour of the day") + ylab("Average number of rides") + ggtitle("Casual users")

# registered
p2 <- hr11 %>%
  group_by(hr) %>%
  summarise(mean_reg=mean(registered), sd_reg=sd(registered), n=n()) %>%
  mutate(se=sd_reg/sqrt(n)) %>%
  ggplot() +
  geom_linerange( aes(x=hr, ymin=mean_reg-se, ymax=mean_reg+se), colour="skyblue", alpha=0.9, size=1.3)+
  geom_bar(aes(x=hr,y=mean_reg), stat="identity", fill="orange", alpha=0.6)+
  xlab("Hour of the day") + ylab("Average number of  rides") + ggtitle("Registered users")
#
library(gridExtra)
pdf("hour-dec.pdf",7,4)
grid.arrange(p1,p2,nrow = 1)
dev.off()
## Day
# casual
p3 <- hr11 %>%
  group_by(day) %>%
  summarise(mean_cas=mean(casual), sd_cas=sd(casual), n=n()) %>%
  mutate(se=sd_cas/sqrt(n)) %>%
  ggplot() + geom_linerange( aes(x=day, ymin=mean_cas-se, ymax=mean_cas+se), colour="skyblue", alpha=0.9, size=1.3)+
  geom_bar(aes(x=day,y=mean_cas), stat="identity", fill="orange", alpha=0.6)+
  xlab("Day") + ylab("Average number of rides") + ggtitle("Casual users")

# registered
p4 <- hr11 %>%
  group_by(day) %>%
  summarise(mean_reg=mean(registered), sd_reg=sd(registered), n=n()) %>%
  mutate(se=sd_reg/sqrt(n)) %>%
  ggplot() +geom_linerange( aes(x=day, ymin=mean_reg-se, ymax=mean_reg+se), colour="skyblue", alpha=0.9, size=1.3)+
  geom_bar(aes(x=day,y=mean_reg), stat="identity", fill="orange", alpha=0.6)+
  xlab("Day") + ylab("Average number of  rides") + ggtitle("Registered users")
#
pdf("day-dec.pdf",7,4)
grid.arrange(p3,p4,nrow = 1)
dev.off()

## weathersit
# casual
p5 <- hr11 %>%
  group_by(weathersit) %>%
  summarise(mean_cas=mean(casual), sd_cas=sd(casual), n=n()) %>%
  mutate(se=sd_cas/sqrt(n)) %>%
  ggplot() +geom_linerange( aes(x=weathersit, ymin=mean_cas-se, ymax=mean_cas+se), colour="skyblue", alpha=0.9, size=1.3)+
  geom_bar(aes(x=weathersit,y=mean_cas), stat="identity", fill="orange", alpha=0.6)+
  xlab("Weather type") + ylab("Average number of rides") + ggtitle("Casual users")

# registered
p6 <- hr11 %>%
  group_by(weathersit) %>%
  summarise(mean_reg=mean(registered), sd_reg=sd(registered), n=n()) %>%
  mutate(se=sd_reg/sqrt(n)) %>%
  ggplot() +geom_linerange( aes(x=weathersit, ymin=mean_reg-se, ymax=mean_reg+se), colour="skyblue", alpha=0.9, size=1.3)+
  geom_bar(aes(x=weathersit,y=mean_reg), stat="identity", fill="orange", alpha=0.6)+
  xlab("Weather type") + ylab("Average number of  rides") + ggtitle("Registered users")

#
pdf("weather-dec.pdf",7,4)
grid.arrange(p5,p6,nrow = 1)
dev.off()
## Temperature
# casual
p7 <- hr11 %>%
  ggplot(aes(x=temp,y=casual)) + geom_point(aes(colour=weekday)) +
  xlab("Temperature") + ylab("Number of rides") + ggtitle("Casual users")

# registered
p8 <- hr11 %>%
  ggplot(aes(x=temp,y=registered)) + geom_point(aes(colour=weekday)) +
  xlab("Temperature") + ylab("Number of rides") + ggtitle("Registered users")

### windspeed

p9 <- hr11 %>%
  ggplot(aes(x=windspeed,y=casual)) + geom_point(aes(colour=weathersit,shape=weathersit)) +
  xlab("Wind speed") + ylab("Number of rides") + ggtitle("Casual users") +
  guides(col=guide_legend(title="Weather type"), shape=guide_legend(title="Weather type"))

# registered
p10 <- hr11 %>%
  ggplot(aes(x=windspeed,y=registered)) + geom_point(aes(colour=weathersit,shape=weathersit)) +
  xlab("Wind speed") + ylab("Number of rides") + ggtitle("Registered users")+
  guides(col=guide_legend(title="Weather type"), shape=guide_legend(title="Weather type"))

pdf("windspeed-dec.pdf",9,4)
grid.arrange(p9,p10,nrow = 1)
dev.off()

### Humidity

p11 <- hr11 %>%
  ggplot(aes(x=atemp,y=casual)) + geom_point(aes(colour=weathersit,shape=weathersit)) +
  xlab("Temperature") + ylab("Number of rides") + ggtitle("Casual users")+
  guides(col=guide_legend(title="Weather type"), shape=guide_legend(title="Weather type"))

# registered
p12 <- hr11 %>%
  ggplot(aes(x=atemp,y=registered)) + geom_point(aes(colour=weathersit,shape=weathersit)) +
  xlab("Temperature") + ylab("Number of rides") + ggtitle("Registered users")+
  guides(col=guide_legend(title="Weather type"), shape=guide_legend(title="Weather type"))

pdf("temp-dec.pdf",9,4)
grid.arrange(p11,p12,nrow = 1)
dev.off()



