#Z score
data <- read.csv("~data.csv",header = T)
data_scale = scale(data,center = F, scale = TRUE)
data_scale = data.frame(data_scale)
write.csv(data_scale,"~data_scale.csv")



#ICC
library(psych)
T1 <- read.csv("~R#1_1.csv",header = T)
T2 <- read.csv("~R#1_2.csv",header = T)
T3 <- read.csv("~R#2.csv",header = T)
x <- dim(T1)[1]
y <- dim(T1)[2]
T12 <- cbind(T1,T2);dim(T12)
T13 <- cbind(T1,T3);dim(T13)

t= 2
icc1 <- c(1:y)
for(i in 1:y) {icc1[i] <- ICC(T12[,c(i,i+y)])$results$ICC[t]}
icc1
mean(icc1);median(icc1)
l1 <- length(which(icc1 >= 0.75)); 
s1 <- length(which(icc1 < 0.75)); 
icc2 <- c(1:y)
for(i in 1:y) {icc2[i] <- ICC(T13[,c(i,i+y)])$results$ICC[t]}
icc2
mean(icc2);median(icc2)
l2 <- length(which(icc2 >= 0.75)); 



#Comparison
train <- read.csv("~train.csv",header = T)
ps  <- apply(train[-1],2,function(x)  shapiro.test(x)$p.value)
pb  <- apply(train[-1],2,function(x)  bartlett.test(x,train$MTM)$p.value)
w1  <- which(ps>0.05 & pb>0.05)
train1  <- train[,names(w1)]

pt <- apply(train1,2,function(x)  t.test(x,train$MTM)$p.value)
fdr1 = p.adjust(pt, "BH")
w2  <- which(fdr1<0.05) 
l1  <- length(w2);l1 
A1  <- train[!(names(train) %in% names(train1))]
pwm <- apply(A1[-1],2,function(x)  wilcox.test(x,A1$MTM)$p.value)
fdr2 = p.adjust(pwm, "BH")
w3  <- which(fdr2<0.05) 
l2  <- length(w3);l2

(a <- names(sort(c(w2,w3))))
A2 <- train[a];dim(A2)  
A3 <- cbind(train[,1],A2)



#LASSO
library(pROC)
library(glmnet)
train <- read.csv("~train.csv",header = T)
x<- as.matrix(train[,-1])
y<- as.double(as.matrix(train[, 1]))
set.seed(999)
fit <- cv.glmnet(x, y, family='binomial', alpha=1, parallel=TRUE,standardize=TRUE, type.measure='auc')
fit
plot(fit,ylim = c(0.35,0.60))
plot(fit$glmnet.fit, xvar="lambda", label=TRUE)
fit$lambda.min
fit$lambda.1se
coef(fit, s=fit$lambda.min)
fit.best <- fit$glmnet.fit 
fit.coef <- coef(fit$glmnet.fit, s = fit$lambda.min) 
fit.coef[which(fit.coef != 0)] 
summary(fit.coef)
coef = coef(fit, s = fit$lambda.min) 
index = which(coef != 0) 
actCoef = coef[index] 
lassoFeature = row.names(coef)[index] 
FeatureCoef = cbind(Feature=lassoFeature,Coef=actCoef) 
FeatureCoef



#SVM
library(pROC)
train <- read.csv("~train.csv",header = T)
train[,1] = as.factor(train[,1])  
performance<- function(table,n=2){
  if(!all(dim(table)==c(2,2)))
    stop("Must be a 2*2 table")
  tn=table[1,1]
  fn=table[2,1]
  tp=table[2,2]
  fp=table[1,2]
  sensitivity=tp/(tp+fn)
  specificity=tn/(tn+fp)
  ppp=tp/(tp+fp)
  npp=tn/(tn+fn)
  hitrate=(tp+tn)/(tp+tn+fp+fn)
  F1=2*sensitivity*ppp/(ppp+sensitivity)
  result<- rbind(sensitivity,specificity,ppp,npp,hitrate,F1)
  rownames(result)<- c("sensitivity","specificity","positivive predictive value","negtive predictive value","accuracy","F1")
  colnames(result)<- c("model")
  return(result)
}

library(e1071)
tuned<- tune.svm(MTM~.,data = train,gamma = 10^(-3:0),cost = 10^(-0:3))
fit.svm<- svm(MTM~.,data = train,gamma=0.1,cost=10)
fit.svm
pred.svm.train<- predict(fit.svm,na.omit(train))
pref.svm.train<-table(na.omit(train$MTM),pred.svm.train,dnn=c("Actual","Predicted"))
pref.svm.train
per.svm.train<-performance(pref.svm.train)
pre.svm.train<- predict(fit.svm,train,decision.values=TRUE)
pre.svmprobs.train<-attr(pre.svm.train,"decision.values")
preroc.svm.train <- roc(train$MTM,pre.svmprobs.train[,1])
ci_de.svm.train<-ci(preroc.svm.train,method="delong")
ci_de.svm.train
round(ci_de.svm.train,3)
plot(preroc.svm.train,
     print.auc=T,
     print.thres=T,
     auc.polygon=T,
     max.auc.polygon=T,
     grid=c(0.2), 
     grid.col=c("green", "red"), 
     auc.polygon.col="skyblue",
     main="AUC", 
     col="1")



#nomogram
library(foreign)
library(rms)
library(forestmodel)
train <- read.csv('~train.csv', header = T)
train <- as.data.frame(train)
head(train)
train$MTM <- ifelse(train$MTM=="1",1,0)

attach(train)
dd<-datadist(train)
options(datadist='dd')

fit1<-glm(MTM~age+gender+viral+afp+size+bclc+rs,data=train,family = binomial())
fit1
fit.result<-summary(fit1)
forest_model(fit1)
car::vif(fit1)

fit2<-lrm(MTM~afp+rs,data=train,x=T,y=T)
fit2
nom <- nomogram(fit2, fun=plogis,fun.at=c(.001, .01, .05,seq(.1,.9, by=.1), .95, .99, .999),lp=F,funlabel="MTM")  
plot(nom)
validate(fit2,method="boot",B=1000,dxy=T)
cal <- calibrate(fit2, method='boot', B=1000)
plot(cal,xlim=c(0,1.0),ylim=c(0,1.0))

library(AICcmodavg)
AICc(fit2, return.K = FALSE, second.ord = TRUE,nobs = NULL)



#DCA
library(rmda)
train <- read.csv("~train.csv",header = T)
validation <- read.csv("~validation.csv",header = T)
Radiomics <- decision_curve(MTM~rs,data = train, 
                            family = binomial(link ='logit'),thresholds= seq(0,1, by = 0.01), 
                            confidence.intervals =0.95,study.design = 'case-control',population.prevalence = 0.5)
Combined<- decision_curve(MTM~afp+rs,data = train, 
                          family = binomial(link ='logit'),thresholds= seq(0,1, by = 0.01), 
                          confidence.intervals =0.95,study.design = 'case-control',population.prevalence = 0.5)
List<- list(Radiomics,Combined)
plot_decision_curve(List,curve.names= c('Radiomics signature','Combined model'),
                    cost.benefit.axis = FALSE,col = c('red','blue'),
                    confidence.intervals = FALSE,standardize = FALSE, lty = c(2,1),
                    lwd = c(3,2, 2),legend.position = "topright")
summary(Radiomics,measure= 'NB')
summary(Combined,measure= 'NB')

Radiomics1 <- decision_curve(MTM~rs,data = validation, 
                             family = binomial(link ='logit'),thresholds= seq(0,1, by = 0.01), 
                             confidence.intervals =0.95,study.design = 'case-control',population.prevalence = 0.5)
Combined1<- decision_curve(MTM~afp+rs,data = validation, 
                           family = binomial(link ='logit'),thresholds= seq(0,1, by = 0.01), 
                           confidence.intervals =0.95,study.design = 'case-control',population.prevalence = 0.5)
List<- list(Radiomics1,Combined1)
plot_decision_curve(List,curve.names= c('Radiomics signature','Combined model'),
                    cost.benefit.axis = FALSE,col = c('red','blue'),
                    confidence.intervals = FALSE,standardize = FALSE, lty = c(2,1),
                    lwd = c(3,2, 2),legend.position = "topright")
summary(Radiomics1,measure= 'NB')
summary(Combined1,measure= 'NB')