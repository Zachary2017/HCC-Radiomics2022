#ICC
T1 <- read.csv("~R#1_1.csv",header = T)
T2 <- read.csv("~R#1_2.csv",header = T)
T3 <- read.csv("~R#2.csv",header = T)
dim(T1);dim(T2);dim(T3)
x <- dim(T1)[1]
y <- dim(T1)[2]
T12 <- cbind(T1,T2);dim(T12)
T13 <- cbind(T1,T3);dim(T13)

library(psych)
t= 2
icc <- c(1:y)
for(i in 1:y) {icc[i] <- ICC(T12[,c(i,i+y)])$results$ICC[t]}
icc
mean(icc);median(icc)
l <- length(which(icc >= 0.75)); 
s <- length(which(icc < 0.75)); 
icc1 <- c(1:y)
for(i in 1:y) {icc1[i] <- ICC(T13[,c(i,i+y)])$results$ICC[t]}
icc1
mean(icc1);median(icc1)
l1 <- length(which(icc1 >= 0.75)); 
l1

#LASSO
library(pROC)
library(glmnet)
m <- read.csv("~train.csv",header = T)
x<- as.matrix(m[,-1])
y<- as.double(as.matrix(m[, 1]))
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

#DCA
library(rmda)
train <- read.csv("~train_combined.csv",header = T)
validation <- read.csv("~validation_combined.csv",header = T)
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