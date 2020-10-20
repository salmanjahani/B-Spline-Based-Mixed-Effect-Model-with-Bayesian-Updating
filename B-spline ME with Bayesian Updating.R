#' Code for B-Spline Based Mixed Effects Model With Bayesian Updating 
#' Simulation
list.of.packages <- c("Matrix","nlme","rootSolve","splines")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load the packages
library(Matrix)                  
library(nlme)                    # Nonlinear Mixed Effects Model
library(rootSolve)               # Nonlinear Root Finding
library(splines)                 # Regression spline functions and classes

rm(list=ls())                    # Clean-up the memory
########################################### System Input ##############################################

setwd("D:/Job Hunt/Sample Code/B-Spline Based Mixed Effect Model with Bayesian Updating")  #please set the correct directory 
load("IMdata.RData")


######################################### Construct a groupedData Object #################################################
data <- data.frame(id = numeric(0), x = numeric(0), y = numeric(0))
for (i in 1:200) {
  x <- seq(1,11)
  y <- cn[[4]][i,]               #Choosing the 4th signal (Compressor Discharge Temperature) in the dataset
  id <- rep(i, 11)
  data <- rbind(data, cbind(id = id, x = x, y = y))
}
data <- groupedData(y ~ x | id, data)
########################################## Train Test Mean ###############################################
n=123                             # first n-1  available signals are selected to train the model
# The n th signal is considered the new signal for which we are interested to do prediction
train=data[data$id%in%(1:n-1),]
# test=data[data$id%in%n,]
trains=lapply(1:(n-1),function(i){train$x[train$id==i]})
trainy=lapply(1:(n-1),function(i){train$y[train$id==i]})
########################################## Exploratory data visualization ###############################################
dev.off()                        #Remove any previous plot
for(i in 1:(20))                 #Plot the first 20 Evaporator outlet temperature signals
{
  plot(trains[[i]],trainy[[i]],xlim=c(0,11),ylim=c(0,200),type="l",pch=2,lwd=2,cex.axis=1.2,xlab='Time',ylab='Temperature (fahrenheit)',main='Compressor Discharge Temperature')
  par(new=T)
}
# for(i in 1:(20))                 #Plot the first 20 Evaporator outlet temperature signals
# {
#   plot(trains[[4]],80*bs(x, df=8,degree=4)[,i],xlim=c(0,11),ylim=c(0,200),type="l",pch=2,lwd=2,cex.axis=1.2,xlab='Time',ylab='Temperature (fahrenheit)',main='Compressor Discharge Temperature')
#   par(new=T)
# }

##################################### B-Spline Mixed effects model ##########################
#Fit B-spline mixed effect model
fitlme <- tryCatch(lme(y~bs(x, df=8,degree=4),random=~bs(x, df=8,degree=4)|id,data=train,control=lmeControl(returnObject=TRUE,opt="optim",optimMethod = "SANN")), error = function(e) e)
if(any(class(fitlme) == "error")==T){cat("LME fit Prob")}  #Check for any error in fitting
# Check the parameter estimates
summary(fitlme)


# Extract parameters
sigma2f = (fitlme$sigma)^2                      # Signal noise
MUBf = fitlme$coefficients$fixed                # Mean vector
SIGMABf = var(fitlme$coefficients$random$id)    # Covariance matrix
#####Online Update#####
####Similar to Junbo, I use a set of randomly-generated signal for the sake of illustration
new_t = rmnorm(1,MUBf,SIGMABf)                  # b vector for this hypothetical new unit
noisevec = rnorm(11,0,sigma2f)                  # Siganl noise
ftt0=lm(y~bs(x,df=8,degree=4),data=train)       #Dummy model
ftt0$coefficients=as.vector(new_t)
test=predict(ftt0,newdata = data.frame(x))+noisevec
test=data.frame(x,y=test)

dev.off()
par(mfrow=c(1,2))

for(i in 1:11){
  tstar=i                       # Time of Bayesian Update
  tests=test$x[test$x<=tstar];m1=length(tests);testy=test$y[1:m1]  # Test Signal
  Rp=testy          # new observed temperature signal
  Rp = t(t(Rp))    
  Zp = matrix(0,m1,9)
  Zp[,1]=rep(1,m1)
  Zp[,(2:9)]=bs(x, df=8,degree=4)[1:m1,]
  
  ########### Bayesian parameters' update  (Posterior Distribution)
  a <- tryCatch(chol2inv(chol(SIGMABf)), error = function(e) e) #a=solve(SIGMABf)
  if(any(class(a) == "error")==T){cat("No inverse")}
  b=a+t(Zp)%*%Zp/sigma2f
  SIGMABpf <- tryCatch(chol2inv(chol(b)), error = function(e) e)
  if(any(class(SIGMABpf) == "error")==T){cat("No inverse");iit=iit-1;next}
  MUBpf=SIGMABpf%*%(t(Zp)%*%Rp/sigma2f+a%*%MUBf)
 
  #### Comparison of Prior and Posterior
  comparison = cbind(MUBf,MUBpf)
  colnames(comparison) = c("Prior","Bayesian Updated")
  comparison
  
  ########## Visual illustration #####
  ftt=lm(y~bs(x,df=8,degree=4),data=train)        #Dummy model to be used for b-spline predictions
  ####Prior Fit
  ftt$coefficients=as.vector(MUBf)
  priorfit=predict(ftt,newdata = data.frame(x))
  ####Posterior Fit
  ftt$coefficients=as.vector(MUBpf)
  posteriorfit=predict(ftt,newdata = data.frame(x))
  ##True Signal
  # dev.off()
  plot(test$x,test$y,xlim=c(0,11),ylim=c(90,160),type="l",lwd=2,xlab='Time',ylab='temperature (fahrenheit)',main='New Compressor Discharge Temperature Signal',cex.axis=1.2) ## True
  ##Current observations of new signal
  points(tests,testy,xlim=c(0,11),ylim=c(90,160),pch=16,cex=2,xlab="",ylab="",main="",cex.axis=1.2)
  par(new=TRUE)
  ## Posterior and Prior predictions
  plot(test$x,posteriorfit,type="l",lwd=2,lty=2,col="red",xlab="",ylab="",main="",xlim=c(0,11),ylim=c(90,160),cex.axis=1.2)
  par(new=TRUE)
  plot(test$x,priorfit,type="l",lwd=2,lty=3,col="blue",xlab="",ylab="",main="",xlim=c(0,11),ylim=c(90,160),cex.axis=1.2)
  # legend("bottomright",col=c("black","blue","red"),legend=c("True","Prior","Posterior"),lty=1:3,lwd=3.5,cex=1.2)
}