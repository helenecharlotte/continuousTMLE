### run.R --- 
#----------------------------------------------------------------------
## Author: Helene Charlotte Rytgaard
## Created: May, 2020 
#----------------------------------------------------------------------
## 
### Commentary:
## Run TMLE across different sample sizes and different
## numbers of time points. 
##  
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#-------------------------------------------------------------------------------------------#
## packages and generic functions
#-------------------------------------------------------------------------------------------#

library(data.table)
library(zoo)
library(stringr)
library(ltmle)
library(parallel)
library(foreach)
library(doParallel)

numextract <- function(string){ 
    as.numeric(str_extract(string, "\\-*\\d+\\.*\\d*"))
}

logit <- function(p) log(p/(1-p)) #qlogis :) 
mse <- function(x, x0=NULL) {
    if (length(x0)==0) {
        return(sqrt(mean((x-mean(x))^2)))
    } else {
        return(sqrt(mean((x-x0)^2)))
    }
}

#-------------------------------------------------------------------------------------------#
## set working directory
#-------------------------------------------------------------------------------------------#

if (system("echo $USER",intern=TRUE)%in%c("jhl781")){
    setwd("/home/ifsv/jhl781/research/phd/berkeley/continuousTMLE/")
} else {
    setwd("~/research/phd/berkeley/continuousTMLE/")
}

#-------------------------------------------------------------------------------------------#
## source relevant scripts
#-------------------------------------------------------------------------------------------#

source("./R/sim-data-continuous.R")
source("./R/fit-poisson-hal.R")

#-------------------------------------------------------------------------------------------#
## set parameters
#-------------------------------------------------------------------------------------------#

betaA <- 0.1
betaL <- 0.3
nu <- 1.3#1.1#0.9#1.2
eta <- 1#1.2#1#4/sqrt(2)*(1/8)
tau <- 3#5

#-------------------------------------------------------------------------------------------#
## true intensity
#-------------------------------------------------------------------------------------------#

true.lambda <- function(t, A, L, betaA, betaL, nu, eta) {
    return(nu*eta*t^{nu-1}*exp(A*betaA + L*betaL))
}

true.Lambda <- function(t, A, L, betaA, betaL, nu, eta) {
    return(eta*t^{nu}*exp(A*betaA + L*betaL))
}

#-------------------------------------------------------------------------------------------#
## generate data
#-------------------------------------------------------------------------------------------#

dt <- sim.data(1000, betaA=betaA, betaL=betaL, nu=nu, eta=eta, tau=tau)

#-------------------------------------------------------------------------------------------#
## fit poisson hal
#-------------------------------------------------------------------------------------------#

dt.hal <- fit.poisson.hal(dt, type=c("equal_range"), nbins=40, verbose=TRUE, lambda.cv=0.05)

#-------------------------------------------------------------------------------------------#
## compare to truth
#-------------------------------------------------------------------------------------------#

dt.hal[, true.lambda:=true.lambda(time, A, L, betaA=betaA, betaL=betaL, nu=nu, eta=eta)]
dt.hal[, true.Lambda:=true.Lambda(time, A, L, betaA=betaA, betaL=betaL, nu=nu, eta=eta)]
dt.hal[, true.density:=true.lambda*exp(-true.Lambda)]

AL <- unique(dt.hal[, c("A", "L"), with=FALSE])[sample(1:nrow(unique(dt.hal[, c("A", "L"), with=FALSE])), 4), c("A", "L"), with=FALSE]

par(mfrow=c(3,4))
#-- lambda;
plot(dt.hal[D>0  & A==AL[1,A] & L==AL[1,L], time], dt.hal[D>0  & A==AL[1,A] & L==AL[1,L], true.lambda], col="red", type="l")
lines(dt.hal[D>0  & A==AL[1,A] & L==AL[1,L], time], dt.hal[D>0  & A==AL[1,A] & L==AL[1,L], fit.lambda])
plot(dt.hal[D>0  & A==AL[2,A] & L==AL[2,L], time], dt.hal[D>0  & A==AL[2,A] & L==AL[2,L], true.lambda], col="red", type="l")
lines(dt.hal[D>0  & A==AL[2,A] & L==AL[2,L], time], dt.hal[D>0  & A==AL[2,A] & L==AL[2,L], fit.lambda])
plot(dt.hal[D>0  & A==AL[3,A] & L==AL[3,L], time], dt.hal[D>0  & A==AL[3,A] & L==AL[3,L], true.lambda], col="red", type="l")
lines(dt.hal[D>0  & A==AL[3,A] & L==AL[3,L], time], dt.hal[D>0  & A==AL[3,A] & L==AL[3,L], fit.lambda])
plot(dt.hal[D>0  & A==AL[4,A] & L==AL[4,L], time], dt.hal[D>0  & A==AL[4,A] & L==AL[4,L], true.lambda], col="red", type="l")
lines(dt.hal[D>0  & A==AL[4,A] & L==AL[4,L], time], dt.hal[D>0  & A==AL[4,A] & L==AL[4,L], fit.lambda])
#-- Lambda;
plot(dt.hal[D>0  & A==AL[1,A] & L==AL[1,L], time], dt.hal[D>0  & A==AL[1,A] & L==AL[1,L], true.Lambda], col="red", type="l")
lines(dt.hal[D>0  & A==AL[1,A] & L==AL[1,L], time], dt.hal[D>0  & A==AL[1,A] & L==AL[1,L], fit.Lambda])
plot(dt.hal[D>0  & A==AL[2,A] & L==AL[2,L], time], dt.hal[D>0  & A==AL[2,A] & L==AL[2,L], true.Lambda], col="red", type="l")
lines(dt.hal[D>0  & A==AL[2,A] & L==AL[2,L], time], dt.hal[D>0  & A==AL[2,A] & L==AL[2,L], fit.Lambda])
plot(dt.hal[D>0  & A==AL[3,A] & L==AL[3,L], time], dt.hal[D>0  & A==AL[3,A] & L==AL[3,L], true.Lambda], col="red", type="l")
lines(dt.hal[D>0  & A==AL[3,A] & L==AL[3,L], time], dt.hal[D>0  & A==AL[3,A] & L==AL[3,L], fit.Lambda])
plot(dt.hal[D>0  & A==AL[4,A] & L==AL[4,L], time], dt.hal[D>0  & A==AL[4,A] & L==AL[4,L], true.Lambda], col="red", type="l")
lines(dt.hal[D>0  & A==AL[4,A] & L==AL[4,L], time], dt.hal[D>0  & A==AL[4,A] & L==AL[4,L], fit.Lambda])
#-- density;
plot(dt.hal[D>0  & A==AL[1,A] & L==AL[1,L], time], dt.hal[D>0  & A==AL[1,A] & L==AL[1,L], true.density], col="red", type="l")
lines(dt.hal[D>0  & A==AL[1,A] & L==AL[1,L], time], dt.hal[D>0  & A==AL[1,A] & L==AL[1,L], fit.density])
plot(dt.hal[D>0  & A==AL[2,A] & L==AL[2,L], time], dt.hal[D>0  & A==AL[2,A] & L==AL[2,L], true.density], col="red", type="l")
lines(dt.hal[D>0  & A==AL[2,A] & L==AL[2,L], time], dt.hal[D>0  & A==AL[2,A] & L==AL[2,L], fit.density])
plot(dt.hal[D>0  & A==AL[3,A] & L==AL[3,L], time], dt.hal[D>0  & A==AL[3,A] & L==AL[3,L], true.density], col="red", type="l")
lines(dt.hal[D>0  & A==AL[3,A] & L==AL[3,L], time], dt.hal[D>0  & A==AL[3,A] & L==AL[3,L], fit.density])
plot(dt.hal[D>0  & A==AL[4,A] & L==AL[4,L], time], dt.hal[D>0  & A==AL[4,A] & L==AL[4,L], true.density], col="red", type="l")
lines(dt.hal[D>0  & A==AL[4,A] & L==AL[4,L], time], dt.hal[D>0  & A==AL[4,A] & L==AL[4,L], fit.density])
        
