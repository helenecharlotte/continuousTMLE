### show.R --- 
#----------------------------------------------------------------------
## Author: Helene Charlotte Rytgaard
## Created: May, 2020 
#----------------------------------------------------------------------
## 
### Commentary:
## Look at output from TMLE estimation.
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

source("./R/sim-data.R")
source("./R/est-fun.R")
source("./R/fit-density-by-hazard.R")
source("./simulation/coverage-fun.R")
source("./simulation/repeat-fun.R")
source("./simulation/compute-true.R")

#-------------------------------------------------------------------------------------------#
## parameters of interest
#-------------------------------------------------------------------------------------------#

K <- 100
run.ltmle <- FALSE
run.ctmle2 <- TRUE
misspecify.Q <- FALSE
misspecify.Q <- FALSE
only.A0 <- FALSE
M <- 500
n <- 1000

#-------------------------------------------------------------------------------------------#
## don't change here
#-------------------------------------------------------------------------------------------#

run.ctmle <- TRUE

#-------------------------------------------------------------------------------------------#
## get output from simulations
#-------------------------------------------------------------------------------------------#

length(out <- readRDS(file=paste0("./simulation/output/",
                                  "outlist-est",
                                  ifelse(run.ltmle, "-ltmle", ""),
                                  ifelse(run.ctmle, "-ctmle", ""),
                                  ifelse(run.ctmle2, "-ctmle2", ""),
                                  ifelse(n==1000, "", paste0("-n", n)), 
                                  "-2020", "-K", K, ifelse(only.A0, "-A0", ""),
                                  ifelse(misspecify.Q, "-Q", ""), 
                                  "-M", M, ".rds")))

#-------------------------------------------------------------------------------------------#
## get true values
#-------------------------------------------------------------------------------------------#

(psi0.A0 <- readRDS(file=paste0("./simulation/output/",
                                "outlist-est-true-0-2020",
                                "-K", K, 
                                ifelse(only.A0, "-A0", ""),
                                "-M", "1", ".rds")))

(psi0.A1 <- readRDS(file=paste0("./simulation/output/",
                                "outlist-est-true-1-2020",
                                "-K", K, 
                                ifelse(only.A0, "-A0", ""),
                                "-M", "1", ".rds")))

#-------------------------------------------------------------------------------------------#
## extract results of interest
#-------------------------------------------------------------------------------------------#

psiA0 <- unlist(lapply(out, function(xout, A=0, which=1) {
    if (xout[1]!="ERROR") {
        yout <- xout[[2-A]]
        if (run.ltmle) yout[1] else yout[[length(yout)]][which]
    }
}))

psiA0.init <- unlist(lapply(out, function(xout, A=0, which=1) {
    if (xout[1]!="ERROR") {
        yout <- xout[[2-A]]
        yout[[1]][which]
    }
}))

psiA1 <- unlist(lapply(out, function(xout, A=1, which=1) {
    if (xout[1]!="ERROR") {
        yout <- xout[[2-A]]
        if (run.ltmle) yout[1] else yout[[length(yout)]][which]
    }
}))

psiA1.init <- unlist(lapply(out, function(xout, A=1, which=1) {
    if (xout[1]!="ERROR") {
        yout <- xout[[2-A]]
        yout[[1]][which]
    }
}))

sdA0 <- unlist(lapply(out, function(xout, A=0, which=3) {
    if (xout[1]!="ERROR") {
        yout <- xout[[2-A]]
        if (run.ltmle) yout[2] else yout[[1]][which]
    }
}))

sdA1 <- unlist(lapply(out, function(xout, A=1, which=3) {
    if (xout[1]!="ERROR") {
        yout <- xout[[2-A]]
        if (run.ltmle) yout[2] else yout[[1]][which]
    }
}))

(count.na <- mean(unlist(lapply(out, function(xout, A=1, which=3) {
    1*(xout[1]=="ERROR")
}))))

if (FALSE) {

    hist(sdA1)
    hist(sdA0)

    par(mfrow=c(2,2))
    hist(psiA0.init)
    abline(v=psi0.A0, col="red")
    hist(psiA1.init)
    abline(v=psi0.A1, col="red")        
    hist(psiA0)
    abline(v=psi0.A0, col="red")
    hist(psiA1)
    abline(v=psi0.A1, col="red")        
    par(mfrow=c(1, 1))
    
}


#-------------------------------------------------------------------------------------------#
## extract results of interest
#-------------------------------------------------------------------------------------------#

message("----------------------------")
message("look at estimates from tmle:") 
message("------------------")
message(paste0("init (A=0): ", round(mean(psiA0.init), 4)))
message(paste0("tmle (A=0): ", round(mean(psiA0), 4)))
message(paste0("true (A=0): ", round(psi0.A0, 4)))
message("------------------")
message(paste0("init (A=1): ", round(mean(psiA1.init), 4)))
message(paste0("tmle (A=1): ", round(mean(psiA1), 4)))
message(paste0("true (A=1): ", round(psi0.A1, 4)))

message("-----------------------------------")
message("look at coverage of tmle estimator:")
message(paste0("coverage (A=0): ", round(cov.fun(psiA0, sdA0, psi0.A0), 4)))
message(paste0("coverage (A=1): ", round(cov.fun(psiA1, sdA1, psi0.A1), 4)))
message(paste0("coverage (diff): ", round(cov.fun(psiA1-psiA0, sqrt(sdA1^2+sdA0^2), psi0.A1-psi0.A0), 4)))

message("----------------------------")
message("look at se estimates (A=0):")
message(paste0("mean sigma : ", round(mean(sdA0), 4)))
message(paste0("mse        : ", round(mse(psiA0), 4)))

message("----------------------------")
message("look at se estimates (A=1):")
message(paste0("mean sigma : ", round(mean(sdA1), 4)))
message(paste0("mse        : ", round(mse(psiA1), 4)))

message("-------------------------------")
message("look at efficiency (A=0):")
message(paste0("MSE / sd       : ", round(mse(psiA0)/mean(sdA0), 4)))

message("-------------------------------")
message("look at efficiency (A=1):")
message(paste0("MSE / sd       : ", round(mse(psiA1)/mean(sdA1), 4)))


