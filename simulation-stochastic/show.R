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

K <- 15#15#10#5#5#100#100#100#100#80#100
run.ltmle <- FALSE##TRUE#FALSE
run.ctmle <- TRUE#FALSE#FALSE
run.ctmle2 <- FALSE#FALSE#FALSE
compute.true.eic <- FALSE
compute.true.psi <- FALSE
misspecify.Q <- TRUE#FALSE#TRUE
M <- 600#250#600#250#100#250#350#600#551#300#300#551
n <- 1000
rescue.itt <- FALSE#TRUE
estimate.Gstar <- FALSE#TRUE#FALSE

#-------------------------------------------------------------------------------------------#
## get output from simulations
#-------------------------------------------------------------------------------------------#

length(out <- readRDS(file=paste0("./simulation-stochastic/output/",
                                  "outlist-est-rescue", ifelse(rescue.itt, "-itt", ""),
                                  ifelse(run.ltmle, "-ltmle", ""),
                                  ifelse(run.ctmle, "-ctmle", ""),
                                  ifelse(run.ctmle2, "-ctmle2", ""),
                                  ifelse(estimate.Gstar, "-estimateG", ""),
                                  ifelse(n==1000, "", paste0("-n", n)), 
                                  "-2020", "-K", K, ifelse(misspecify.Q, "-Q", ""), 
                                  "-M", M, ".rds")))

#-------------------------------------------------------------------------------------------#
## get true values
#-------------------------------------------------------------------------------------------#

(psi0.A0 <- readRDS(file=paste0("./simulation-stochastic/output/",
                                "outlist-est-true-rescue", ifelse(rescue.itt, "-itt", ""), 
                                "-0-2020",
                                "-K", K, #ifelse(misspecify.Q, "-Q", ""),
                                "-M", "2", ".rds")))

(psi0.A1 <- readRDS(file=paste0("./simulation-stochastic/output/",
                                "outlist-est-true-rescue", ifelse(rescue.itt, "-itt", ""), 
                                "-1-2020",
                                "-K", K, #ifelse(misspecify.Q, "-Q", ""),
                                "-M", "2", ".rds")))

psi0.A1 - psi0.A0

if (FALSE) { #-- FIXME: not adapted to stochastic/rescue yet.
    (true.eic.A0 <- readRDS(file=paste0("./simulation-stochastic/output/",
                                        "outlist-est-true-sd-0-2020",
                                        "-K", K, #ifelse(misspecify.Q, "-Q", ""),
                                        "-M", "1", ".rds"))*sqrt(ifelse(n==1000, 1, 1000/n)))

    (true.eic.A1 <- readRDS(file=paste0("./simulation-stochastic/output/",
                                        "outlist-est-true-sd-1-2020",
                                        "-K", K, #ifelse(misspecify.Q, "-Q", ""),
                                        "-M", "1", ".rds"))*sqrt(ifelse(n==1000, 1, 1000/n)))
}
    
#-------------------------------------------------------------------------------------------#
## extract results of interest
#-------------------------------------------------------------------------------------------#

mean(psiA0 <- unlist(lapply(out, function(xout, A=0, which=1) {
    if (xout[1]!="ERROR") {
        yout <- xout[[2-A]]
        if (run.ltmle) yout[1] else yout[[length(yout)]][which]
    }
})))

mean(psiA0.init <- unlist(lapply(out, function(xout, A=0, which=1) {
    if (xout[1]!="ERROR") {
        yout <- xout[[2-A]]
        yout[[1]][which]
    }
})))

mean(psiA1 <- unlist(lapply(out, function(xout, A=1, which=1) {
    if (xout[1]!="ERROR") {
        yout <- xout[[2-A]]
        if (run.ltmle) yout[1] else yout[[length(yout)]][which]
    }
})))

mean(psiA1.init <- unlist(lapply(out, function(xout, A=1, which=1) {
    if (xout[1]!="ERROR") {
        yout <- xout[[2-A]]
        yout[[1]][which]
    }
})))

par(mfrow=c(1,2))
hist(sdA0 <- unlist(lapply(out, function(xout, A=0, which=3) {
    if (xout[1]!="ERROR") {
        yout <- xout[[2-A]]
        if (run.ltmle) yout[2] else yout[[1]][which]
    }
})))

hist(sdA1 <- unlist(lapply(out, function(xout, A=1, which=3) {
    if (xout[1]!="ERROR") {
        yout <- xout[[2-A]]
        if (run.ltmle) yout[2] else yout[[1]][which]
    }
})))
par(mfrow=c(1,1))

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


if (FALSE) {
    mean(psiA1.init <- unlist(lapply(out, function(xout, A=1, which=1) {
        browser()
        if (xout[1]!="ERROR") {
            yout <- xout[[2-A]]
            yout[[1]][which]
        }
    })))
}


#-------------------------------------------------------------------------------------------#
## extract results of interest
#-------------------------------------------------------------------------------------------#

message("----------------------------")
message("look at estimates from tmle:") 
message("------------------")
message(paste0("init (A=0): ", round(mean(psiA0.init), 4), ", bias: ", round(mean(psiA0.init-psi0.A0), 4)))
message(paste0("tmle (A=0): ", round(mean(psiA0), 4), ", bias: ", round(mean(psiA0-psi0.A0), 4)))
message(paste0("true (A=0): ", round(psi0.A0, 4)))
message("------------------")
message(paste0("init (A=1): ", round(mean(psiA1.init), 4), ", bias: ", round(mean(psiA1.init-psi0.A1), 4)))
message(paste0("tmle (A=1): ", round(mean(psiA1), 4), ", bias: ", round(mean(psiA1-psi0.A1), 4)))
message(paste0("true (A=1): ", round(psi0.A1, 4)))
message("------------------")
message(paste0("init (diff): ", round(mean(psiA1.init-psiA0.init), 4), ", bias: ", round(mean(psiA1.init-psi0.A1-(psiA0.init-psi0.A0)), 4)))
message(paste0("tmle (diff): ", round(mean(psiA1-psiA0), 4), ", bias: ", round(mean(psiA1-psi0.A1-(psiA0-psi0.A0)), 4)))
message(paste0("true (diff): ", round(psi0.A1-psi0.A0, 4)))


message("-----------------------------------")
message("look at coverage of tmle estimator:")
message(paste0("coverage (A=0): ", round(cov.fun(psiA0, sdA0, psi0.A0), 4)))
message(paste0("coverage (A=1): ", round(cov.fun(psiA1, sdA1, psi0.A1), 4)))
message(paste0("coverage (diff): ", round(cov.fun(psiA1-psiA0, sqrt(sdA1^2+sdA0^2), psi0.A1-psi0.A0), 4)))


message(paste0("coverage (A=0): ", round(cov.fun(psiA0, mean(sdA0), psi0.A0), 4)))
message(paste0("coverage (A=1): ", round(cov.fun(psiA1, mean(sdA1), psi0.A1), 4)))


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
#message(paste0("MSE / true eic : ", round(mse(psiA0)/true.eic.A0, 4)))

message("-------------------------------")
message("look at efficiency (A=1):")
message(paste0("MSE / sd       : ", round(mse(psiA1)/mean(sdA1), 4)))
#message(paste0("MSE / true eic : ", round(mse(psiA1)/true.eic.A1, 4)))



#-------------------------------------------------------------------------------------------#
## check convergence and estimated eps
#-------------------------------------------------------------------------------------------#

if (!run.ltmle) {
    no.iterations.A0 <- unlist(lapply(out, function(xout, A=0) {
        yout <- xout[[2-A]]
        yout[[1]]["no.iteration"]
    }))
    no.iterations.A1 <- unlist(lapply(out, function(xout, A=1) {
        yout <- xout[[2-A]]
        yout[[1]]["no.iteration"]
    }))
    table(no.iterations.A0)
    table(no.iterations.A1)

    eps.hat.A0 <- unlist(lapply(out, function(xout, A=0) {
        yout <- xout[[2-A]]
        yout[[length(yout)]]["eps.hat.H"]
    }))

    eps.hat.A1 <- unlist(lapply(out, function(xout, A=1) {
        yout <- xout[[2-A]]
        yout[[length(yout)]]["eps.hat.H"]
    }))
    
    cbind(no.iterations.A1, eps.hat.A1, psiA1.init, psiA1)
    cbind(no.iterations.A0, eps.hat.A0, psiA0.init, psiA0)

    out[[4]][[2]]
}
