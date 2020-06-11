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

source("./R/sim-data-continuous.R")
source("./R/fit-poisson-hal.R")
source("./R/fit-coxph.R")
source("./R/estimate-target-from-intensity.R")
source("./R/log-linear-targeting.R")

#-------------------------------------------------------------------------------------------#
## parameters of interest
#-------------------------------------------------------------------------------------------#

M <- 250

#-------------------------------------------------------------------------------------------#
## get output from simulations
#-------------------------------------------------------------------------------------------#

length(out <- readRDS(file=paste0("./simulation-cox-targeting/output/",
                                  "outlist-est",
                                  ifelse(interaction.AL, "-interactionAL", ""),
                                  "-M", M, ".rds")))

#-------------------------------------------------------------------------------------------#
## get true value
#-------------------------------------------------------------------------------------------#

(psi0 <- readRDS(file=paste0("./simulation-cox-targeting/output/",
                             "outlist-psi0",
                             ifelse(interaction.AL, "-interactionAL", ""),
                             ifelse(interaction.Atime, "-interactionAtime", ""),
                             ".rds")))

psi0.A1 <- psi0["psi0.A1"]
psi0.A0 <- psi0["psi0.A0"]
psi0 <- psi0["psi0"]

#-------------------------------------------------------------------------------------------#
## look at results of interest
#-------------------------------------------------------------------------------------------#

M <- 450#500#500#200#250#100#500
misspecify.Y <- FALSE
interaction.AL <- FALSE
interaction.Atime <- TRUE
length(out <- readRDS(file=paste0("./simulation-cox-targeting/output/",
                                  "outlist-est",
                                  ifelse(interaction.AL, "-interactionAL", ""),
                                  ifelse(interaction.Atime, "-interactionAtime", ""),
                                  ifelse(misspecify.Y, "-misspecifyY", ""),
                                  "-M", M, ".rds")))

#-- initial and tmle1 estimator of psi0.A1: 
mean(unlist(lapply(out, function(x) x[[1]][1])))
mean(unlist(lapply(out, function(x) x[[length(x)]][1])))
psi0.A1

par(mfrow=c(1,2))
hist(unlist(lapply(out, function(x) x[[1]][1])), main="init", xlab="psi.hat", ylab="")
abline(v=psi0.A1, col="red")
hist(unlist(lapply(out, function(x) x[[length(x)]][1])), main="tmle (max. 5 iterations)", xlab="psi.hat", ylab="")
abline(v=psi0.A1, col="red")
par(mfrow=c(1,1))

#mean(unlist(lapply(out, function(x) x[[2]][2])))
mean(unlist(lapply(out, function(x) x[[length(x)]][2])))

mse(unlist(lapply(out, function(x) x[[1]][1])))
mse(unlist(lapply(out, function(x) x[[length(x)]][1])))
sd2 <- sd(unlist(lapply(out, function(x) x[[length(x)]][1])))

#-- check coverage: 
(coverage <- mean(unlist(lapply(out, function(x, psi0=psi0.A1) {
    c(x[[length(x)]][1]-1.96*x[[length(x)]][2]<=psi0 &
      x[[length(x)]][1]+1.96*x[[length(x)]][2]>=psi0)
}))))

#-- check "oracle" coverage: 
(coverage2 <- mean(unlist(lapply(out, function(x, psi0=psi0.A1) {
    c(x[[length(x)]][1]-1.96*sd2<=psi0 &
      x[[length(x)]][1]+1.96*sd2>=psi0)
}))))


#--- look at other;

#-- naive hr:
mean(unlist(lapply(out, function(x) x[[1]]["hr.A"])))
mean(unlist(lapply(out, function(x) x[[1]]["hr.pval"]))<=0.05)

#-- unadjusted km:
c("km.fit"=mean(unlist(lapply(out, function(x) x[[1]]["km.est"]))), "truth"=psi0.A1)
c("km.se"=sd(unlist(lapply(out, function(x) x[[1]]["km.est"]))), "tmle.se"=sd2)
c("km.se.est"=mean(unlist(lapply(out, function(x) x[[1]]["km.se"]))), "tmle.se.est"=mean(unlist(lapply(out, function(x) x[[length(x)]][2]))))

#-- model check; estimation of time-varying coef:
if (!misspecify.Y) {
    c(mean(unlist(lapply(out, function(x) x[[1]][5]))), betaA)
    c(mean(unlist(lapply(out, function(x) x[[1]][6]))), -0.35*betaA)
    c(mean(unlist(lapply(out, function(x) x[[1]][7]))), -betaL)
}




#------ old; 

#-- initial and tmle1 estimator of diff: 
mean(unlist(lapply(out, function(x) x$init["init.diff"])))
mean(unlist(lapply(out, function(x) x$tmle1["tmle.diff"])))
psi0

#-- initial and tmle1 estimator of psi0.A1: 
mean(unlist(lapply(out, function(x) x$init["init.A1"])))
mean(unlist(lapply(out, function(x) x$tmle1["tmle.A1"])))
psi0.A1

#-- initial and tmle1 estimator of psi0.A0: 
mean(unlist(lapply(out, function(x) x$init["init.A0"])))
mean(unlist(lapply(out, function(x) x$tmle1["tmle.A0"])))
psi0.A0







#-- tmle estimator: 
mean(na.omit(unlist(lapply(out, function(x) x[length(x)]))))
#mse(unlist(lapply(out, function(x) x[length(x)])))

#-- cox estimator: 
mean(unlist(lapply(out, function(x) x[1])))
#mse(unlist(lapply(out, function(x) x[1])))

par(mfrow=c(1,3))
hist(unlist(lapply(out, function(x) x[2])), main="initial")
abline(v=psi0, col="red")
hist(unlist(lapply(out, function(x) x[length(x)])), main="TMLE")
abline(v=psi0, col="red")
hist(unlist(lapply(out, function(x) x[1])), main="Cox")
abline(v=psi0, col="red")
par(mfrow=c(1,1))


