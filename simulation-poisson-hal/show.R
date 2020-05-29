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

length(out <- readRDS(file=paste0("./simulation-poisson-hal/output/",
                                  "outlist-est",
                                  ifelse(interaction.AL, "-interactionAL", ""),
                                  "-M", M, ".rds")))

#-------------------------------------------------------------------------------------------#
## get true value
#-------------------------------------------------------------------------------------------#

(psi0 <- readRDS(file=paste0("./simulation-poisson-hal/output/",
                             "outlist-psi0",
                             ifelse(interaction.AL, "-interactionAL", ""),
                             ".rds")))

#-------------------------------------------------------------------------------------------#
## look at results of interest
#-------------------------------------------------------------------------------------------#

#-- initial estimator: 
mean(unlist(lapply(out, function(x) x[2])))

#-- tmle estimator: 
mean(unlist(lapply(out, function(x) x[length(x)])))
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


