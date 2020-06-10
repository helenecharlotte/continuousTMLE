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
library(glmnet)
library(survival)
library(stringr)
library(ltmle)
library(nleqslv)
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
source("./R/cox-targeting.R")
source("./simulation-cox-targeting/repeat-fun.R")

#-------------------------------------------------------------------------------------------#
## set parameters
#-------------------------------------------------------------------------------------------#

betaA <- 0#0
betaL <- 0.3
nu <- 1.3#1.1#0.9#1.2
eta <- 1#1.2#1#4/sqrt(2)*(1/8)
tau <- 1#5
M <- 500#250#250
get.truth <- FALSE
interaction.AL <- TRUE
misspecify.Y <- TRUE

#-------------------------------------------------------------------------------------------#
## true value of target parameter
#-------------------------------------------------------------------------------------------#

if (get.truth) {
    
    print(psi0.A1 <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                              categorical=FALSE,
                              intervention.A=1, tau=tau,
                              interaction.AL=interaction.AL))
    print(psi0.A0 <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                              categorical=FALSE,
                              intervention.A=0, tau=tau,
                              interaction.AL=interaction.AL))

    print(psi0 <- psi0.A1 - psi0.A0)

    saveRDS(c(psi0=psi0, psi0.A1=psi0.A1, psi0.A0=psi0.A0),
            file=paste0("./simulation-cox-targeting/output/",
                        "outlist-psi0",
                        ifelse(interaction.AL, "-interactionAL", ""),
                        ".rds"))
}

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
## repeat simulations (parallelize)
#-------------------------------------------------------------------------------------------#

if (system("echo $USER",intern=TRUE)%in%c("jhl781")){ 
    no_cores <- 30
} else {
    no_cores <- detectCores() - 1
}

registerDoParallel(no_cores)

out <- foreach(m=1:M, .errorhandling="pass"#, #.combine=list, .multicombine = TRUE
               ) %dopar% {
                   repeat.fun(m, betaA=betaA, betaL=betaL, nu=nu, eta=eta, tau=tau,
                              misspecify.Y=misspecify.Y,
                              interaction.AL=interaction.AL, verbose=FALSE)
               }

stopImplicitCluster()


saveRDS(out,
        file=paste0("./simulation-cox-targeting/output/",
                    "outlist-est",
                    ifelse(interaction.AL, "-interactionAL", ""),
                    ifelse(misspecify.Y, "-misspecifyY", ""),
                    "-M", M, ".rds"))




