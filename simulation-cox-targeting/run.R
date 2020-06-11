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
library(survival)
library(riskRegression)

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
betaL <- 0.6
nu <- 1.7#1.1#0.9#1.2
eta <- 0.7#1#1.2#1#4/sqrt(2)*(1/8)
tau <- 1#5
M <- 500#500#100#250#250
get.truth <- FALSE
interaction.AL <- FALSE
misspecify.Y <- FALSE
interaction.Atime <- TRUE
fit.km <- TRUE

if (interaction.Atime) betaA <- -0.65
t0 <- 0.9
tau <- 1.2

#-------------------------------------------------------------------------------------------#
## true value of target parameter
#-------------------------------------------------------------------------------------------#

if (get.truth) {

    source("./R/sim-data-continuous.R")

    par(mfrow=c(1,2))
    print(psi0.A1 <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                              categorical=FALSE,
                              intervention.A=1, t0=t0, tau=tau, verbose=TRUE,
                              interaction.Atime=interaction.Atime,
                              interaction.AL=interaction.AL))
    print(psi0.A0 <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                              categorical=FALSE,
                              intervention.A=0, t0=t0, tau=tau,
                              interaction.Atime=interaction.Atime, verbose=TRUE,
                              interaction.AL=interaction.AL))
    print(psi0 <- psi0.A1 - psi0.A0)

    saveRDS(c(psi0=psi0, psi0.A1=psi0.A1, psi0.A0=psi0.A0),
            file=paste0("./simulation-cox-targeting/output/",
                        "outlist-psi0",
                        ifelse(interaction.AL, "-interactionAL", ""),
                        ifelse(interaction.Atime, "-interactionAtime", ""),
                        ".rds"))


    surv.A1.list <- list()
    tlist <- seq(0.5, 5, length=10)
    for (jj in 1:length(tlist)) {
        psi0.A1.jj <- sim.data(1e4, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                                  categorical=FALSE,
                                  intervention.A=1, t0=tlist[jj], tau=tau, verbose=TRUE,
                                  interaction.Atime=interaction.Atime,
                                  interaction.AL=interaction.AL)
        psi0.A0.jj <- sim.data(1e4, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                                  categorical=FALSE,
                                  intervention.A=0, t0=tlist[jj], tau=tau,
                                  interaction.Atime=interaction.Atime, verbose=TRUE,
                                  interaction.AL=interaction.AL)
        print(surv.A1.list[[jj]] <- psi0.A1.jj - psi0.A0.jj)
    }

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
                   repeat.fun(m, betaA=betaA, betaL=betaL, nu=nu, eta=eta, tau=tau, t0=t0,
                              misspecify.Y=misspecify.Y,
                              interaction.Atime=interaction.Atime, fit.km=fit.km,
                              interaction.AL=interaction.AL, verbose=TRUE, browse=FALSE)
               }

stopImplicitCluster()


saveRDS(out,
        file=paste0("./simulation-cox-targeting/output/",
                    "outlist-est",
                    ifelse(interaction.AL, "-interactionAL", ""),
                    ifelse(interaction.Atime, "-interactionAtime", ""),
                    ifelse(misspecify.Y, "-misspecifyY", ""),
                    "-M", M, ".rds"))




