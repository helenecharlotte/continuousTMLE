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

betaA <- 0.5#0
betaL <- 1.1#0.6#1.1
nu <- 1.7#1.1#0.9#1.2
eta <- 0.7#1#1.2#1#4/sqrt(2)*(1/8)
tau <- 1#5
M <- 1000#350#1000#250#1000#402#600#250#M <- 600#500#100#250#250
get.truth <- FALSE
interaction.AL <- FALSE#TRUE
misspecify.Y <- TRUE
interaction.Atime <- TRUE#TRUE#TRUE
randomize.A <- TRUE
censoring.informative <- TRUE#TRUE#TRUE#FALSE
censoring.high <- FALSE
fit.km <- TRUE
centered <- TRUE

#--  which effect estimated? 
a <- 0

if (interaction.Atime) betaA <- -0.75
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

    surv.list <- list()
    surv.A1.list <- list()
    surv.A0.list <- list()
    tlist <- seq(0, 1.5, length=19)
    for (jj in 1:length(tlist)) {
        psi0.A1.jj <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                               categorical=FALSE,
                               intervention.A=1, t0=t0, tau=tlist[jj], verbose=TRUE,
                               interaction.Atime=interaction.Atime,
                               interaction.AL=interaction.AL)
        psi0.A0.jj <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                               categorical=FALSE,
                               intervention.A=0, t0=t0, tau=tlist[jj],
                               interaction.Atime=interaction.Atime, verbose=TRUE,
                               interaction.AL=interaction.AL)
        print(surv.list[[jj]] <- psi0.A1.jj - psi0.A0.jj)
        surv.A1.list[[jj]] <- psi0.A1.jj
        surv.A0.list[[jj]] <- psi0.A0.jj
        saveRDS(cbind(t=tlist[1:jj], surv.list),
                file=paste0("./simulation-cox-targeting/output/",
                            "outlist-psi0-survival-function",
                            ifelse(interaction.AL, "-interactionAL", ""),
                            ifelse(interaction.Atime, "-interactionAtime", ""),
                            ".rds"))
        saveRDS(cbind(t=tlist[1:jj], surv.A1.list),
                file=paste0("./simulation-cox-targeting/output/",
                            "outlist-psi0-survival-A1-function",
                            ifelse(interaction.AL, "-interactionAL", ""),
                            ifelse(interaction.Atime, "-interactionAtime", ""),
                            ".rds"))
        saveRDS(cbind(t=tlist[1:jj], surv.A0.list),
                file=paste0("./simulation-cox-targeting/output/",
                            "outlist-psi0-survival-A0-function",
                            ifelse(interaction.AL, "-interactionAL", ""),
                            ifelse(interaction.Atime, "-interactionAtime", ""),
                            ".rds"))
    }

}

#-------------------------------------------------------------------------------------------#
## true intensity
#-------------------------------------------------------------------------------------------#

true.lambda <- function(t, A, L, betaA, betaL, nu, eta) {
    #return(nu*eta*t^{nu-1}*exp(A*betaA + L*betaL))
    return(exp(#-0.45+#0.55*A*(t<=tau/3)-0.65*A*(t>=tau/3)+
    (t<=t0)*betaA*A+
    (t>t0)*(-0.45)*betaA*A-
    L1*betaL-0.6*L2+0.8*L3-0.3*L3*L1
    + 0.3))
}

true.Lambda <- function(t, A, L, betaA, betaL, nu, eta) {
    return(eta*t^{nu}*exp(A*betaA + L*betaL))
}

#-------------------------------------------------------------------------------------------#
## repeat simulations (parallelize)
#-------------------------------------------------------------------------------------------#

if (system("echo $USER",intern=TRUE)%in%c("jhl781")){ 
    no_cores <- 35
} else {
    no_cores <- detectCores() - 1
}

registerDoParallel(no_cores)

out <- foreach(m=1:M, .errorhandling="pass"#, #.combine=list, .multicombine = TRUE
               ) %dopar% {
                   repeat.fun(m, betaA=betaA, betaL=betaL, nu=nu, eta=eta, tau=tau, t0=t0, a=a,
                              misspecify.Y=misspecify.Y, randomize.A=randomize.A,
                              censoring.informative=censoring.informative,
                              censoring.high=censoring.high, centered=TRUE, 
                              interaction.Atime=interaction.Atime, fit.km=fit.km,
                              interaction.AL=interaction.AL, verbose=TRUE, browse=FALSE)
               }

stopImplicitCluster()


saveRDS(out,
        file=paste0("./simulation-cox-targeting/output/",
                    "outlist-est",
                    ifelse(a==0, "-A=0", "-A=1"),
                    ifelse(interaction.AL, "-interactionAL", ""),
                    ifelse(interaction.Atime, "-interactionAtime", ""),
                    ifelse(randomize.A, "-randomizeA", ""),
                    ifelse(censoring.informative, "", "-almostUninformativeCensoring"),
                    ifelse(censoring.high, "-highLevelOfCensoring", ""),
                    ifelse(misspecify.Y, "-misspecifyY", ""),
                    ifelse(centered, "-centered", ""),
                    "-M", M, ".rds"))




