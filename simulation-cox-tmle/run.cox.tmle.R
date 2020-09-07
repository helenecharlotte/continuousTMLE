### run.R --- 
#----------------------------------------------------------------------
## Author: Helene Charlotte Rytgaard
## Created: June, 2020 
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
library(Matrix)
library(coefplot)

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

source("./R/sim.data.continuous.R")
source("./R/cox.tmle.R")
source("./R/cox.super.learner.R")
source("./R/poisson.hal.R")
source("./R/poisson.hal.sl.R")

#-------------------------------------------------------------------------------------------#
## set parameters for simulations
#-------------------------------------------------------------------------------------------#

M <- 500
get.truth <- FALSE

betaA <- -0.15
betaL <- 1.1
nu <- 1.7
nu1 <- FALSE
eta <- 0.7
tau <- 1

interaction.AL <- FALSE
square.effect <- FALSE
square.effect2 <- TRUE
interaction.Atime <- FALSE

randomize.A <- TRUE
censoring.informative <- TRUE

misspecify.Y <- FALSE
fit.km <- TRUE
output.hr <- TRUE
SL <- FALSE
poisson.initial <- FALSE
lambda.cv <- NULL
penalize.time <- FALSE
adjust.penalization <- TRUE
SL.poisson <- FALSE

cens.misspecified <- FALSE
km.cens <- FALSE
SL.cens <- FALSE
poisson.cens <- FALSE

#-------------------------------------------------------------------------------------------#
## a bit of reprocessing
#-------------------------------------------------------------------------------------------#

if (interaction.Atime) betaA <- -0.7
if (nu1) nu <- 1

t0 <- 0.9
tau <- 1.2
change.point <- NULL

if (cens.misspecified) {
    cens.model <- Surv(time, delta==0)~L1
} else {
    cens.model <- Surv(time, delta==0)~L1+L2+L3+A*L1
}

if (interaction.Atime) {
    if (misspecify.Y) {
        change.point <- NULL
        outcome.model <- Surv(time, delta==1)~A*L3+L1+L2+L3*L1
        outcome.model <- Surv(time, delta==1)~A+L1+L2+L3
    } else {
        change.point <- t0
        outcome.model <- Surv(time, delta==1)~A+L1+L2+L3
    }
    if (interaction.AL) mod.period2 <- "*L3" else mod.period2 <- ""
} else if (interaction.AL) {
    mod.period2 <- ""
    if (misspecify.Y) {
        outcome.model <- Surv(time, delta==1)~A
    } else {
        outcome.model <- Surv(time, delta==1)~A*L1.squared+A+L1.squared+L1*L2+L3
    }
} else if (square.effect | square.effect2) {
    mod.period2 <- ""
    if (misspecify.Y) {
        outcome.model <- Surv(time, delta==1)~A+L1
    } else {
        outcome.model <- Surv(time, delta==1)~A+L1.squared
    }
} else {
    mod.period2 <- ""
    if (misspecify.Y) {
        outcome.model <- Surv(time, delta==1)~A+L1+L2+L3
    } else { 
        outcome.model <- Surv(time, delta==1)~A+L1.squared+L1*L2.root+L3.sin
    }
}

#-------------------------------------------------------------------------------------------#
## get the true value of target parameter(s)
#-------------------------------------------------------------------------------------------#

if (get.truth) {

    source("./R/sim.data.continuous.R")

    par(mfrow=c(1,2))
    
    print(psi0.A1 <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                              categorical=FALSE,
                              intervention.A=1, t0=t0, tau=tau,
                              square.effect=square.effect,
                              interaction.Atime=interaction.Atime, verbose=TRUE,
                              square.effect2=square.effect2,
                              interaction.AL=interaction.AL))
    print(psi0.A0 <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                              categorical=FALSE,
                              intervention.A=0, t0=t0, tau=tau,
                              square.effect=square.effect,
                              square.effect2=square.effect2,
                              interaction.Atime=interaction.Atime, verbose=TRUE,
                              interaction.AL=interaction.AL))
    print(psi0 <- psi0.A1 - psi0.A0) 

    saveRDS(c(psi0=psi0, psi0.A1=psi0.A1, psi0.A0=psi0.A0),
            file=paste0("./simulation-cox-tmle/output/",
                        "outlist-psi0-tau", tau, 
                        ifelse(interaction.AL, "-interactionAL", ""),
                        ifelse(interaction.Atime, "-interactionAtime", ""),
                        ifelse(square.effect, "-squareL2", ""),
                        ifelse(square.effect2, "-squareL2unif", ""),
                        ifelse(nu1, "-nu1", ""),
                        ".rds"))

}


#-------------------------------------------------------------------------------------------#
## repeat simulations (parallelize)
#-------------------------------------------------------------------------------------------#

if (system("echo $USER",intern=TRUE)%in%c("jhl781")){ 
    no_cores <- 12
} else {
    no_cores <- detectCores() - 1
}

registerDoParallel(no_cores)

out <- foreach(m=1:M, .errorhandling="pass"
               ) %dopar% {

                   dt <- sim.data(1000, betaA=betaA, betaL=betaL, nu=nu, eta=eta, t0=t0,
                                  seed=m+100,
                                  categorical=FALSE, randomize.A=randomize.A,
                                  censoring.informative=censoring.informative,
                                  interaction.Atime=interaction.Atime,
                                  square.effect=square.effect,
                                  square.effect2=square.effect2,
                                  interaction.AL=interaction.AL)
                   
                   print(paste0("m=", m))
                   
                   return(
                       list(#-- estimate difference effect:
                           "diff"=cox.tmle(dt, change.point=change.point, browse=FALSE,
                                           poisson.initial=poisson.initial,
                                           adjust.penalization=adjust.penalization,
                                           cut.time=12,
                                           cut.time.A=8,
                                           poisson.cens=poisson.cens,
                                           lambda.cv=lambda.cv,
                                           treat.effect="both",
                                           cut.covars=8, cut.L1.A=8,
                                           cut.L.interaction=6,
                                           tau=tau, km.cens=km.cens,
                                           cens.model=cens.model,
                                           SL.poisson=SL.poisson, SL=SL,
                                           SL.cens=SL.cens,
                                           output.km=fit.km, output.hr=output.hr,
                                           penalize.time=penalize.time, 
                                           outcome.model=outcome.model, mod.period2=mod.period2),
                           #-- estimate effect A=1: 
                           "A=1"=cox.tmle(dt, change.point=change.point, browse=FALSE,
                                          poisson.initial=poisson.initial,
                                          adjust.penalization=adjust.penalization,
                                          cut.time=12,
                                          cut.time.A=8,
                                          poisson.cens=poisson.cens,
                                          lambda.cv=lambda.cv,
                                          cut.covars=8, cut.L1.A=8,
                                          cut.L.interaction=6,
                                          tau=tau, km.cens=km.cens,
                                          cens.model=cens.model,
                                          SL.poisson=SL.poisson,
                                          SL.cens=SL.cens,
                                          output.km=fit.km, output.hr=output.hr,
                                          penalize.time=penalize.time, 
                                          outcome.model=outcome.model, mod.period2=mod.period2,
                                          SL=SL),
                           #-- estimate effect A=0: 
                           "A=0"=cox.tmle(dt, change.point=change.point, browse2=FALSE,
                                          poisson.initial=poisson.initial,
                                          adjust.penalization=adjust.penalization,
                                          cut.time=8, cut.time.A=8,
                                          poisson.cens=poisson.cens,
                                          lambda.cv=lambda.cv,
                                          cut.covars=6, cut.L1.A=8,
                                          cut.L.interaction=6,
                                          treat.effect="0",
                                          tau=tau, km.cens=km.cens,
                                          cens.model=cens.model,
                                          SL.poisson=SL.poisson,
                                          SL.cens=SL.cens,
                                          output.km=fit.km, output.hr=output.hr,
                                          penalize.time=penalize.time,
                                          outcome.model=outcome.model, mod.period2=mod.period2,
                                          SL=SL))
                   )

               }

stopImplicitCluster()


saveRDS(out,
        file=paste0("./simulation-cox-tmle/output/",
                    "outlist-run-cox-tmle-tau", tau,
                    ifelse(interaction.AL, "-interactionAL", ""),
                    ifelse(interaction.Atime, "-interactionAtime", ""),
                    ifelse(square.effect, "-squareL2", ""),
                    ifelse(square.effect2, "-squareL2unif", ""),
                    ifelse(nu1, "-nu1", ""),
                    ifelse(randomize.A, "-randomizeA", ""),
                    ifelse(censoring.informative, "-informativeCensoring"),
                    ifelse(misspecify.Y, "-misspecifyY", ""),
                    ifelse(length(lambda.cv)>0, "-no-penalization", ""),
                    ifelse(km.cens, "-KMcensoring", ""),
                    ifelse(cens.misspecified, "-censoringMisspecified", ""),
                    ifelse(SL, "-useSL", ""),
                    ifelse(SL.cens, "-useSLforCensoring", ""),
                    ifelse(poisson.initial, "-testPoisson", ""),
                    ifelse(SL.poisson, "-useSLwithPoisson", ""),
                    ifelse(poisson.initial & !adjust.penalization, "-noAdjustPenalization", ""),
                    ifelse(poisson.initial & !penalize.time, "-noPenalizationTime", ""),
                    ifelse(poisson.cens, "-testPoissonCensoring", ""),
                    "-M", M, ".rds"))


if (FALSE) {

    out <- readRDS(file=paste0("./simulation-cox-tmle/output/",
                               "outlist-run-cox-tmle-tau", tau,
                               ifelse(interaction.AL, "-interactionAL", ""),
                               ifelse(interaction.Atime, "-interactionAtime", ""),
                               ifelse(square.effect, "-squareL2", ""),
                               ifelse(square.effect2, "-squareL2unif", ""),
                               ifelse(nu1, "-nu1", ""),
                               ifelse(randomize.A, "-randomizeA", ""),
                               ifelse(censoring.informative, "-informativeCensoring"),
                               ifelse(misspecify.Y, "-misspecifyY", ""),
                               ifelse(length(lambda.cv)>0, "-no-penalization", ""),
                               ifelse(km.cens, "-KMcensoring", ""),
                               ifelse(cens.misspecified, "-censoringMisspecified", ""),
                               ifelse(SL, "-useSL", ""),
                               ifelse(SL.cens, "-useSLforCensoring", ""),
                               ifelse(poisson.initial, "-testPoisson", ""),
                               ifelse(SL.poisson, "-useSLwithPoisson", ""),
                               ifelse(poisson.initial & !adjust.penalization, "-noAdjustPenalization", ""),
                               ifelse(poisson.initial & !penalize.time, "-noPenalizationTime", ""),
                               ifelse(poisson.cens, "-testPoissonCensoring", ""),
                               "-M", M, ".rds"))

    (psi0 <- readRDS(file=paste0("./simulation-cox-tmle/output/",
                                 "outlist-psi0-tau", tau, 
                                 ifelse(interaction.AL, "-interactionAL", ""),
                                 ifelse(interaction.Atime, "-interactionAtime", ""),
                                 ifelse(square.effect, "-squareL2", ""),
                                 ifelse(square.effect2, "-squareL2unif", ""),
                                 ifelse(nu1, "-nu1", ""),
                                 ".rds")))

    psi0.A1 <- psi0["psi0.A1"]
    psi0.A0 <- psi0["psi0.A0"]
    psi0 <- psi0["psi0"]

    par(mfrow=c(2,4))
    hist(psiA1.est <- unlist(lapply(out, function(out1) {
        out2 <- out1[["A=1"]]
        if (length(out2)>1) return(out2[[1]][1])
    })), xlab="", ylab="", main="initial: A=1")
    abline(v=psi0.A1, col="red")
    hist(psiA0.est <- unlist(lapply(out, function(out1) {
        out2 <- out1[["A=0"]]
        if (length(out2)>1) return(out2[[1]][1])
    })), xlab="", ylab="", main="initial: A=0")
    abline(v=psi0.A0, col="red")
    hist(psiA1.est-psiA0.est, main="initial: diff")
    abline(v=psi0.A1-psi0.A0, col="red")
    hist(psi.est <- unlist(lapply(out, function(out1) {
        out2 <- out1[["diff"]]
        if (length(out2)>1) return(out2[[1]][1])
    })), xlab="", ylab="", main="initial: diff (directly targeted")
    abline(v=psi0.A1-psi0.A0, col="red")
    hist(psiA1.est <- unlist(lapply(out, function(out1) {
        out2 <- out1[["A=1"]]
        if (length(out2)>1) return(out2[[length(out2)]][1])
    })), xlab="", ylab="", main="tmle: A=1")
    abline(v=psi0.A1, col="red")
    hist(psiA0.est <- unlist(lapply(out, function(out1) {
        if (length(out1)>1) { 
            out2 <- out1[["A=0"]]
            out3 <- try(out2[[length(out2)]][1])
            if (!inherits(out3, "try-error")) return(out3)
        }
    })), xlab="", ylab="", main="tmle: A=0")
    abline(v=psi0.A0, col="red")
    hist(psiA1.est-psiA0.est, main="tmle: diff")
    abline(v=psi0.A1-psi0.A0, col="red")
    hist(psi.est <- unlist(lapply(out, function(out1) {
        if (length(out1)>1) { 
            out2 <- out1[["diff"]]
            out3 <- try(out2[[length(out2)]][1])
            if (!inherits(out3, "try-error")) return(out3)
        }
    })), xlab="", ylab="", main="tmle: diff (directly targeted)")
    abline(v=psi0.A1-psi0.A0, col="red")
    par(mfrow=c(1,1))

    (cov.A1 <- mean(unlist(lapply(out, function(out1, true=psi0.A1) {
        out2 <- out1[["A=1"]]
        if (length(out2)>1) {
            est <- out2[[length(out2)]][1]
            sd <- out2[[length(out2)]][2]
            return(est-1.96*sd<=true & true<=est+1.96*sd)
        }
    }))))
    (cov.A0 <- mean(unlist(lapply(out, function(out1, true=psi0.A0) {
        out2 <- out1[["A=0"]]
        if (length(out2)>1) {
            est <- try(out2[[length(out2)]][1])
            if (!inherits(est, "try-error")) {
                sd <- out2[[length(out2)]][2]
                return(est-1.96*sd<=true & true<=est+1.96*sd)
            }
        }
    }))))

    (bias.A1 <- mean(psiA1.est - psi0.A1))
    (bias.A0 <- mean(psiA0.est - psi0.A0))
    (bias.diff <- mean(psiA1.est-psiA0.est -
                       psi0.A1+psi0.A0))
    (bias.diff <- mean(psi.est -
                       psi0.A1+psi0.A0))

    (mse.A1 <- mse(psiA1.est))
    (mse.A0 <- mse(psiA0.est))
    (mse.diff <- mse(psiA1.est-psiA0.est))
    (mse.diff.target <- mse(psi.est))


}






