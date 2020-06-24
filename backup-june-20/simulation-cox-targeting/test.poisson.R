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
source("./R/fit-poisson-hal2.R")
source("./R/fit-coxph.R")
source("./R/estimate-target-from-intensity.R")
source("./R/log-linear-targeting.R")
source("./R/cox-targeting.R")
source("./R/cox.tmle.R")
source("./simulation-cox-targeting/repeat-fun.R")

#-------------------------------------------------------------------------------------------#
## set parameters
#-------------------------------------------------------------------------------------------#

betaA <- -0.15#0
betaL <- 1.1#0.6#1.1
nu <- 1.7#1.1#0.9#1.2
eta <- 0.7#1#1.2#1#4/sqrt(2)*(1/8)
tau <- 1#5
M <- 1000#1000#1000#410#320#1000#350#1000#250#1000#402#600#250#M <- 600#500#100#250#250
get.truth <- FALSE
interaction.AL <- FALSE#FALSE#TRUE
misspecify.Y <- TRUE#FALSE
interaction.Atime <- TRUE#TRUE#TRUE#TRUE#TRUE
randomize.A <- TRUE
censoring.informative <- TRUE#TRUE#FALSE
censoring.high <- FALSE
fit.km <- TRUE
centered <- TRUE
test.poisson <- TRUE#TRUE

#--  which effect estimated? 
a <- 1

if (interaction.Atime) betaA <- -0.75
t0 <- 0.9
tau <- 1.2
change.point <- NULL

if (interaction.Atime) {
    if (misspecify.Y) {
        change.point <- NULL
        outcome.model <- Surv(time, delta==1)~A*L3+L1+L2+L3*L1
    } else {
        change.point <- t0
        outcome.model <- Surv(time, delta==1)~A+L1+L2+L3
    }
    if (interaction.AL) mod.period2 <- "*L3" else mod.period2 <- ""
} else {
    mod.period2 <- ""
    if (misspecify.Y) {
        outcome.model <- Surv(time, delta==1)~A
    } else {
        outcome.model <- Surv(time, delta==1)~A*L1.squared+A+L1.squared+L1*L2+L3
    }
}

#-------------------------------------------------------------------------------------------#
## repeat simulations (parallelize)
#-------------------------------------------------------------------------------------------#

if (system("echo $USER",intern=TRUE)%in%c("jhl781")){ 
    no_cores <- 25
} else {
    no_cores <- detectCores() - 1
}

registerDoParallel(no_cores)

out <- foreach(m=1:M, .errorhandling="pass"#, #.combine=list, .multicombine = TRUE
               ) %dopar% {

                   dt <- sim.data(1000, betaA=betaA, betaL=betaL, nu=nu, eta=eta, t0=t0,
                                  seed=m+100,
                                  categorical=FALSE, randomize.A=randomize.A,
                                  censoring.informative=censoring.informative,
                                  censoring.high=censoring.high,
                                  interaction.Atime=interaction.Atime,
                                  interaction.AL=interaction.AL)

                   print(paste0("m=", m))

                   list(#-- estimate difference effect:
                       #"diff"=cox.tmle(dt, treat.effect="both", change.point=change.point,
                       #                outcome.model=outcome.model, mod.period2=mod.period2),
                       #-- estimate effect A=1: 
                       "A=1"=cox.tmle(dt, change.point=change.point,browse=FALSE,
                                      poisson.initial=TRUE, cut.time=12, cut.time.A=12,
                                      lambda.cv=NULL,
                                      cut.covars=6,
                                      outcome.model=outcome.model, mod.period2=mod.period2),
                       #-- estimate effect A=0: 
                       "A=0"=cox.tmle(dt, change.point=change.point,browse=FALSE,
                                      poisson.initial=TRUE, cut.time=12, cut.time.A=12,
                                      #lambda.cv=0.001,
                                      cut.covars=6,
                                      treat.effect="0",
                                      outcome.model=outcome.model, mod.period2=mod.period2))

               }

stopImplicitCluster()


saveRDS(out,
        file=paste0("./simulation-cox-targeting/output/",
                    "outlist-cox-tmle-est-use-poisson",
                    #ifelse(a==0, "-A=0", "-A=1"),
                    ifelse(interaction.AL, "-interactionAL", ""),
                    ifelse(interaction.Atime, "-interactionAtime", ""),
                    ifelse(randomize.A, "-randomizeA", ""),
                    ifelse(censoring.informative, "", "-almostUninformativeCensoring"),
                    ifelse(censoring.high, "-highLevelOfCensoring", ""),
                    #ifelse(misspecify.Y, "-misspecifyY", ""),
                    ifelse(centered, "-centered", ""),
                    ifelse(test.poisson, "-testPoisson", ""),
                    "-M", M, ".rds"))


if (FALSE) {

    cox.tmle(dt, change.point=change.point, #browse=TRUE,
             poisson.initial=TRUE, cut.time=5, cut.time.A=5,
             #lambda.cv=0.001,
             cut.covars=5,
             treat.effect="1",
             outcome.model=outcome.model, mod.period2=mod.period2)
    cox.tmle(dt, change.point=change.point, #browse=TRUE,
             poisson.initial=TRUE, cut.time=5, cut.time.A=5,
             #lambda.cv=0.001,
             cut.covars=5,
             treat.effect="0",
             outcome.model=outcome.model, mod.period2=mod.period2)
    
    M <- 500#35#150
    test.poisson <- TRUE
    out <- readRDS(file=paste0("./simulation-cox-targeting/output/",
                               "outlist-cox-tmle-est-use-poisson",
                               #ifelse(a==0, "-A=0", "-A=1"),
                               ifelse(interaction.AL, "-interactionAL", ""),
                               ifelse(interaction.Atime, "-interactionAtime", ""),
                               ifelse(randomize.A, "-randomizeA", ""),
                               ifelse(censoring.informative, "", "-almostUninformativeCensoring"),
                               ifelse(censoring.high, "-highLevelOfCensoring", ""),
                               #ifelse(misspecify.Y, "-misspecifyY", ""),
                               ifelse(centered, "-centered", ""),
                               ifelse(test.poisson, "-testPoisson", ""),
                               "-M", M, ".rds"))

     (psi0 <- readRDS(file=paste0("./simulation-cox-targeting/output/",
                             "outlist-psi0",
                             ifelse(interaction.AL, "-interactionAL", ""),
                             ifelse(interaction.Atime, "-interactionAtime", ""),
                             ".rds")))

    psi0.A1 <- psi0["psi0.A1"]
    psi0.A0 <- psi0["psi0.A0"]
    psi0 <- psi0["psi0"]

    par(mfrow=c(1,2))
    #hist(diff.est <- unlist(lapply(out, function(out1) {
    #    out2 <- out1[[1]]
    #    out2[[length(out2)]][1]
    #})), xlab="", ylab="", main="estimation: diff")
    #abline(v=psi0, col="red")
    hist(psiA1.est <- unlist(lapply(out, function(out1) {
        out2 <- out1[[1]]
        out2[[length(out2)]][1]
    })), xlab="", ylab="", main="estimation: A=1")
    abline(v=psi0.A1, col="red")
    hist(psiA0.est <- unlist(lapply(out, function(out1) {
        out2 <- out1[[2]]
        out2[[length(out2)]][1]
    })), xlab="", ylab="", main="estimation: A=0")
    abline(v=psi0.A0, col="red")
    par(mfrow=c(1,1))

    (cov.A1 <- mean(unlist(lapply(out, function(out1, true=psi0.A1) {
        out2 <- out1[[1]]
        est <- out2[[length(out2)]][1]
        sd <- out2[[length(out2)]][2]
        return(est-1.96*sd<=true & true<=est+1.96*sd)
    }))))
    (cov.A0 <- mean(unlist(lapply(out, function(out1, true=psi0.A0) {
        out2 <- out1[[2]]
        est <- out2[[length(out2)]][1]
        sd <- out2[[length(out2)]][2]
        return(est-1.96*sd<=true & true<=est+1.96*sd)
    }))))

    (mse.A1 <- mse(unlist(lapply(out, function(out1) {
        out2 <- out1[[1]]
        out2[[length(out2)]][1]
    }))))
    (mse.A0 <- mse(unlist(lapply(out, function(out1) {
        out2 <- out1[[2]]
        out2[[length(out2)]][1]
    }))))


}








if (FALSE) {

    m <- 2#m+1
    dt <- sim.data(200, betaA=betaA, betaL=betaL, nu=nu, eta=eta, t0=t0,
                   seed=m+100,
                   categorical=FALSE, randomize.A=randomize.A,
                   censoring.informative=censoring.informative,
                   censoring.high=censoring.high,
                   interaction.Atime=interaction.Atime,
                   interaction.AL=interaction.AL)

    source("./R/cox.tmle.R")
    (out1 <- cox.tmle(dt, change.point=change.point,browse=FALSE, one.step=TRUE, deps.size=0.01,
                      estimate.curve=TRUE, no.small.steps=250, browse2=TRUE,
                      outcome.model=outcome.model, mod.period2=mod.period2))
    (out0 <- cox.tmle(dt, change.point=change.point,browse=FALSE, one.step=TRUE, deps.size=0.01,
                      estimate.curve=TRUE, no.small.steps=250, treat.effect="0",
                      outcome.model=outcome.model, mod.period2=mod.period2))

    tau2 <- tau
    out1[tau==max(tau[tau<=tau2])]
    out0[tau==max(tau[tau<=tau2])]

    surv.list.A1 <- cbind(readRDS(file=paste0("./simulation-cox-targeting/output/",
                                              "outlist-psi0-survival-A1-function",
                                              ifelse(interaction.AL, "-interactionAL", ""),
                                              ifelse(interaction.Atime, "-interactionAtime", ""),
                                              ".rds")), Treatment="a=1")
    surv.list.A0 <- cbind(readRDS(file=paste0("./simulation-cox-targeting/output/",
                                              "outlist-psi0-survival-A0-function",
                                              ifelse(interaction.AL, "-interactionAL", ""),
                                              ifelse(interaction.Atime, "-interactionAtime", ""),
                                              ".rds")), Treatment="a=1")
    
    par(mfrow=c(1,2))
    plot(out1[tau<=max(unlist(surv.list.A1[, 1])), tau],
         out1[tau<=max(unlist(surv.list.A1[, 1])), tmle.fit], type="l")
    lines(unlist(surv.list.A1[, 1]), unlist(surv.list.A1[, 2]), col="blue")
    plot(out0[tau<=max(unlist(surv.list.A0[, 1])), tau],
         out0[tau<=max(unlist(surv.list.A0[, 1])), tmle.fit], type="l")
    lines(unlist(surv.list.A0[, 1]), unlist(surv.list.A0[, 2]), col="blue")
    par(mfrow=c(1,1))

    
    #-- only one time-point

    cox.tmle(dt, change.point=change.point,browse=FALSE, one.step=TRUE, deps.size=0.01,
             outcome.model=outcome.model, mod.period2=mod.period2)

    cox.tmle(dt, change.point=change.point,browse=FALSE, one.step=FALSE,
             outcome.model=outcome.model, mod.period2=mod.period2)

    cox.tmle(dt, change.point=change.point,browse=FALSE, one.step=TRUE, deps.size=0.01,
             treat.effect="0",
             outcome.model=outcome.model, mod.period2=mod.period2)

    cox.tmle(dt, change.point=change.point,browse=FALSE, one.step=FALSE, treat.effect="0",
             outcome.model=outcome.model, mod.period2=mod.period2)
    
}

if (FALSE) {

    M <- 1000#1000#300#50#1000
    misspecify.Y <- FALSE#FALSE#TRUE
    out <- readRDS(file=paste0("./simulation-cox-targeting/output/",
                               "outlist-cox-tmle-est",
                               #ifelse(a==0, "-A=0", "-A=1"),
                               ifelse(interaction.AL, "-interactionAL", ""),
                               ifelse(interaction.Atime, "-interactionAtime", ""),
                               ifelse(randomize.A, "-randomizeA", ""),
                               ifelse(censoring.informative, "", "-almostUninformativeCensoring"),
                               ifelse(censoring.high, "-highLevelOfCensoring", ""),
                               ifelse(misspecify.Y, "-misspecifyY", ""),
                               ifelse(centered, "-centered", ""),
                               ifelse(test.poisson, "-testPoisson", ""),
                               "-M", M, ".rds"))

    (psi0 <- readRDS(file=paste0("./simulation-cox-targeting/output/",
                             "outlist-psi0",
                             ifelse(interaction.AL, "-interactionAL", ""),
                             ifelse(interaction.Atime, "-interactionAtime", ""),
                             ".rds")))

    psi0.A1 <- psi0["psi0.A1"]
    psi0.A0 <- psi0["psi0.A0"]
    psi0 <- psi0["psi0"]

    par(mfrow=c(1,3))
    hist(diff.est <- unlist(lapply(out, function(out1) {
        out2 <- out1[[1]]
        out2[[length(out2)]][1]
    })), xlab="", ylab="", main="estimation: diff")
    abline(v=psi0, col="red")
    hist(psiA1.est <- unlist(lapply(out, function(out1) {
        out2 <- out1[[2]]
        out2[[length(out2)]][1]
    })), xlab="", ylab="", main="estimation: A=1")
    abline(v=psi0.A1, col="red")
    hist(psiA0.est <- unlist(lapply(out, function(out1) {
        out2 <- out1[[3]]
        out2[[length(out2)]][1]
    })), xlab="", ylab="", main="estimation: A=0")
    abline(v=psi0.A0, col="red")
    par(mfrow=c(1,1))

    (cov.diff <- mean(unlist(lapply(out, function(out1, true=psi0) {
        out2 <- out1[[1]]
        est <- out2[[length(out2)]][1]
        sd <- out2[[length(out2)]][2]
        return(est-1.96*sd<=true & true<=est+1.96*sd)
    }))))
    (cov.A1 <- mean(unlist(lapply(out, function(out1, true=psi0.A1) {
        out2 <- out1[[2]]
        est <- out2[[length(out2)]][1]
        sd <- out2[[length(out2)]][2]
        return(est-1.96*sd<=true & true<=est+1.96*sd)
    }))))
    (cov.A0 <- mean(unlist(lapply(out, function(out1, true=psi0.A0) {
        out2 <- out1[[3]]
        est <- out2[[length(out2)]][1]
        sd <- out2[[length(out2)]][2]
        return(est-1.96*sd<=true & true<=est+1.96*sd)
    }))))

    (mse.diff <- mse(unlist(lapply(out, function(out1) {
        out2 <- out1[[1]]
        out2[[length(out2)]][1]
    }))))
}
