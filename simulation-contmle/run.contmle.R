### run.contmle.R --- 
#----------------------------------------------------------------------
## Author: Helene Charlotte Wiese Rytgaard
## Created: October, 2020 
#----------------------------------------------------------------------
## 
### Commentary:
## Run continuous tmle across different sample sizes and different
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
source("./R/contmle.R")
source("./R/cox.super.learner.R")
source("./R/poisson.hal.R")
source("./R/poisson.hal.sl.R")

#-------------------------------------------------------------------------------------------#
## set parameters for simulations
#-------------------------------------------------------------------------------------------#

#M <- 503; 504; 505; 507
M         <- 500#1000#1001#600#5#501#1000#10#500#500#500#500#500#508#100#506#25#503#501#25#100#5#100#100#500#500
n         <- 1000
get.truth <- FALSE#FALSE

betaA <- -0.15
betaL <- 1.1
nu    <- 1.7
eta   <- 0.7
t0    <- 0.9

treat.effect <- "both"
tau          <- 1.2

square.effect1        <- TRUE#FALSE#TRUE
square.effect2        <- FALSE#TRUE#FALSE
interaction.Atime     <- TRUE#FALSE#TRUE
reversed.setting      <- TRUE#FALSE#TRUE
competing.risk        <- FALSE#TRUE#FALSE
randomize.A           <- TRUE
censoring.informative <- TRUE
censoring.high        <- FALSE
no.censoring          <- FALSE

misspecify.outcome <- FALSE
misspecify.cens    <- FALSE

fit.outcome  <- "cox"
fit.cens     <- "cox"
fit.cr       <- "cox"
one.step     <- FALSE


#-------------------------------------------------------------------------------------------#
## a bit of reprocessing
#-------------------------------------------------------------------------------------------#

if (interaction.Atime) betaA <- -0.7
change.point <- NULL
if (competing.risk) tau <- 0.7

if (reversed.setting) {
    betaA <- 0.5
    t0 <- 0.7
}

if (square.effect2) {
    if (misspecify.outcome) {
        outcome.model <- Surv(time, delta==1)~A+L1
    } else {
        outcome.model <- Surv(time, delta==1)~A+L1.squared
    }
}

if (interaction.Atime) {
    if (misspecify.outcome) {
        change.point <- NULL
        if ((square.effect2 | square.effect1) & reversed.setting) {
            outcome.model <- Surv(time, delta==1)~A+L1
        } else {
            outcome.model <- Surv(time, delta==1)~A+L1+L2+L3
        }
    } else {
        change.point <- t0
        if (square.effect2 | square.effect1) {
            if (reversed.setting) {
                outcome.model <- Surv(time, delta==1)~A+L1.squared+L2+L3
            } else {
                outcome.model <- Surv(time, delta==1)~A+L1+L2.squared+L3
            }
        } else {
            outcome.model <- Surv(time, delta==1)~A+L1+L2+L3
        }
    }
}

if (misspecify.cens) {
    cens.model <- Surv(time, delta==0)~A+L1+L2+L3
} else {
    if (censoring.informative) {
        if (square.effect1) {
            cens.model <- Surv(time, delta==0)~L1+L2+L3+A*L1.squared
        } else {
            cens.model <- Surv(time, delta==0)~L1+L2+L3+A*L1
        }
    } else {
        cens.model <- Surv(time, delta==0)~1
    }
}

#-------------------------------------------------------------------------------------------#
## get the true value of target parameter(s)
#-------------------------------------------------------------------------------------------#

if (get.truth) {

    source("./R/sim.data.continuous.R")

    par(mfrow=c(1,2))

    print(psi0.A1 <- sim.data(1e7, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                              categorical=FALSE, competing.risk=competing.risk,
                              intervention.A=1, tau=tau, t0=t0,
                              verbose=TRUE,
                              reversed.setting=reversed.setting,
                              square.effect2=square.effect2,
                              square.effect1=square.effect1,
                              interaction.Atime=interaction.Atime))
    print(psi0.A0 <- sim.data(1e7, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                              categorical=FALSE, competing.risk=competing.risk,
                              intervention.A=0, tau=tau, t0=t0,
                              verbose=TRUE,
                              reversed.setting=reversed.setting,
                              square.effect2=square.effect2,
                              square.effect1=square.effect1,
                              interaction.Atime=interaction.Atime))
    print(psi0 <- psi0.A1 - psi0.A0) 

    saveRDS(c(psi0=psi0, psi0.A1=psi0.A1, psi0.A0=psi0.A0),
            file=paste0("./simulation-contmle/output/",
                        "save-psi0",
                        ifelse(competing.risk, "-competingrisk", ""), 
                        "-tau", tau, 
                        ifelse(square.effect2, "-squareL2unif", ""),
                        ifelse(square.effect1, "-squareL1unif", ""),
                        ifelse(interaction.Atime, "-interactionAtime", ""),
                        ifelse(reversed.setting, "-reversedsetting", ""),
                        ".rds"))

    
    surv.list <- list()
    surv.A1.list <- list()
    surv.A0.list <- list()
    tlist <- seq(0, 1.68, length=200)#250)
    for (jj in 1:length(tlist)) {
        psi0.A1.jj <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                               categorical=FALSE, competing.risk=competing.risk,
                               intervention.A=1, t0=t0, tau=tlist[jj],
                               reversed.setting=reversed.setting,
                               square.effect2=square.effect2,
                               square.effect1=square.effect1,
                               interaction.Atime=interaction.Atime)
        psi0.A0.jj <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                               categorical=FALSE, competing.risk=competing.risk,
                               intervention.A=0, t0=t0, tau=tlist[jj],
                               reversed.setting=reversed.setting,
                               square.effect2=square.effect2,
                               square.effect1=square.effect1,
                               interaction.Atime=interaction.Atime)
        print(surv.list[[jj]] <- psi0.A1.jj - psi0.A0.jj)
        surv.A1.list[[jj]] <- psi0.A1.jj
        surv.A0.list[[jj]] <- psi0.A0.jj
        saveRDS(cbind(t=tlist[1:jj], surv.list),
                file=paste0("./simulation-contmle/output/",
                            "outlist-psi0-survival-function",
                            ifelse(square.effect2, "-squareL2unif", ""),
                            ifelse(square.effect1, "-squareL1unif", ""),
                            ifelse(interaction.Atime, "-interactionAtime", ""),
                            ifelse(reversed.setting, "-reversedsetting", ""),
                            ".rds"))
        saveRDS(cbind(t=tlist[1:jj], surv.A1.list),
                file=paste0("./simulation-contmle/output/",
                            "outlist-psi0-survival-A1-function",
                            ifelse(square.effect2, "-squareL2unif", ""),
                            ifelse(square.effect1, "-squareL1unif", ""),
                            ifelse(interaction.Atime, "-interactionAtime", ""),
                            ifelse(reversed.setting, "-reversedsetting", ""),
                            ".rds"))
        saveRDS(cbind(t=tlist[1:jj], surv.A0.list),
                file=paste0("./simulation-contmle/output/",
                            "outlist-psi0-survival-A0-function",
                            ifelse(square.effect2, "-squareL2unif", ""),
                            ifelse(square.effect1, "-squareL1unif", ""),
                            ifelse(interaction.Atime, "-interactionAtime", ""),
                            ifelse(reversed.setting, "-reversedsetting", ""),
                            ".rds"))
    }

}

if (FALSE) {

    source("./R/sim.data.continuous.R")
    km.est <- list()
    cox.out <- list()
    pval.out <- list()
    for (m in 1:100) {
        dt <- sim.data(n, betaA=betaA, betaL=betaL, nu=nu, eta=eta, t0=t0,
                       seed=m+100,
                       competing.risk=competing.risk,
                       categorical=FALSE, randomize.A=randomize.A,
                       censoring.informative=censoring.informative,
                       censoring=!no.censoring,
                       square.effect2=square.effect2,
                       square.effect1=square.effect1,
                       verbose=FALSE,
                       reversed.setting=reversed.setting,
                       interaction.Atime=interaction.Atime)

        hr.mod <- coxph(Surv(time, delta==1)~A, data=dt[time<=1.68])
        print(hr.mod)
        hr.pval <- summary(hr.mod)$coefficients["A", 5]
        print(hr <- exp(coef(hr.mod)["A"]))
        cox.out[[m]] <- hr
        pval.out[[m]] <- hr.pval
        

        km.est[[m]] <- contmle(dt,
                               change.point=change.point,
                               outcome.model=outcome.model,
                               cr=competing.risk, verbose=FALSE,
                               one.step=one.step, deps.size=0.001, no.small.steps=500,
                               fit.outcome=fit.outcome,
                               fit.cens=fit.cens,
                               fit.cr=fit.cr,
                               treat.effect=treat.effect,
                               tau=tau,
                               cut.time=12, V=5,
                               cut.time.A=12,
                               cut.covars=10, cut.L1.A=10, #3,
                               cut.L.interaction=3,
                               lambda.cvs=seq(0.0000001, 0.005, length=25),
                               output.km=TRUE, output.hr=TRUE, only.km=TRUE,
                               sl.models=list(mod1c=c(Surv(time, delta==1)~A+L1+L2+L3, t0=0.9),
                                              mod2c=c(Surv(time, delta==1)~A*L1.squared+L1*L2+L3, t0=NULL),
                                              #
                                              mod3=c(Surv(time, delta==1)~A*L1.squared+L1*L2+L3, t0=0.7),
                                              mod4=c(Surv(time, delta==1)~A*L1.squared, t0=0.7),
                                              mod5a=c(Surv(time, delta==1)~A+L1.squared+L2+L3, t0=0.7),
                                              mod6a=c(Surv(time, delta==1)~A+L1.squared, t0=NULL),
                                              mod7=c(Surv(time, delta==1)~A+L1+L2+L3, t0=NULL),
                                              mod8b=c(Surv(time, delta==1)~A*L1+L2+L3, t0=NULL),
                                              mod9b=c(Surv(time, delta==1)~A*L1.squared+L2+L3, t0=NULL),
                                              mod10=c(Surv(time, delta==1)~A+L1.squared, t0=0.7)
                                              ))
    }
    
    mean(unlist(km.est)-psi0)
    mean((unlist(cox.out)))
    mean(log(unlist(cox.out)))
    table(unlist(pval.out)[log(unlist(cox.out))<=0]<0.05)
    table(unlist(pval.out)[log(unlist(cox.out))>0]<0.05)
    

    
}

#-------------------------------------------------------------------------------------------#
## repeat simulations (parallelize)
#-------------------------------------------------------------------------------------------#

if (system("echo $USER",intern=TRUE)%in%c("jhl781")){ 
    no_cores <- 35
    if (fit.cens=="hal" | fit.outcome=="hal") {
        no_cores <- 5
    }
} else {
    no_cores <- detectCores() - 1
}

registerDoParallel(no_cores)

out <- foreach(m=1:M, .errorhandling="pass"
               ) %dopar% {

                   dt <- sim.data(n, betaA=betaA, betaL=betaL, nu=nu, eta=eta, t0=t0,
                                  seed=m+100,
                                  competing.risk=competing.risk,
                                  categorical=FALSE, randomize.A=randomize.A,
                                  censoring.informative=censoring.informative,
                                  censoring.high=censoring.high,
                                  censoring=!no.censoring,
                                  square.effect2=square.effect2,
                                  square.effect1=square.effect1,
                                  reversed.setting=reversed.setting,
                                  interaction.Atime=interaction.Atime)
                   
                   print(paste0("m=", m))

                   print(dt[, table(time<=tau, delta, A)])

                   return(list("est"=contmle(dt,
                                             change.point=change.point,
                                             outcome.model=outcome.model,
                                             cens.model=cens.model,
                                             cr.model=Surv(time, delta==2) ~A+L1,
                                             cr=competing.risk,
                                             one.step=one.step, #deps.size=0.0001,
                                             fit.outcome=fit.outcome,
                                             fit.cens=fit.cens,
                                             treat.effect=treat.effect,
                                             tau=tau,
                                             cut.time=10, V=5,#10,
                                             cut.time.A=10,
                                             cut.covars=8, cut.L1.A=8,#3,
                                             cut.L.interaction=3,
                                             #lambda.cvs=seq(0.0000001, 0.01, length=50),
                                             lambda.cvs=seq(0.0000001, 0.01, length=50),
                                             #lambda.cvs.cens=seq(0.001, 0.07, length=50),
                                             output.km=TRUE, output.hr=TRUE,
                                             sl.models=list(mod1d=c(Surv(time, delta==1)~A+L1+L2+L3, t0=0.9),
                                                            mod2d=c(Surv(time, delta==1)~A+L2.squared+L1*L2+L3, t0=NULL),
                                                            mod3=c(Surv(time, delta==1)~A+L1.squared+L1*L2+L3, t0=0.7),
                                                            mod4=c(Surv(time, delta==1)~A+L2.squared, t0=0.7),
                                                            mod5a=c(Surv(time, delta==1)~A+L1.squared, t0=0.3),
                                                            mod6a=c(Surv(time, delta==1)~A+L1.squared, t0=0.9),
                                                            mod5b=c(Surv(time, delta==1)~A+L1.squared+L2+L3, t0=0.7),
                                                            mod6b=c(Surv(time, delta==1)~A+L2.squared, t0=NULL),
                                                            mod7=c(Surv(time, delta==1)~A+L1+L2+L3, t0=NULL),
                                                            mod8c=c(Surv(time, delta==1)~A*L1+L2+L3, t0=NULL),
                                                            mod9c=c(Surv(time, delta==1)~A*L1.squared+L2+L3, t0=NULL),
                                                            mod10=c(Surv(time, delta==1)~A+L1.squared, t0=0.7)
                                                            ))))

               }

stopImplicitCluster()


saveRDS(out,
        file=paste0("./simulation-contmle/output/",
                    "outlist-continuous-tmle",
                    ifelse(competing.risk, "-competingrisk", ""), 
                    "-tau", tau,
                    paste0("-effect-", ifelse(treat.effect=="both", "ate", paste0("A", treat.effect))),
                    ifelse(n!=1000, paste0("-n", n), ""),
                    ifelse(one.step, "-onestep", ""),
                    ifelse(square.effect2, "-squareL2unif", ""),
                    ifelse(square.effect1, "-squareL1unif", ""),
                    ifelse(interaction.Atime, "-interactionAtime", ""),
                    ifelse(reversed.setting, "-reversedsetting", ""),
                    ifelse(randomize.A, "-randomizeA", ""),
                    ifelse(censoring.informative, "-informativeCensoring", ""),
                    ifelse(censoring.high, "-highCensoring", ""),
                    ifelse(no.censoring, "-noCensoring", ""),
                    ifelse(misspecify.outcome, "-misspecify-outcome", ""),
                    ifelse(misspecify.cens, "-misspecify-cens", ""),
                    paste0("-outcome-fit", fit.outcome),
                    paste0("-censoring-fit", fit.cens),
                    ifelse(competing.risk, paste0("-cr-fit", fit.cr), ""), 
                    "-M", M, ".rds"))


if (FALSE) {

    out <- readRDS(file=paste0("./simulation-contmle/output/",
                               "outlist-continuous-tmle",
                               ifelse(competing.risk, "-competingrisk", ""), 
                               "-tau", tau,
                               paste0("-effect-", ifelse(treat.effect=="both", "ate", paste0("A", treat.effect))),
                               ifelse(n!=1000, paste0("-n", n), ""),
                               ifelse(one.step, "-onestep", ""),
                               ifelse(square.effect2, "-squareL2unif", ""),
                               ifelse(square.effect1, "-squareL1unif", ""),
                               ifelse(interaction.Atime, "-interactionAtime", ""),
                               ifelse(reversed.setting, "-reversedsetting", ""),
                               ifelse(randomize.A, "-randomizeA", ""),
                               ifelse(censoring.informative, "-informativeCensoring", ""),
                               ifelse(censoring.high, "-highCensoring", ""),
                               ifelse(no.censoring, "-noCensoring", ""),
                               ifelse(misspecify.outcome, "-misspecify-outcome", ""),
                               ifelse(misspecify.cens, "-misspecify-cens", ""),
                               paste0("-outcome-fit", fit.outcome),
                               paste0("-censoring-fit", fit.cens),
                               ifelse(competing.risk, paste0("-cr-fit", fit.cr), ""), 
                               "-M", M, ".rds"))

    (psi0 <- readRDS(file=paste0("./simulation-contmle/output/",
                                 "save-psi0",
                                 ifelse(competing.risk, "-competingrisk", ""), 
                                 "-tau", tau, 
                                 ifelse(square.effect2, "-squareL2unif", ""),
                                 ifelse(square.effect1, "-squareL1unif", ""),
                                 ifelse(interaction.Atime, "-interactionAtime", ""),
                                 ifelse(reversed.setting, "-reversedsetting", ""),
                                 ".rds")))

    psi0.A1 <- psi0["psi0.A1"]
    psi0.A0 <- psi0["psi0.A0"]
    psi0 <- psi0["psi0"]

    if (treat.effect=="0") psi0 <- psi0.A0 else if (treat.effect=="1") psi0 <- psi0.A1
    
    init.fit <- unlist(lapply(out, function(o1) {
        o2 <- o1[["est"]]
        return(o2[[1]][1])
    }))

    mean(init.fit - psi0)
    
    tmle.fit <- unlist(lapply(out, function(o1) {
        o2 <- o1[["est"]]
        return(o2[[length(o2)]][1])
    }))

    mean(tmle.fit - psi0)

    km.fit <- unlist(lapply(out, function(o1) {
        o2 <- o1[["est"]]
        return(o2[[1]]["km.est"])
    }))

    mean(km.fit - psi0)

    mse(km.fit)
    mse(tmle.fit)

    (mse(tmle.fit) / mse(km.fit))^2

    par(mfrow=c(1,3))
    hist(init.fit); abline(v=psi0, col="red")
    hist(tmle.fit); abline(v=psi0, col="red")
    hist(km.fit); abline(v=psi0, col="red")
    par(mfrow=c(1,1))
    
    (cov <- mean(unlist(lapply(out, function(o1, true=psi0) {
        o2 <- o1[["est"]]
        if (length(o2)>1) {
            est <- try(o2[[length(o2)]][1])
            sd <- o2[[length(o2)]][2]
            return(est-1.96*sd<=true & true<=est+1.96*sd)
        } 
    }))))

    (sd.tmle <- sd(tmle.fit))

    (oracle.cov <- mean(unlist(lapply(out, function(o1, true=psi0, sd=sd.tmle) {
        o2 <- o1[["est"]]
        if (length(o2)>1) {
            est <- try(o2[[length(o2)]][1])
            return(est-1.96*sd<=true & true<=est+1.96*sd)
        }
    }))))

    (cov.km <- mean(unlist(lapply(out, function(o1, true=psi0) {
        o2 <- o1[["est"]]
        if (length(o2)>1) {
            est <- try(o2[[1]]["km.est"])
            sd <- o2[[1]]["km.se"]
            return(est-1.96*sd<=true & true<=est+1.96*sd)
        } 
    }))))

    sd.km <- sd(km.fit)
    (oracle.cov.km <- mean(unlist(lapply(out, function(o1, true=psi0, sd=sd.km) {
        o2 <- o1[["est"]]
        if (length(o2)>1) {
            est <- try(o2[[1]]["km.est"])
            #sd <- o2[[1]]["km.se"]
            return(est-1.96*sd<=true & true<=est+1.96*sd)
        } 
    }))))




}




if (FALSE) {

    out1 <- readRDS(file=paste0("./simulation-contmle/output/",
                               "outlist-continuous-tmle",
                               ifelse(competing.risk, "-competingrisk", ""), 
                               "-tau", tau,
                               paste0("-effect-", paste0("A", "1")),
                               ifelse(n!=1000, paste0("-n", n), ""),
                               ifelse(one.step, "-onestep", ""),
                               ifelse(square.effect2, "-squareL2unif", ""),
                               ifelse(square.effect1, "-squareL1unif", ""),
                               ifelse(interaction.Atime, "-interactionAtime", ""),
                               ifelse(reversed.setting, "-reversedsetting", ""),
                               ifelse(randomize.A, "-randomizeA", ""),
                               ifelse(censoring.informative, "-informativeCensoring", ""),
                               ifelse(censoring.high, "-highCensoring", ""),
                               ifelse(no.censoring, "-noCensoring", ""),
                               ifelse(misspecify.outcome, "-misspecify-outcome", ""),
                               ifelse(misspecify.cens, "-misspecify-cens", ""),
                               paste0("-outcome-fit", fit.outcome),
                               paste0("-censoring-fit", fit.cens),
                               ifelse(competing.risk, paste0("-cr-fit", fit.cr), ""), 
                               "-M", M, ".rds"))

    out0 <- readRDS(file=paste0("./simulation-contmle/output/",
                                "outlist-continuous-tmle",
                                ifelse(competing.risk, "-competingrisk", ""), 
                                "-tau", tau,
                                paste0("-effect-", paste0("A", "0")),
                                ifelse(n!=1000, paste0("-n", n), ""),
                                ifelse(one.step, "-onestep", ""),
                                ifelse(square.effect2, "-squareL2unif", ""),
                                ifelse(square.effect1, "-squareL1unif", ""),
                                ifelse(interaction.Atime, "-interactionAtime", ""),
                                ifelse(reversed.setting, "-reversedsetting", ""),
                                ifelse(randomize.A, "-randomizeA", ""),
                                ifelse(censoring.informative, "-informativeCensoring", ""),
                                ifelse(censoring.high, "-highCensoring", ""),
                                ifelse(no.censoring, "-noCensoring", ""),
                                ifelse(misspecify.outcome, "-misspecify-outcome", ""),
                                ifelse(misspecify.cens, "-misspecify-cens", ""),
                                paste0("-outcome-fit", fit.outcome),
                                paste0("-censoring-fit", fit.cens),
                                ifelse(competing.risk, paste0("-cr-fit", fit.cr), ""), 
                                "-M", M, ".rds"))

    (psi0 <- readRDS(file=paste0("./simulation-contmle/output/",
                                 "save-psi0",
                                 ifelse(competing.risk, "-competingrisk", ""), 
                                 "-tau", tau, 
                                 ifelse(square.effect2, "-squareL2unif", ""),
                                 ifelse(square.effect1, "-squareL1unif", ""),
                                 ifelse(interaction.Atime, "-interactionAtime", ""),
                                 ifelse(reversed.setting, "-reversedsetting", ""),
                                 ".rds")))

    psi0.A1 <- psi0["psi0.A1"]
    psi0.A0 <- psi0["psi0.A0"]
    psi0 <- psi0["psi0"]

    tmle.fit1 <- unlist(lapply(out1, function(o1) {
        o2 <- o1[["est"]]
        return(o2[[length(o2)]][1])
    }))

    mean(tmle.fit1 - psi0.A1)

    km.fit1 <- unlist(lapply(out1, function(o1) {
        o2 <- o1[["est"]]
        return(o2[[1]]["km.est"])
    }))

    mean(km.fit1 - psi0.A1)

    sd.tmle.fit1 <- unlist(lapply(out1, function(o1) {
        o2 <- o1[["est"]]
        return(o2[[length(o2)]]["sd.eic"])
    }))
    
    sd.km.fit1 <- unlist(lapply(out1, function(o1) {
        o2 <- o1[["est"]]
        return(o2[[1]]["km.se"])
    }))
    
    tmle.fit0 <- unlist(lapply(out0, function(o1) {
        o2 <- o1[["est"]]
        return(o2[[length(o2)]][1])
    }))

    mean(tmle.fit0 - psi0.A0)

    km.fit0 <- unlist(lapply(out0, function(o1) {
        o2 <- o1[["est"]]
        return(o2[[1]]["km.est"])
    }))

    mean(km.fit0 - psi0.A0)

    sd.tmle.fit0 <- unlist(lapply(out0, function(o1) {
        o2 <- o1[["est"]]
        return(o2[[length(o2)]]["sd.eic"])
    }))
    
    sd.km.fit0 <- unlist(lapply(out0, function(o1) {
        o2 <- o1[["est"]]
        return(o2[[1]]["km.se"])
    }))

    mse(tmle.fit1) / mse(km.fit1)
    mse(tmle.fit0) / mse(km.fit0)

    cov.fun <- function(fit=tmle.fit1, sd=sd.tmle.fit1, true=psi0.A1) {
        return(fit-1.96*sd<=true & true<=fit+1.96*sd)
    }

    (cov.tmle1 <- mean(cov.fun(fit=tmle.fit1, sd=sd.tmle.fit1, true=psi0.A1)))
    (cov.tmle0 <- mean(cov.fun(fit=tmle.fit0, sd=sd.tmle.fit0, true=psi0.A0)))
    (cov.km1 <- mean(cov.fun(fit=km.fit1, sd=sd.km.fit1, true=psi0.A1)))
    (cov.km0 <- mean(cov.fun(fit=km.fit0, sd=sd.km.fit0, true=psi0.A0)))

    
    (cov.tmle <- mean(cov.fun(fit=tmle.fit1-tmle.fit0, sd=sqrt(sd.tmle.fit1^2+sd.tmle.fit0^2),
                              true=psi0)))
    (cov.km <- mean(cov.fun(fit=km.fit1-km.fit0, sd=sqrt(sd.km.fit1^2+sd.km.fit0^2),
                            true=psi0)))

}
