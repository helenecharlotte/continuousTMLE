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

#source("./oct-7/R/sim.data.continuous.R")
source("./R/sim.data.continuous.R")
#source("./sep-24-1/R/cox.tmle.R")
source("./R/cox.tmle.R")
source("./R/cox.super.learner.R")
source("./R/poisson.hal.R")
source("./R/poisson.hal.sl.R")

#-------------------------------------------------------------------------------------------#
## set parameters for simulations
#-------------------------------------------------------------------------------------------#

M <- 291#290#500#250#350#310#10#500#300#250#500#300#500#100#501#492#530#520#500#500#500#1#500#100#510#220#1#220#250#100#200#101#100#300#100#500#100#500#500#500#502#500#120#502#150#100#501#500#501#500
n <- 1000
get.truth <- FALSE#FALSE#TRUE

betaA <- -0.15
betaL <- 1.1
nu <- 1.7
nu1 <- FALSE
eta <- 0.7
t0 <- 0.9#0.4
tau <- 1.2

interaction.AL <- FALSE#TRUE
square.effect <- FALSE
square.effect2 <- TRUE#FALSE#TRUE#TRUE#FALSE#TRUE#FALSE#TRUE#FALSE#TRUE#FALSE#TRUE#FALSE#TRUE#FALSE#TRUE#FALSE#TRUE#TRUE#TRUE#TRUE#TRUE#FALSE#TRUE
interaction.Atime <- TRUE#TRUE#TRUE#FALSE#TRUE#FALSE#TRUE#FALSE#TRUE
new.setting <- TRUE#FALSE#TRUE#TRUE#TRUE

competing.risk <- FALSE

randomize.A <- TRUE
censoring.informative <- TRUE#TRUE#FALSE#TRUE#FALSE

only.cox.sl <- TRUE#FALSE#TRUE#FALSE#TRUE#FALSE#TRUE
only.stochastic <- FALSE#TRUE#FALSE#FALSE#TRUE
one.step <- FALSE#TRUE#FALSE#TRUE#FALSE
save.hals.for.plotting <- TRUE#FALSE#FALSE#TRUE#FALSE
only.diff <- TRUE#FALSE#TRUE

misspecify.Y <- TRUE#FALSE#FALSE#FALSE#FALSE#TRUE#FALSE#TRUE#TRUE#TRUE#TRUE#FALSE#TRUE
fit.km <- TRUE
output.hr <- TRUE
sl <- TRUE#TRUE#FALSE#TRUE#FALSE#TRUE
sl.method <- 2#2#1
poisson.initial <- TRUE#FALSE#TRUE#FALSE#TRUE
lambda.cv <- NULL
penalize.time <- FALSE
adjust.penalization <- TRUE
sl.poisson <- TRUE#TRUE#FALSE#FALSE#TRUE#FALSE

cens.misspecified <- FALSE
km.cens <- FALSE
sl.cens <- FALSE#TRUE
poisson.cens <- FALSE

#-------------------------------------------------------------------------------------------#
## a bit of reprocessing
#-------------------------------------------------------------------------------------------#

if (interaction.Atime) betaA <- -0.7
if (nu1) nu <- 1

change.point <- NULL

if (new.setting) {
    betaA <- 0.5#0.3
    t0 <- 0.7
    #eta <- 0.5
    cvot.setting <- TRUE
    #; t0 <- 0.8; tau <- 1.2; nu <- 1.4; eta <- 0.4;
    #if (interaction.Atime & square.effect2) {
    #    nu <- 1.7
    #    eta <- 0.3
    #}
} else {
    cvot.setting <- FALSE
}

if (cens.misspecified) {
    cens.model <- Surv(time, delta==0)~1#L1
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
        if (square.effect2) {
            outcome.model <- Surv(time, delta==1)~A+L1+L2.squared+L2+L3
            outcome.model <- Surv(time, delta==1)~A+L1+L2.squared+L3
        } else {
            outcome.model <- Surv(time, delta==1)~A+L1+L2+L3
        }
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

    #print(psi0.stochastic <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
    #                                  categorical=FALSE, competing.risk=competing.risk,
    #                                  intervention.A=function(L) 0.5*(L[,1]<0.5),#0.3*(L[,2]>0.3),
    #                                  t0=t0, tau=tau,
    #                                  square.effect=square.effect,
    #                                  interaction.Atime=interaction.Atime, verbose=TRUE,
    #                                  square.effect2=square.effect2,
    #                                  interaction.AL=interaction.AL))
    
    print(psi0.A1 <- sim.data(1e7, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                              categorical=FALSE, #competing.risk=competing.risk,
                              intervention.A=1, t0=t0, tau=tau,
                              square.effect=square.effect,
                              interaction.Atime=interaction.Atime, verbose=TRUE,
                              square.effect2=square.effect2,
                              cvot.setting=cvot.setting,
                              interaction.AL=interaction.AL))
    print(psi0.A0 <- sim.data(1e7, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                              categorical=FALSE, #competing.risk=competing.risk,
                              intervention.A=0, t0=t0, tau=tau,
                              square.effect=square.effect,
                              square.effect2=square.effect2,
                              cvot.setting=cvot.setting,
                              interaction.Atime=interaction.Atime, verbose=TRUE,
                              interaction.AL=interaction.AL))
    print(psi0 <- psi0.A1 - psi0.A0) 

    saveRDS(c(psi0=psi0, psi0.A1=psi0.A1, psi0.A0=psi0.A0),
            file=paste0("./simulation-cox-tmle/output/",
                        "outlist-psi0",
                        ifelse(competing.risk, "-competingrisk", ""), 
                        "-tau", tau, 
                        ifelse(interaction.AL, "-interactionAL", ""),
                        ifelse(interaction.Atime, "-interactionAtime", ""),
                        ifelse(square.effect, "-squareL2", ""),
                        ifelse(square.effect2, "-squareL2unif", ""),
                        ifelse(new.setting, "-newsetting", ""),
                        ifelse(nu1, "-nu1", ""),
                        ".rds"))

    surv.list <- list()
    surv.A1.list <- list()
    surv.A0.list <- list()
    tlist <- seq(0, 1.68, length=100)#250)
    for (jj in 1:length(tlist)) {
        psi0.A1.jj <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                               categorical=FALSE, competing.risk=competing.risk,
                               intervention.A=1, t0=t0, tau=tlist[jj],
                               square.effect=square.effect,
                               interaction.Atime=interaction.Atime, verbose=TRUE,
                               cvot.setting=cvot.setting,
                               square.effect2=square.effect2,
                               interaction.AL=interaction.AL)
        psi0.A0.jj <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                               categorical=FALSE, competing.risk=competing.risk,
                               intervention.A=0, t0=t0, tau=tlist[jj],
                               square.effect=square.effect,
                               interaction.Atime=interaction.Atime, verbose=TRUE,
                               cvot.setting=cvot.setting,
                               square.effect2=square.effect2,
                               interaction.AL=interaction.AL)
        print(surv.list[[jj]] <- psi0.A1.jj - psi0.A0.jj)
        surv.A1.list[[jj]] <- psi0.A1.jj
        surv.A0.list[[jj]] <- psi0.A0.jj
        saveRDS(cbind(t=tlist[1:jj], surv.list),
                file=paste0("./simulation-cox-targeting/output/",
                            "outlist-psi0-survival-function",
                            ifelse(interaction.AL, "-interactionAL", ""),
                            ifelse(interaction.Atime, "-interactionAtime", ""),
                            ifelse(square.effect, "-squareL2", ""),
                            ifelse(square.effect2, "-squareL2unif", ""),
                            ifelse(new.setting, "-newsetting", ""),
                            "oct6.rds"))
        saveRDS(cbind(t=tlist[1:jj], surv.A1.list),
                file=paste0("./simulation-cox-targeting/output/",
                            "outlist-psi0-survival-A1-function",
                            ifelse(interaction.AL, "-interactionAL", ""),
                            ifelse(interaction.Atime, "-interactionAtime", ""),
                            ifelse(square.effect2, "-squareL2unif", ""),
                            ifelse(new.setting, "-newsetting", ""),
                            "oct6.rds"))
        saveRDS(cbind(t=tlist[1:jj], surv.A0.list),
                file=paste0("./simulation-cox-targeting/output/",
                            "outlist-psi0-survival-A0-function",
                            ifelse(interaction.AL, "-interactionAL", ""),
                            ifelse(interaction.Atime, "-interactionAtime", ""),
                            ifelse(square.effect, "-squareL2", ""),
                            ifelse(square.effect2, "-squareL2unif", ""),
                            ifelse(new.setting, "-newsetting", ""),
                            "oct6.rds"))
    }


}

if (FALSE) {

    cox.out <- list()
    pval.out <- list()
    #nu <- 1.7; eta <- 0.3; censoring.informative <- TRUE
    for (m in 1:30) {
        source("./R/sim.data.continuous.R")
        dt <- sim.data(n, betaA=betaA, betaL=betaL, nu=nu, eta=eta, t0=t0,
                       seed=m+100,
                       competing.risk=competing.risk,
                       categorical=FALSE, randomize.A=randomize.A,
                       censoring.informative=censoring.informative,
                       interaction.Atime=interaction.Atime,
                       square.effect=square.effect,
                       cvot.setting=cvot.setting,
                       square.effect2=square.effect2,
                       interaction.AL=interaction.AL)

        hr.mod <- coxph(Surv(time, delta==1)~A, data=dt)
        print(hr.mod)
        hr.pval <- summary(hr.mod)$coefficients["A", 5]
        print(hr <- exp(coef(hr.mod)["A"]))
        cox.out[[m]] <- hr
        pval.out[[m]] <- hr.pval
    }

    mean(unlist(cox.out))
    mean(log(unlist(cox.out)))
    table(unlist(pval.out)<0.05)

}

#-------------------------------------------------------------------------------------------#
## repeat simulations (parallelize)
#-------------------------------------------------------------------------------------------#

if (system("echo $USER",intern=TRUE)%in%c("jhl781")){ 
    no_cores <- 5#35#30#35#30#25#30#25#30#20#5#23#12
    if (poisson.initial) {
        no_cores <- 25
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
                                  interaction.Atime=interaction.Atime,
                                  square.effect=square.effect,
                                  cvot.setting=cvot.setting,
                                  square.effect2=square.effect2,
                                  interaction.AL=interaction.AL)
                   
                   print(paste0("m=", m))

                   if (only.cox.sl) {
                       return(cox.tmle(dt, change.point=change.point, browse=FALSE, only.cox.sl=TRUE,
                                       poisson.initial=poisson.initial,
                                       adjust.penalization=adjust.penalization,
                                       cut.time=12,
                                       cut.time.A=8,
                                       sl.method=sl.method,
                                       poisson.cens=poisson.cens,
                                       lambda.cv=lambda.cv,
                                       cut.covars=8, cut.L1.A=8,
                                       cut.L.interaction=6,
                                       tau=tau, km.cens=km.cens,
                                       cens.model=cens.model,
                                       sl.poisson=sl.poisson,
                                       sl.cens=sl.cens,
                                       output.km=fit.km, output.hr=output.hr,
                                       penalize.time=penalize.time, 
                                       outcome.model=outcome.model, mod.period2=mod.period2,
                                       sl=sl))
                   }

                   if (only.stochastic) {
                       return(list("stochastic"=cox.tmle(dt, change.point=change.point, browse=FALSE,
                                                         one.step=one.step, #deps.size=0.0001,
                                                         poisson.initial=poisson.initial,
                                                         adjust.penalization=adjust.penalization,
                                                         cut.time=12, 
                                                         cut.time.A=8, 
                                                         poisson.cens=poisson.cens,
                                                         lambda.cv=lambda.cv,
                                                         treat.effect="stochastic",
                                                         pi.star.fun=function(A, L) 0.5*(L[,1]<0.5),#0.3*(L[,2]>0.3),#0.7*(L[,2]>0),
                                                         cut.covars=8, cut.L1.A=8,
                                                         cut.L.interaction=6,
                                                         tau=tau, km.cens=km.cens,
                                                         cens.model=cens.model,
                                                         sl.poisson=sl.poisson, sl=sl, sl.method=sl.method,
                                                         sl.cens=sl.cens,
                                                         output.km=fit.km, output.hr=output.hr,
                                                         penalize.time=penalize.time, 
                                                         outcome.model=outcome.model, mod.period2=mod.period2)))
                   }

                   if (save.hals.for.plotting) {
                       cox.tmle(dt, change.point=change.point, browse=FALSE,
                                one.step=one.step, #deps.size=0.0001,
                                poisson.initial=TRUE,
                                adjust.penalization=adjust.penalization,
                                cut.time=30, 
                                cut.time.A=15,
                                poisson.cens=poisson.cens,
                                lambda.cv=lambda.cv,
                                treat.effect="both",
                                cut.covars=10, cut.L1.A=8,
                                cut.L.interaction=6,
                                tau=tau, km.cens=km.cens,
                                cens.model=cens.model,
                                #lambda.cvs=seq(0.0001, 0.01, length=25),
                                #lambda.cvs=seq(0.00001, 0.001, length=25),#15),#10),
                                #lambda.cvs=seq(0.00001, 0.01, length=5),#15),#10),
                                lambda.cvs=seq(0.00001, 0.01, length=50)[1:13],#15),#10),
                                #lambda.cvs=c(0.05, 0.01, 0.001, 0.004, 0.003, 0.0001, 0.00001),#15),#10),
                                #lambda.cvs=rev(c(0.01, 0.0075, 0.005, 0.0025, 0.001, 0.00075,
                                #                 0.0005, 0.00025, 0.00075, 0.0005, 0.00025,
                                #                 0.0001, 0.000075, 0.00005, 0.000025, 0.00001)),#15),#10),
                                #lambda.cvs=rev(c(0.008, 0.005, 0.003, 0.001, 0.0008, 0.0006, 0.0004, 0.0002, 0.0001, 0.00005, 0.00001)),#15),#10),
                                sl.poisson=sl.poisson, sl=sl, sl.method=sl.method,
                                sl.cens=sl.cens,
                                save.hals.for.plotting=TRUE,
                                output.km=fit.km, output.hr=output.hr,
                                penalize.time=FALSE,#TRUE, 
                                outcome.model=outcome.model, mod.period2=mod.period2)
                   }

                   (list("diff"=cox.tmle(dt, change.point=change.point, browse=FALSE,
                                                   one.step=one.step, #deps.size=0.0001,
                                                   poisson.initial=poisson.initial,
                                                   adjust.penalization=adjust.penalization,
                                                   cut.time=3, V=3,
                                                   cut.time.A=3,
                                                   poisson.cens=poisson.cens,
                                                   lambda.cv=lambda.cv, #0.0000001000
                                                   treat.effect="both",
                                                   cut.covars=3, cut.L1.A=3,
                                                   cut.L.interaction=6,
                                                   tau=tau, km.cens=km.cens,
                                                   cens.model=cens.model,
                                                   #lambda.cvs=seq(0.0001, 0.008, length=8),
                                                   #lambda.cvs=seq(0.00005, 0.005, length=10),
                                                   #lambda.cvs=seq(0.000001, 0.001, length=10),
                                                   #lambda.cvs=seq(0.0000001, 0.01, length=50),
                                                   #lambda.cvs=seq(0.001, 0.1, length=10),
                                                   #lambda.cvs=seq(0.0000001, 0.01, length=30),
                                                   lambda.cvs=seq(0.01, 0.002, length=5),
                                                   sl.poisson=sl.poisson, sl=sl, sl.method=sl.method,
                                                   sl.cens=sl.cens,
                                                   output.km=fit.km, output.hr=output.hr,
                                                   penalize.time=penalize.time, 
                                                   outcome.model=outcome.model, mod.period2=mod.period2)))

                   if (only.diff) {
                       return(list("diff"=cox.tmle(dt, change.point=change.point, browse=FALSE,
                                                   one.step=one.step, #deps.size=0.0001,
                                                   poisson.initial=poisson.initial,
                                                   adjust.penalization=adjust.penalization,
                                                   cut.time=15, V=5,
                                                   cut.time.A=12,
                                                   poisson.cens=poisson.cens,
                                                   lambda.cv=lambda.cv, #0.0000001000
                                                   treat.effect="both",
                                                   cut.covars=10, cut.L1.A=8,
                                                   cut.L.interaction=6,
                                                   tau=tau, km.cens=km.cens,
                                                   cens.model=cens.model,
                                                   #lambda.cvs=seq(0.0001, 0.008, length=8),
                                                   #lambda.cvs=seq(0.00005, 0.005, length=10),
                                                   #lambda.cvs=seq(0.000001, 0.001, length=10),
                                                   #lambda.cvs=seq(0.0000001, 0.01, length=50),
                                                   #lambda.cvs=seq(0.001, 0.1, length=10),
                                                   #lambda.cvs=seq(0.0000001, 0.01, length=30),
                                                   lambda.cvs=seq(0.0000001, 0.002, length=30),
                                                   sl.poisson=sl.poisson, sl=sl, sl.method=sl.method,
                                                   sl.cens=sl.cens,
                                                   output.km=fit.km, output.hr=output.hr,
                                                   penalize.time=penalize.time, 
                                                   outcome.model=outcome.model, mod.period2=mod.period2)))
                   }
                   
                   return(
                       list(#-- estimate difference:
                           "diff"=cox.tmle(dt, change.point=change.point, browse=FALSE,
                                           one.step=one.step, #deps.size=0.0001,
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
                                           #lambda.cvs=seq(0.0015, 0.008, length=51),
                                           sl.poisson=sl.poisson, sl=sl, sl.method=sl.method,
                                           sl.cens=sl.cens,
                                           output.km=fit.km, output.hr=output.hr,
                                           penalize.time=penalize.time, 
                                           outcome.model=outcome.model, mod.period2=mod.period2),
                           #-- estimate effect A=1: 
                           "A=1"=cox.tmle(dt, change.point=change.point, browse=FALSE,
                                          one.step=one.step, #deps.size=0.0001,
                                          poisson.initial=poisson.initial,
                                          adjust.penalization=adjust.penalization,
                                          cut.time=12,
                                          cut.time.A=8,
                                          #lambda.cvs=seq(0.0015, 0.008, length=51),
                                          poisson.cens=poisson.cens,
                                          lambda.cv=lambda.cv,
                                          cut.covars=8, cut.L1.A=8,
                                          cut.L.interaction=6,
                                          tau=tau, km.cens=km.cens,
                                          cens.model=cens.model,
                                          sl.poisson=sl.poisson,
                                          sl.cens=sl.cens,
                                          output.km=fit.km, output.hr=output.hr,
                                          penalize.time=penalize.time, 
                                          outcome.model=outcome.model, mod.period2=mod.period2,
                                          sl=sl, sl.method=sl.method),
                           #-- estimate effect A=0: 
                           "A=0"=cox.tmle(dt, change.point=change.point, browse2=FALSE,
                                          one.step=one.step, #deps.size=0.0001,
                                          poisson.initial=poisson.initial,
                                          adjust.penalization=adjust.penalization,
                                          cut.time=8, cut.time.A=8,
                                          #lambda.cvs=seq(0, 0.008, length=6)[-1],
                                          poisson.cens=poisson.cens, 
                                          lambda.cv=lambda.cv, 
                                          cut.covars=6, cut.L1.A=8,
                                          #lambda.cvs=seq(0.0055, 0.02, length=51)[-(1:5)],
                                          cut.L.interaction=6,
                                          treat.effect="0",
                                          tau=tau, km.cens=km.cens,
                                          cens.model=cens.model,
                                          sl.poisson=sl.poisson,
                                          sl.cens=sl.cens,
                                          output.km=fit.km, output.hr=output.hr,
                                          penalize.time=penalize.time,
                                          outcome.model=outcome.model, mod.period2=mod.period2,
                                          sl=sl, sl.method=sl.method))
                   )

               }

stopImplicitCluster()


saveRDS(out,
        file=paste0("./simulation-cox-tmle/output/",
                    "outlist-run-cox-tmle",
                    ifelse(competing.risk, "-competingrisk", ""), 
                    "-tau", tau, 
                    ifelse(n!=1000, paste0("-n", n), ""),
                    ifelse(only.stochastic, "-stochasticeffect", ""),
                    ifelse(one.step, "-onestep", ""),
                    ifelse(interaction.AL, "-interactionAL", ""),
                    ifelse(interaction.Atime, "-interactionAtime", ""),
                    ifelse(square.effect, "-squareL2", ""),
                    ifelse(square.effect2, "-squareL2unif", ""),
                    ifelse(new.setting, "-newsetting", ""),
                    ifelse(nu1, "-nu1", ""),
                    ifelse(randomize.A, "-randomizeA", ""),
                    ifelse(censoring.informative, "-informativeCensoring", ""),
                    ifelse(misspecify.Y, "-misspecifyY", ""),
                    ifelse(length(lambda.cv)>0, "-no-penalization", ""),
                    ifelse(km.cens, "-KMcensoring", ""),
                    ifelse(cens.misspecified, "-censoringMisspecified", ""),
                    ifelse(sl, "-useSL", ""),
                    ifelse(sl & only.cox.sl, "-onlySL", ""),
                    ifelse(sl & !only.cox.sl & sl.method==1, "-devianceResiduals", ""), 
                    ifelse(sl.cens, "-useSLforCensoring", ""),
                    ifelse(poisson.initial, "-testPoisson", ""),
                    ifelse(sl.poisson, "-useSLwithPoisson", ""),
                    ifelse(poisson.initial & !adjust.penalization, "-noAdjustPenalization", ""),
                    ifelse(poisson.initial & !penalize.time, "-noPenalizationTime", ""),
                    ifelse(poisson.cens, "-testPoissonCensoring", ""),
                    "-M", M, ".rds"))


if (FALSE) {

    out <- readRDS(file=paste0("./simulation-cox-tmle/output/",
                               "outlist-run-cox-tmle",
                               ifelse(competing.risk, "-competingrisk", ""), 
                               "-tau", tau, 
                               ifelse(n!=1000, paste0("-n", n), ""),
                               ifelse(only.stochastic, "-stochasticeffect", ""),
                               ifelse(one.step, "-onestep", ""),
                               ifelse(interaction.AL, "-interactionAL", ""),
                               ifelse(interaction.Atime, "-interactionAtime", ""),
                               ifelse(square.effect, "-squareL2", ""),
                               ifelse(square.effect2, "-squareL2unif", ""),
                               ifelse(new.setting, "-newsetting", ""),
                               ifelse(nu1, "-nu1", ""),
                               ifelse(randomize.A, "-randomizeA", ""),
                               ifelse(censoring.informative, "-informativeCensoring", ""),
                               ifelse(misspecify.Y, "-misspecifyY", ""),
                               ifelse(length(lambda.cv)>0, "-no-penalization", ""),
                               ifelse(km.cens, "-KMcensoring", ""),
                               ifelse(cens.misspecified, "-censoringMisspecified", ""),
                               ifelse(sl, "-useSL", ""),
                               ifelse(sl & only.cox.sl, "-onlySL", ""),
                               ifelse(sl & !only.cox.sl & sl.method==1, "-devianceResiduals", ""), 
                               ifelse(sl.cens, "-useSLforCensoring", ""),
                               ifelse(poisson.initial, "-testPoisson", ""),
                               ifelse(sl.poisson, "-useSLwithPoisson", ""),
                               ifelse(poisson.initial & !adjust.penalization, "-noAdjustPenalization", ""),
                               ifelse(poisson.initial & !penalize.time, "-noPenalizationTime", ""),
                               ifelse(poisson.cens, "-testPoissonCensoring", ""),
                               "-M", M, ".rds"))

    (psi0 <- readRDS(file=paste0("./simulation-cox-tmle/output/",
                                 "outlist-psi0",
                                 ifelse(competing.risk, "-competingrisk", ""), 
                                 "-tau", tau, 
                                 ifelse(interaction.AL, "-interactionAL", ""),
                                 ifelse(interaction.Atime, "-interactionAtime", ""),
                                 ifelse(square.effect, "-squareL2", ""),
                                 ifelse(square.effect2, "-squareL2unif", ""),
                                 ifelse(new.setting, "-newsetting", ""),
                                 ifelse(nu1, "-nu1", ""),
                                 ".rds")))

    psi0.A1 <- psi0["psi0.A1"]
    psi0.A0 <- psi0["psi0.A0"]
    psi0 <- psi0["psi0"]

    if (FALSE) { 
        par(mfrow=c(1,2))
        hist(psiA1.est <- unlist(lapply(out, function(out1) {
            out2 <- out1[["stochastic"]]
            if (length(out2)>1) return(out2[[1]][1])
        })), xlab="", ylab="", main="initial: stochastic")
        abline(v=psi0.stochastic, col="red")
        hist(psiA1.est <- unlist(lapply(out, function(out1) {
            out2 <- out1[["stochastic"]]
            if (length(out2)>1) return(out2[[length(out2)]][1])
        })), xlab="", ylab="", main="tmle: stochastic")
        abline(v=psi0.stochastic, col="red")
        (cov.stochastic <- mean(unlist(lapply(out, function(out1, true=psi0.stochastic) {
            out2 <- out1[["stochastic"]]
            if (length(out2)>1) {
                est <- out2[[length(out2)]][1]
                sd <- out2[[length(out2)]][2]
                #sd <- sd(psiA1.est)
                return(est-1.96*sd<=true & true<=est+1.96*sd)
            }
        }))))
    }

    if (FALSE) {
        par(mfrow=c(1,2))
        hist(psiA1.one <- unlist(lapply(out.one, function(out1) {
            out2 <- out1[["A=1"]]
            if (length(out2)>1) return(out2[[length(out2)]][1])
        })), xlab="", ylab="", main="one step: A=1")
        abline(v=psi0.A1, col="red")
        hist(psiA1.iter <- unlist(lapply(out.iter, function(out1) {
            out2 <- out1[["A=1"]]
            if (length(out2)>1) return(out2[[length(out2)]][1])
        })), xlab="", ylab="", main="iterative tmle: A=1")
        abline(v=psi0.A1, col="red")
    }
    
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





