
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
library(prodlim)
library(survival)
library(riskRegression)
library(Matrix)
library(coefplot)
library(hdnom)


if (system("echo $USER",intern=TRUE)%in%c("jhl781")){ 
    no_cores <- 30
} else {
    no_cores <- detectCores() - 1
}

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

source("./simulation-contmle/run.fun.R")
source("./R/sim.data.continuous.R")
source("./R/contmle.R")
source("./R/cox.super.learner.R")
source("./R/poisson.hal.R")
source("./R/poisson.hal.sl.R")

#-------------------------------------------------------------------------------------------#
## run: competing risks
#-------------------------------------------------------------------------------------------#

M <- 500#300#563#562#500#60#100#600#550#250
n <- 1000

setting <- 2
censoring.informative <- FALSE
tau <- #c(0.2,
    c(0.35, 0.5) 
iterative <- TRUE
target <- 1:2
competing.risk <- TRUE

if (FALSE) {
    out <- run.fun(M=1, n=n, competing.risk=competing.risk, target=target, tau=tau,
                   setting=setting,
                   censoring.informative=!censoring.informative,
                   verbose=TRUE,
                   weighted.norm="sigma",
                   iterative=iterative, no_cores=no_cores, save.output=FALSE)
    out <- run.fun(M=1, n=n, competing.risk=competing.risk, target=target[1], tau=tau,
                   setting=setting,
                   censoring.informative=censoring.informative,
                   verbose=TRUE,
                   iterative=iterative, no_cores=no_cores, save.output=FALSE)
    out <- run.fun(M=1, n=n, competing.risk=competing.risk, target=target[1], tau=tau[1],
                   setting=setting,
                   censoring.informative=censoring.informative,
                   verbose=TRUE,
                   iterative=iterative, no_cores=no_cores, save.output=FALSE)
    out <- run.fun(M=1, n=n, competing.risk=competing.risk, target=target, tau=tau,
                   setting=setting,
                   censoring.informative=censoring.informative,
                   verbose=TRUE,
                   iterative=iterative, no_cores=no_cores, save.output=FALSE)
}
if (FALSE) { #compute true!
    run.fun(M=1, n=n, competing.risk=competing.risk, target=target, tau=tau, setting=setting,
            get.truth=TRUE, no_cores=no_cores)
}

if (TRUE) {
    run.fun(M=M, n=n, competing.risk=competing.risk, target=target, tau=tau, setting=setting,
            censoring.informative=censoring.informative,
            iterative=iterative, no_cores=no_cores)

    run.fun(M=M, n=n, competing.risk=competing.risk, target=target, tau=tau, setting=setting,
            censoring.informative=censoring.informative, misspecify.outcome=TRUE,
            iterative=iterative, no_cores=no_cores)

    run.fun(M=M, n=n, competing.risk=competing.risk, target=target, tau=tau, setting=setting,
            censoring.informative=censoring.informative,
            fit.outcome="sl", fit.cens="sl", fit.cr="sl",
            iterative=iterative, no_cores=no_cores)
}

setting <- 2
censoring.informative <- TRUE
tau <- 0.5#c(0.2, 0.35, 0.5) 
iterative <- TRUE
target <- 1:2
competing.risk <- TRUE

if (TRUE) {
    run.fun(M=M, n=n, competing.risk=competing.risk, target=target, tau=tau, setting=setting,
            censoring.informative=censoring.informative,
            iterative=iterative, no_cores=no_cores)

    run.fun(M=M, n=n, competing.risk=competing.risk, target=target, tau=tau, setting=setting,
            censoring.informative=censoring.informative, misspecify.outcome=TRUE,
            iterative=iterative, no_cores=no_cores)

    run.fun(M=M, n=n, competing.risk=competing.risk, target=target, tau=tau, setting=setting,
            censoring.informative=censoring.informative,
            fit.outcome="sl", fit.cens="sl", fit.cr="sl",
            iterative=iterative, no_cores=no_cores)
}



if (FALSE) {
    out <- run.fun(M=1, n=n, competing.risk=TRUE, target=target, tau=tau, setting=setting,
                   iterative=iterative, no_cores=no_cores,
                   save.output=FALSE)
}


if (FALSE) {


    if (FALSE) {

        run.fun(M=1, n=n, competing.risk=TRUE, tau=tau, setting=setting,
                get.truth=TRUE, no_cores=no_cores)

    }

    if (FALSE) {
        run1 <- run.fun(M=1, n=n, competing.risk=TRUE,
                        tau=c(0.2,0.3), setting=setting, censoring.informative=TRUE,
                        save.output=FALSE, verbose=TRUE,
                        one.step=one.step, no_cores=no_cores)

    }

    #----- informative censoring ------#

    #-- cox models for all parts
    if (FALSE) {
        run.fun(M=M, n=n, competing.risk=TRUE, tau=tau, setting=setting, censoring.informative=TRUE,
                one.step=one.step, no_cores=no_cores)
    }


    if (FALSE) {
        run.fun(M=M, n=n, competing.risk=TRUE, tau=tau, setting=setting, censoring.informative=TRUE,
                one.step=one.step, no_cores=no_cores,
                misspecify.outcome=TRUE)
    }

    if (TRUE) {
        #-- sl for all parts
        run.fun(M=M, n=n, competing.risk=TRUE, tau=tau, setting=setting, censoring.informative=TRUE,
                one.step=one.step, no_cores=no_cores,
                fit.outcome="sl",
                fit.cens="sl",
                fit.cr="sl")
    }

    #----- uninformative censoring ------#

    if (FALSE) {

        #-- cox models for all parts
        run.fun(M=M, n=n, competing.risk=TRUE, tau=tau, setting=setting, censoring.informative=FALSE,
                one.step=one.step, no_cores=no_cores)
        run.fun(M=M, n=n, competing.risk=TRUE, tau=tau, setting=setting, censoring.informative=FALSE,
                one.step=one.step, no_cores=no_cores,
                misspecify.outcome=TRUE)

    }

    if (TRUE) {
        #-- sl for all parts
        run.fun(M=M, n=n, competing.risk=TRUE, tau=tau, setting=setting, censoring.informative=FALSE,
                one.step=one.step, no_cores=no_cores,
                fit.outcome="sl",
                fit.cens="sl",
                fit.cr="sl")



    }
}
