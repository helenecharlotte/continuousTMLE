
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
    no_cores <- 15
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
#source("/home/ifsv/jhl781/research/phd/berkeley/backups-continuousTMLE/nov-29/R/contmle.R")
source("./R/cox.super.learner.R")
source("./R/poisson.hal.R")
source("./R/poisson.hal.sl.R")

#-------------------------------------------------------------------------------------------#
## run: competing risks
#-------------------------------------------------------------------------------------------#

M <- 50#563#562#500#60#100#600#550#250
n <- 1000

setting <- 2
tau <- 0.5#c(0.2, 0.35, 0.5)
one.step <- FALSE


if (FALSE) {

    run.fun(M=1, n=n, competing.risk=TRUE, tau=tau, setting=setting,
            get.truth=TRUE, no_cores=no_cores)

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

if (FALSE) {
#-- sl for all parts
run.fun(M=M, n=n, competing.risk=TRUE, tau=tau, setting=setting, censoring.informative=TRUE,
        one.step=one.step, no_cores=no_cores,
        fit.outcome="sl",
        fit.cens="sl",
        fit.cr="sl")



#----- uninformative censoring ------#

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

