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
## set parameters
#-------------------------------------------------------------------------------------------#

betaA <- 0.1
betaL <- 0.3
nu <- 1.3#1.1#0.9#1.2
eta <- 1#1.2#1#4/sqrt(2)*(1/8)
tau <- 1#5

#-------------------------------------------------------------------------------------------#
## true value of target parameter
#-------------------------------------------------------------------------------------------#

(psi0 <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta, intervention.A=1, tau=tau) -
     sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta, intervention.A=0, tau=tau))

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
## generate data
#-------------------------------------------------------------------------------------------#

dt <- sim.data(1000, betaA=betaA, betaL=betaL, nu=nu, eta=eta)

dt[, summary(time)]

covars <- c("A", "L")

#-------------------------------------------------------------------------------------------#
## fit target parameter with coxph (g-formula)
#-------------------------------------------------------------------------------------------#

(psi.fit.cox <- fit.coxph(dt, tau=tau))

#-------------------------------------------------------------------------------------------#
## fit poisson hal for outcome intensity
#-------------------------------------------------------------------------------------------#

dt.hal <- fit.poisson.hal(dt, type=c("equal_range"), nbins=40, verbose=TRUE, lambda.cv=0.05,
                          covars=covars, tau=tau,
                          browse=FALSE)

#-------------------------------------------------------------------------------------------#
## initial estimation of target parameter
#-------------------------------------------------------------------------------------------#

estimate.target(dt.hal)

#-------------------------------------------------------------------------------------------#
## compare hal estimation of intensities to truth
#-------------------------------------------------------------------------------------------#

dt.hal[, true.lambda:=true.lambda(time, A, L, betaA=betaA, betaL=betaL, nu=nu, eta=eta)]
dt.hal[, true.Lambda:=true.Lambda(time, A, L, betaA=betaA, betaL=betaL, nu=nu, eta=eta)]
dt.hal[, true.density:=true.lambda*exp(-true.Lambda)]

AL <- unique(dt.hal[, c("A", "L"), with=FALSE])[sample(1:nrow(unique(dt.hal[, c("A", "L"), with=FALSE])), 4), c("A", "L"), with=FALSE]

par(mfrow=c(3,4))
#-- lambda;
plot(dt.hal[time<=tau  & A==AL[1,A] & L==AL[1,L], time], dt.hal[time<=tau  & A==AL[1,A] & L==AL[1,L], true.lambda], col="red", type="l")
lines(dt.hal[time<=tau  & A==AL[1,A] & L==AL[1,L], time], dt.hal[time<=tau  & A==AL[1,A] & L==AL[1,L], fit.lambda])
plot(dt.hal[time<=tau  & A==AL[2,A] & L==AL[2,L], time], dt.hal[time<=tau  & A==AL[2,A] & L==AL[2,L], true.lambda], col="red", type="l")
lines(dt.hal[time<=tau  & A==AL[2,A] & L==AL[2,L], time], dt.hal[time<=tau  & A==AL[2,A] & L==AL[2,L], fit.lambda])
plot(dt.hal[time<=tau  & A==AL[3,A] & L==AL[3,L], time], dt.hal[time<=tau  & A==AL[3,A] & L==AL[3,L], true.lambda], col="red", type="l")
lines(dt.hal[time<=tau  & A==AL[3,A] & L==AL[3,L], time], dt.hal[time<=tau  & A==AL[3,A] & L==AL[3,L], fit.lambda])
plot(dt.hal[time<=tau  & A==AL[4,A] & L==AL[4,L], time], dt.hal[time<=tau  & A==AL[4,A] & L==AL[4,L], true.lambda], col="red", type="l")
lines(dt.hal[time<=tau  & A==AL[4,A] & L==AL[4,L], time], dt.hal[time<=tau  & A==AL[4,A] & L==AL[4,L], fit.lambda])
#-- Lambda;
plot(dt.hal[time<=tau  & A==AL[1,A] & L==AL[1,L], time], dt.hal[time<=tau  & A==AL[1,A] & L==AL[1,L], true.Lambda], col="red", type="l")
lines(dt.hal[time<=tau  & A==AL[1,A] & L==AL[1,L], time], dt.hal[time<=tau  & A==AL[1,A] & L==AL[1,L], fit.Lambda])
plot(dt.hal[time<=tau  & A==AL[2,A] & L==AL[2,L], time], dt.hal[time<=tau  & A==AL[2,A] & L==AL[2,L], true.Lambda], col="red", type="l")
lines(dt.hal[time<=tau  & A==AL[2,A] & L==AL[2,L], time], dt.hal[time<=tau  & A==AL[2,A] & L==AL[2,L], fit.Lambda])
plot(dt.hal[time<=tau  & A==AL[3,A] & L==AL[3,L], time], dt.hal[time<=tau  & A==AL[3,A] & L==AL[3,L], true.Lambda], col="red", type="l")
lines(dt.hal[time<=tau  & A==AL[3,A] & L==AL[3,L], time], dt.hal[time<=tau  & A==AL[3,A] & L==AL[3,L], fit.Lambda])
plot(dt.hal[time<=tau  & A==AL[4,A] & L==AL[4,L], time], dt.hal[time<=tau  & A==AL[4,A] & L==AL[4,L], true.Lambda], col="red", type="l")
lines(dt.hal[time<=tau  & A==AL[4,A] & L==AL[4,L], time], dt.hal[time<=tau  & A==AL[4,A] & L==AL[4,L], fit.Lambda])
#-- density;
plot(dt.hal[time<=tau  & A==AL[1,A] & L==AL[1,L], time], dt.hal[time<=tau  & A==AL[1,A] & L==AL[1,L], true.density], col="red", type="l")
lines(dt.hal[time<=tau  & A==AL[1,A] & L==AL[1,L], time], dt.hal[time<=tau  & A==AL[1,A] & L==AL[1,L], fit.density])
plot(dt.hal[time<=tau  & A==AL[2,A] & L==AL[2,L], time], dt.hal[time<=tau  & A==AL[2,A] & L==AL[2,L], true.density], col="red", type="l")
lines(dt.hal[time<=tau  & A==AL[2,A] & L==AL[2,L], time], dt.hal[time<=tau  & A==AL[2,A] & L==AL[2,L], fit.density])
plot(dt.hal[time<=tau  & A==AL[3,A] & L==AL[3,L], time], dt.hal[time<=tau  & A==AL[3,A] & L==AL[3,L], true.density], col="red", type="l")
lines(dt.hal[time<=tau  & A==AL[3,A] & L==AL[3,L], time], dt.hal[time<=tau  & A==AL[3,A] & L==AL[3,L], fit.density])
plot(dt.hal[time<=tau  & A==AL[4,A] & L==AL[4,L], time], dt.hal[time<=tau  & A==AL[4,A] & L==AL[4,L], true.density], col="red", type="l")
lines(dt.hal[time<=tau  & A==AL[4,A] & L==AL[4,L], time], dt.hal[time<=tau  & A==AL[4,A] & L==AL[4,L], fit.density])
        
#-------------------------------------------------------------------------------------------#
## fit poisson hal
#-------------------------------------------------------------------------------------------#

#-- for censoring process; 
dt.hal.cens <- fit.poisson.hal(dt, type=c("equal_range"), which.delta=0, nbins=40, verbose=TRUE, lambda.cv=0.05)
setnames(dt.hal.cens, c("fit.lambda", "fit.Lambda"), paste0(c("fit.lambda", "fit.Lambda"), ".cens"))

dt.hal <- merge(dt.hal, dt.hal.cens[, c("tint", covars, paste0(c("fit.lambda", "fit.Lambda"), ".cens")), with=FALSE], by=c("tint", covars))

dt.hal[, fit.Surv.cens:=exp(-fit.Lambda.cens)]

#-------------------------------------------------------------------------------------------#
## fit treatment propensity
#-------------------------------------------------------------------------------------------#

fit.A <- glm(A ~ L, family=binomial(), data=dt)

dt.hal[, fit.pi:=predict(fit.A, newdata=dt.hal, type="response")]

#-------------------------------------------------------------------------------------------#
## targeting steps
#-------------------------------------------------------------------------------------------#

log.linear.targeting(dt.hal)
estimate.target(dt.hal, iteration=1)
log.linear.targeting(dt.hal, iteration=2)
estimate.target(dt.hal, iteration=2)
log.linear.targeting(dt.hal, iteration=3)
estimate.target(dt.hal, iteration=3)
log.linear.targeting(dt.hal, iteration=4)
estimate.target(dt.hal, iteration=4)


























if (FALSE) {

    #-------------------------------------------------------------------------------------------#
    ## try to do targeting
    #-------------------------------------------------------------------------------------------#


    #dt.hal[, summary(time)]
    #dt.hal[RT>0, summary(time)]
    #dt.hal[RT>0, summary(time)]
    #dt[delta==1, summary(time)]
    #dt.hal <- dt.hal[time <= tau]

    #dt.hal.A1 <- dt.hal[RT>0][A==1]
    #dt.hal.A0 <- dt.hal[RT>0][A==0]

    #setnames(dt.hal.A1, c("fit.lambda", "fit.Lambda"), paste0(c("fit.lambda", "fit.Lambda"), ".A1"))
    #setnames(dt.hal.A0, c("fit.lambda", "fit.Lambda"), paste0(c("fit.lambda", "fit.Lambda"), ".A0"))

    #dt.tst <- merge(dt.hal.A1[, c("tint", "L", paste0(c("fit.lambda", "fit.Lambda"), ".A1"))],
    #                dt.hal.A0[, c("tint", "L", paste0(c("fit.lambda", "fit.Lambda"), ".A0"))],
    #                by=c("tint", "L"))

    #dt.tst[, c("fit.lambda.A1", "fit.lambda.A0", "fit.Lambda.A1", "fit.Lambda.A0")]

    #dt.hal[, fit.Surv:=exp(-fit.Lambda)]

    #-- tjek: hvorfor NA? nu bare fjerne. 
    #dt.hal <- dt.hal[!is.na(fit.Surv)]

    #-- tjek: får jeg hentet den rigtige?  
    #dt.hal[, fit.Surv.tau:=fit.Surv[.N], by=covars]

    #-- ignore pi and Sc for now: 

    dt.hal[, clever.weight:=-( (A==1)/fit.pi - (A==0)/(1-fit.pi) ) * (time <= tau) / fit.Surv.cens]
    dt.hal[, clever.covar:=fit.Surv.tau/fit.Surv]
    dt.hal[, tmle.covar:=clever.weight*clever.covar]

    fit.tmle <- glm(as.formula(paste0("D ~ ", "tmle.covar-1",
                                      "+offset(log(RT*fit.lambda))")),
                    family=poisson(link=log),
                    data=dt.hal[RT>0])
    fit.tmle

    dt.hal[RT>0,
           fit.lambda.1:=exp(predict(fit.tmle, dt.hal[RT>0],
                                     newoffset=dt.hal[RT>0, log(RT*clever.weight)]))/RT]

    dt.hal[, c("fit.lambda", "fit.lambda.1")]

    #-- updated cumulative hazard
    dt.hal[, fit.Lambda.1:=cumsum(diff*fit.lambda.1), by=covars]

    dt.hal[, c("fit.Lambda", "fit.Lambda.1")]

    dt.hal[, fit.Surv.1:=exp(-fit.Lambda.1)]

    dt.hal[, c("fit.Surv", "fit.Surv.1")]

    dt.hal[, fit.Surv.1.tau:=fit.Surv.1[.N], by=covars]

    #--- estimation of target parameter?



    estimate.target(iteration=0)
    estimate.target(iteration=1)


    #--- try again

    dt.hal[, clever.covar.1:=fit.Surv.1.tau/fit.Surv.1]

    dt.hal[, tmle.covar.1:=clever.weight*clever.covar.1]

    fit.tmle <- glm(as.formula(paste0("D ~ ", "tmle.covar.1-1",
                                      "+offset(log(RT*fit.lambda.1))")),
                    family=poisson(link=log),
                    data=dt.hal[RT>0])
    fit.tmle

    dt.hal[RT>0,
           fit.lambda.2:=exp(predict(fit.tmle, dt.hal[RT>0],
                                     newoffset=dt.hal[RT>0, log(RT)]))/RT]

    dt.hal[, c("fit.lambda", "fit.lambda.1", "fit.lambda.2")]

    #-- updated cumulative hazard
    dt.hal[, fit.Lambda.2:=cumsum(diff*fit.lambda.2), by=covars]

    dt.hal[, c("fit.Lambda", "fit.Lambda.1", "fit.Lambda.2")]

    dt.hal[, fit.Surv.2:=exp(-fit.Lambda.2)]

    dt.hal[, c("fit.Surv", "fit.Surv.1", "fit.Surv.2")]

    dt.hal[, fit.Surv.2.tau:=fit.Surv.2[.N], by=covars]

    estimate.target(iteration=0)
    estimate.target(iteration=1)
    estimate.target(iteration=2)



    #--- try again, again

    dt.hal[, clever.covar.2:=fit.Surv.2.tau/fit.Surv.2]

    dt.hal[, tmle.covar.2:=clever.weight*clever.covar.2]

    fit.tmle <- glm(as.formula(paste0("D ~ ", "tmle.covar.2-1",
                                      "+offset(log(RT*fit.lambda.2))")),
                    family=poisson(link=log),
                    data=dt.hal[RT>0])
    fit.tmle

    dt.hal[RT>0,
           fit.lambda.3:=exp(predict(fit.tmle, dt.hal[RT>0],
                                     newoffset=dt.hal[RT>0, log(RT)]))/RT]

    dt.hal[, c("fit.lambda", "fit.lambda.1", "fit.lambda.2", "fit.lambda.3")]

    #-- updated cumulative hazard
    dt.hal[, fit.Lambda.3:=cumsum(diff*fit.lambda.3), by=covars]

    dt.hal[, c("fit.Lambda", "fit.Lambda.1", "fit.Lambda.2", "fit.Lambda.3")]

    dt.hal[, fit.Surv.3:=exp(-fit.Lambda.3)]

    dt.hal[, c("fit.Surv", "fit.Surv.1", "fit.Surv.2", "fit.Surv.3")]

    dt.hal[, fit.Surv.3.tau:=fit.Surv.3[.N], by=covars]

    estimate.target(iteration=0)
    estimate.target(iteration=1)
    estimate.target(iteration=2)
    estimate.target(iteration=3)

    psi0


    #--- hvad med cox?

    library(survival)
    fit.cox <- coxph(Surv(time, delta) ~ A+L, data=dt)

    dt.1 <- copy(dt)
    dt.0 <- copy(dt)
    dt.1[, A:=1]
    dt.1[, time:=tau]
    dt.0[, A:=0]
    dt.0[, time:=tau]

    dt[, surv.cox.A1:=exp(-predict(fit.cox, newdata=dt.1, type="expected"))]
    dt[, surv.cox.A0:=exp(-predict(fit.cox, newdata=dt.0, type="expected"))]

    dt[, mean(surv.cox.A0-surv.cox.A1)]


}
