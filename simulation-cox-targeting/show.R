### show.R --- 
#----------------------------------------------------------------------
## Author: Helene Charlotte Rytgaard
## Created: May, 2020 
#----------------------------------------------------------------------
## 
### Commentary:
## Look at output from TMLE estimation.
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
## parameters of interest
#-------------------------------------------------------------------------------------------#

M <- 250

#-------------------------------------------------------------------------------------------#
## get output from simulations
#-------------------------------------------------------------------------------------------#

length(out <- readRDS(file=paste0("./simulation-cox-targeting/output/",
                                  "outlist-est",
                                  ifelse(interaction.AL, "-interactionAL", ""),
                                  "-M", M, ".rds")))


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
interaction.AL <- TRUE#FALSE
misspecify.Y <- TRUE
interaction.Atime <- TRUE
randomize.A <- TRUE
censoring.informative <- TRUE
fit.km <- TRUE

#--  which effect estimated? 
a <- 1

if (interaction.Atime) betaA <- -0.65
t0 <- 0.9
tau <- 1.2

#-------------------------------------------------------------------------------------------#
## get true value
#-------------------------------------------------------------------------------------------#

(psi0 <- readRDS(file=paste0("./simulation-cox-targeting/output/",
                             "outlist-psi0",
                             ifelse(interaction.AL, "-interactionAL", ""),
                             ifelse(interaction.Atime, "-interactionAtime", ""),
                             ".rds")))

psi0.A1 <- psi0["psi0.A1"]
psi0.A0 <- psi0["psi0.A0"]
psi0 <- psi0["psi0"]

#-------------------------------------------------------------------------------------------#
## look at results of interest
#-------------------------------------------------------------------------------------------#

M <- 350#210#1000#250#402#401#600#600#500#300#600#600#550#1000#1000#450#500#500#200#250#100#500
a <- 1
misspecify.Y <- TRUE#TRUE
interaction.AL <- FALSE#TRUE#FALSE#TRUE#FALSE
interaction.Atime <- TRUE#TRUE#TRUE#TRUE
censoring.informative <- TRUE#TRUE#FALSE
randomize.A <- TRUE
centered <- TRUE
length(out <- readRDS(file=paste0("./simulation-cox-targeting/output/",
                                  "outlist-est",
                                  ifelse(a==0, "-A=0", "-A=1"),
                                  ifelse(interaction.AL, "-interactionAL", ""),
                                  ifelse(interaction.Atime, "-interactionAtime", ""),
                                  ifelse(randomize.A, "-randomizeA", ""),
                                  ifelse(censoring.informative, "", "-almostUninformativeCensoring"),
                                  ifelse(censoring.high, "-highLevelOfCensoring", ""),
                                  ifelse(misspecify.Y, "-misspecifyY", ""),
                                  ifelse(centered, "-centered", ""),
                                  "-M", M, ".rds")))

#-- initial and tmle1 estimator of psi0.A1: 
mean(unlist(lapply(out, function(x) x[[1]][1])))
mean(unlist(lapply(out, function(x) x[[length(x)]][1])))
if (a==1) (psi0.true <- psi0.A1) else (psi0.true <- psi0.A0)

par(mfrow=c(1,2))
hist(unlist(lapply(out, function(x) x[[1]][1])), main="init", xlab="psi.hat", ylab="")
abline(v=psi0.true, col="red")
hist(unlist(lapply(out, function(x) x[[length(x)]][1])), main="tmle (max. 5 iterations)", xlab="psi.hat", ylab="")
abline(v=psi0.true, col="red")
par(mfrow=c(1,1))

#mean(unlist(lapply(out, function(x) x[[2]][2])))
mean(unlist(lapply(out, function(x) x[[length(x)]][2])))

mse(unlist(lapply(out, function(x) x[[1]][1])))
mse(unlist(lapply(out, function(x) x[[length(x)]][1])))
sd2 <- sd(unlist(lapply(out, function(x) x[[length(x)]][1])))

#-- check coverage: 
print(paste0("cov=", (coverage <- mean(unlist(lapply(out, function(x, psi0=psi0.true) {
    c(x[[length(x)]][1]-1.96*x[[length(x)]][2]<=psi0 &
      x[[length(x)]][1]+1.96*x[[length(x)]][2]>=psi0)
}))))))

#-- check "oracle" coverage: 
print(paste0("oracle cov=", (coverage2 <- mean(unlist(lapply(out, function(x, psi0=psi0.true) {
    c(x[[length(x)]][1]-1.96*sd2<=psi0 &
      x[[length(x)]][1]+1.96*sd2>=psi0)
}))))))


#--- look at other;

#-- naive hr:
mean(exp(unlist(lapply(out, function(x) x[[1]]["hr.A"]))))
mean(unlist(lapply(out, function(x) x[[1]]["hr.pval"]))<=0.05)

#-- unadjusted km:
c("km.fit"=mean(unlist(lapply(out, function(x) x[[1]]["km.est"]))), "truth"=psi0.true)
c("km.sd"=sd(unlist(lapply(out, function(x) x[[1]]["km.est"]))), "tmle.sd"=sd2)
c("km.se.est"=mean(unlist(lapply(out, function(x) x[[1]]["km.se"]))), "tmle.se.est"=mean(unlist(lapply(out, function(x) x[[length(x)]][2]))))

c("km.bias"=mean(unlist(lapply(out, function(x) x[[1]]["km.est"])))-psi0.true,
  "tmle.bias"=mean(unlist(lapply(out, function(x) x[[length(x)]][1])))-psi0.true)
c("km.mse"=mse(unlist(lapply(out, function(x) x[[1]]["km.est"]))),
  "tmle.mse"=mse(unlist(lapply(out, function(x) x[[length(x)]][1]))),
  "rel.mse"=tmle.mse/km.mse)

#-- model check; estimation of time-varying coef:
if (!misspecify.Y & interaction.Atime) {
    print(c(mean(unlist(lapply(out, function(x) x[[1]][5]))), betaA))
    print(c(mean(unlist(lapply(out, function(x) x[[1]][6]))), -0.45*betaA))
    print(c(mean(unlist(lapply(out, function(x) x[[1]][7]))), -betaL))
    print(c(mean(unlist(lapply(out, function(x) x[[1]][8]))), -1.2))
    print(c(mean(unlist(lapply(out, function(x) x[[1]][9]))), 0.8))
    #print(c(mean(unlist(lapply(out, function(x) x[[1]][10]))), -0.3))
}


#-------------------------------------------------------------------------------------------#
## look at results combined
#-------------------------------------------------------------------------------------------#

M <- 350#1000#600#500#250#600#300 #600#550#500#1000#450#500#500#200#250#100#500
misspecify.Y <- FALSE
censoring.informative <- TRUE#FALSE#TRUE#FALSE#FALSE
interaction.AL <- FALSE#FALSE
interaction.Atime <- TRUE
length(out1 <- readRDS(file=paste0("./simulation-cox-targeting/output/",
                                   "outlist-est",
                                   ifelse(1==0, "-A=0", "-A=1"),
                                   ifelse(interaction.AL, "-interactionAL", ""),
                                   ifelse(interaction.Atime, "-interactionAtime", ""),
                                   ifelse(randomize.A, "-randomizeA", ""),
                                   ifelse(censoring.informative, "", "-almostUninformativeCensoring"),
                                   ifelse(censoring.high, "-highLevelOfCensoring", ""),
                                   ifelse(misspecify.Y, "-misspecifyY", ""),
                                   ifelse(centered, "-centered", ""),
                                   "-M", M, ".rds")))
length(out0 <- readRDS(file=paste0("./simulation-cox-targeting/output/",
                                   "outlist-est",
                                   ifelse(0==0, "-A=0", "-A=1"),
                                   ifelse(interaction.AL, "-interactionAL", ""),
                                   ifelse(interaction.Atime, "-interactionAtime", ""),
                                   ifelse(randomize.A, "-randomizeA", ""),
                                   ifelse(censoring.informative, "", "-almostUninformativeCensoring"),
                                   ifelse(censoring.high, "-highLevelOfCensoring", ""),
                                   ifelse(misspecify.Y, "-misspecifyY", ""),
                                   ifelse(centered, "-centered", ""),
                                   "-M", M, ".rds")))

psi0

mean(unlist(lapply(out1, function(x) x[[1]][1])))-mean(unlist(lapply(out0, function(x) x[[1]][1])))
mean(unlist(lapply(out1, function(x) x[[length(x)]][1])))-mean(unlist(lapply(out0, function(x) x[[length(x)]][1])))

par(mfrow=c(1,2))
hist(unlist(lapply(out1, function(x) x[[1]][1]))-
     unlist(lapply(out0, function(x) x[[1]][1])), main="init (estimate diff)", xlab="psi.hat", ylab="")
abline(v=psi0, col="red")
hist(unlist(lapply(out1, function(x) x[[length(x)]][1]))-
     unlist(lapply(out0, function(x) x[[length(x)]][1])), main="tmle (estimate diff)", xlab="psi.hat", ylab="")
abline(v=psi0, col="red")
par(mfrow=c(1,1))

#-- check coverage: 
print(paste0("cov=", (coverage <- mean(unlist(lapply(1:length(out1), function(kk, psi0.true=psi0) {
    x1 <- out1[[kk]]
    x0 <- out0[[kk]]
    c(x1[[length(x1)]][1]-x0[[length(x0)]][1]-1.96*sqrt(x1[[length(x1)]][2]^2+x0[[length(x0)]][2]^2)<=psi0.true &
      x1[[length(x1)]][1]-x0[[length(x0)]][1]+1.96*sqrt(x1[[length(x1)]][2]^2+x0[[length(x0)]][2]^2)>=psi0.true)
}))))))

sd0 <- sd(unlist(lapply(out0, function(x) x[[length(x)]][1])))
sd1 <- sd(unlist(lapply(out1, function(x) x[[length(x)]][1])))

#-- check oracle coverage: 
print(paste0("oracle cov=", (oracle.coverage <- mean(unlist(lapply(1:length(out1), function(kk, psi0.true=psi0) {
    x1 <- out1[[kk]]
    x0 <- out0[[kk]]
    c(x1[[length(x1)]][1]-x0[[length(x0)]][1]-1.96*sqrt(sd1^2+sd0^2)<=psi0.true &
      x1[[length(x1)]][1]-x0[[length(x0)]][1]+1.96*sqrt(sd1^2+sd0^2)>=psi0.true)
}))))))

tmle.est <- mean(unlist(lapply(1:length(out1), function(kk, psi0.true=psi0) {
    x1 <- out1[[kk]]
    x0 <- out0[[kk]]
    return(x1[[length(x1)]][1]-x0[[length(x0)]][1])
})))

tmle.se.est <- mean(unlist(lapply(1:length(out1), function(kk, psi0.true=psi0) {
    x1 <- out1[[kk]]
    x0 <- out0[[kk]]
    return(sqrt(x1[[length(x1)]][2]^2+x0[[length(x0)]][2]^2))
})))

tmle.sd <- sd(unlist(lapply(1:length(out1), function(kk, psi0.true=psi0) {
    x1 <- out1[[kk]]
    x0 <- out0[[kk]]
    return(x1[[length(x1)]][1]-x0[[length(x0)]][1])
})))

tmle.mse <- mse(unlist(lapply(1:length(out1), function(kk, psi0.true=psi0) {
    x1 <- out1[[kk]]
    x0 <- out0[[kk]]
    return(x1[[length(x1)]][1]-x0[[length(x0)]][1])
})))

#-- naive hr:
mean(exp(unlist(lapply(out1, function(x) x[[1]]["hr.A"]))))
mean(unlist(lapply(out1, function(x) x[[1]]["hr.pval"]))<=0.05)

#-- unadjusted km:

km.est <- mean(unlist(lapply(1:length(out1), function(kk, psi0.true=psi0) {
    x1 <- out1[[kk]]
    x0 <- out0[[kk]]
    return(x1[[1]]["km.est"]-x0[[1]]["km.est"])
})))

km.se.est <- mean(unlist(lapply(1:length(out1), function(kk, psi0.true=psi0) {
    x1 <- out1[[kk]]
    x0 <- out0[[kk]]
    return(sqrt(x1[[1]]["km.se"]^2+x0[[1]]["km.se"]^2))
})))

km.sd <- sd(unlist(lapply(1:length(out1), function(kk, psi0.true=psi0) {
    x1 <- out1[[kk]]
    x0 <- out0[[kk]]
    return(x1[[1]]["km.est"]-x0[[1]]["km.est"])
})))

km.mse <- sd(unlist(lapply(1:length(out1), function(kk, psi0.true=psi0) {
    x1 <- out1[[kk]]
    x0 <- out0[[kk]]
    return(x1[[1]]["km.est"]-x0[[1]]["km.est"])
})))

c("km.fit"=km.est, "tmle.fit"=tmle.est, "truth"=psi0)
c("km.bias"=km.est-psi0, "tmle.bias"=tmle.est-psi0)
c("km.sd"=km.sd, "tmle.sd"=tmle.sd)
c("km.mse"=km.mse, "tmle.mse"=tmle.mse, "rel.mse"=tmle.mse/km.mse)
c("km.se.est"=km.se.est, "tmle.se.est"=tmle.se.est)


#-------------------------------------------------------------------------------------------#
## true survival difference? 
#-------------------------------------------------------------------------------------------#

surv.list <- readRDS(file=paste0("./simulation-cox-targeting/output/",
                                 "outlist-psi0-survival-function",
                                 ifelse(interaction.AL, "-interactionAL", ""),
                                 ifelse(interaction.Atime, "-interactionAtime", ""),
                                 ".rds"))

library(ggplot2)
ggplot(data.table(t=as.numeric(surv.list[, 1]), true.diff=as.numeric(surv.list[, 2]))) +
    theme_minimal() + 
    geom_line(aes(x=t, y=true.diff)) + xlab("time") + ylab(expression(psi[0])) +
    ggtitle("True value of target parameter (risk difference)")

ggsave("./simulation-cox-targeting/output/plot-true-survival-difference.pdf", width=10, height=5.5)



surv.list.A1 <- cbind(readRDS(file=paste0("./simulation-cox-targeting/output/",
                                          "outlist-psi0-survival-A1-function",
                                          ifelse(interaction.AL, "-interactionAL", ""),
                                          ifelse(interaction.Atime, "-interactionAtime", ""),
                                          ".rds")), treatment="a=1")

surv.list.A0 <- cbind(readRDS(file=paste0("./simulation-cox-targeting/output/",
                                          "outlist-psi0-survival-A0-function",
                                          ifelse(interaction.AL, "-interactionAL", ""),
                                          ifelse(interaction.Atime, "-interactionAtime", ""),
                                          ".rds")), treatment="a=0")



surv.both <- data.table(t=c(unlist(surv.list.A1[, 1]), unlist(surv.list.A0[, 1])),
                        S=c(1-unlist(surv.list.A1[, 2]), 1-unlist(surv.list.A0[, 2])),
                        treatment=c(unlist(surv.list.A1[, 3]), unlist(surv.list.A0[, 3])))

library(ggplot2)
ggplot(surv.both) +
    theme_minimal() + theme_bw(base_size=15) +
    geom_line(aes(x=t, y=S, linetype=treatment)) +
    xlab("time") + ylab(expression(psi[0]^{"A=a"})) +
    ggtitle("True counterfactual survival curves")

ggsave("./simulation-cox-targeting/output/plot-counterfactual-survival-curves.pdf", width=10, height=5.5)






