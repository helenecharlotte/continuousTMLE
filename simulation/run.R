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
library(stringr)
library(ltmle)
library(parallel)


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

source("./R/sim-data.R")
source("./R/est-fun.R")
source("./R/fit-density-by-hazard.R")
source("./simulation/coverage-fun.R")


#-------------------------------------------------------------------------------------------#
## repeat simulations
#-------------------------------------------------------------------------------------------#

K <- 100#100#100#100#100#80#100
run.ltmle <- FALSE##TRUE#FALSE
run.ctmle <- TRUE#FALSE#FALSE
compute.true.eic <- TRUE
misspecify.Q <- TRUE
only.A0 <- FALSE


(psi0.test.multi.M0 <- sim.data(1e6, seed=10011,
                                         only.A0=only.A0,
                                         intervention.A=function(L0, L.prev, A.prev, A1) cbind(logit(1)),
                                         K=K
                                         ))

(psi0.test.multi.M1 <- sim.data(1e6, seed=10011,
                                         only.A0=only.A0,
                                         intervention.A=function(L0, L.prev, A.prev, A1) cbind(logit(0)),
                                         K=K
                                         ))

psi0.test.multi.M1 - psi0.test.multi.M0

if (compute.true.eic) {

    N <- 1e5
    m <- 1
    dt <- sim.data(N, seed=10011+m, censoring=TRUE,
                            only.A0=only.A0,
                            browse=FALSE,
                            K=K)

    true.eic.0 <-  suppressMessages(est.fun(copy(dt), censoring=TRUE,
                                            targeting=1, 
                                            smooth.initial=TRUE,
                                            browse9=FALSE,
                                            compute.true.eic=TRUE,
                                            intervention.A0=function(L0, A0) logit(1*(A0==0)),
                                            intervention.A=function(L0, A0, L.prev, A.prev, A) logit(1*(A==0)),
                                            browse0=FALSE, misspecify.Q=misspecify.Q))
    
    true.eic.1 <- suppressMessages(est.fun(copy(dt), censoring=TRUE,
                                           targeting=1, 
                                           smooth.initial=TRUE,
                                           browse9=FALSE,
                                           compute.true.eic=TRUE,
                                           intervention.A0=function(L0, A0) logit(1*(A0==1)),
                                           intervention.A=function(L0, A0, L.prev, A.prev, A) logit(1*(A==1)),
                                           browse0=FALSE, misspecify.Q=misspecify.Q))

    message("--------------------------------")
    print(true.eic0 <- true.eic.0[[1]][3] * sqrt(N) / sqrt(1000))
    print(true.eic1 <- true.eic.1[[1]][3] * sqrt(N) / sqrt(1000))

    saveRDS(true.eic0,
            file=paste0("./simulations/output/",
                        "outlist-est-maximum-likelihood-multi-true-sd-0-2020",
                        "-K", K, ifelse(misspecify.Q, "-Q", ""),
                        "-M", M, ".rds"))

    saveRDS(true.eic1,
            file=paste0("./simulations/output/",
                        "outlist-est-maximum-likelihood-multi-true-sd-1-2020",
                        "-K", K, ifelse(misspecify.Q, "-Q", ""),
                        "-M", M, ".rds"))


    est.list.test.multi.M0 <- list(true=c(psi0.test.multi.M0, true.eic0))
    est.list.test.multi.M1 <- list(true=c(psi0.test.multi.M1, true.eic1))
    est.list.test.multi.target2.M0 <- list(true=c(psi0.test.multi.M0, true.eic0))
    est.list.test.multi.target2.M1 <- list(true=c(psi0.test.multi.M1, true.eic1))

    ltmle.list.0 <- list(true=c(psi0.test.multi.M0, true.eic0))
    ltmle.list.1 <- list(true=c(psi0.test.multi.M1, true.eic1))
    
} else {

    est.list.test.multi.M0 <- list(true=psi0.test.multi.M0)
    est.list.test.multi.M1 <- list(true=psi0.test.multi.M1)
    est.list.test.multi.target2.M0 <- list(true=psi0.test.multi.M0)
    est.list.test.multi.target2.M1 <- list(true=psi0.test.multi.M1)

    ltmle.list.0 <- list(true=psi0.test.multi.M0)
    ltmle.list.1 <- list(true=psi0.test.multi.M1)
}


M <- 1#302#502#300#252#300#1#500#400#400#800#10#400
for (m in 1:M) { 

    dt <- sim.data(1000, seed=10011+m, censoring=TRUE,
                   only.A0=only.A0,
                   browse=FALSE,
                   K=K)

    if (run.ctmle) {
        
        (est.list.test.multi.M0[[m+1]] <- est.fun(copy(dt), censoring=TRUE,
                                                  targeting=1, 
                                                  smooth.initial=TRUE,
                                                  browse3=FALSE,
                                                  intervention.A0=function(L0, A0) logit(1*(A0==0)),
                                                  intervention.A=function(L0, A0, L.prev, A.prev, A) logit(1*(A==0)),
                                                  browse0=FALSE, misspecify.Q=misspecify.Q))

        (est.list.test.multi.M1[[m+1]] <- est.fun(copy(dt), censoring=TRUE,
                                                  targeting=1, 
                                                  smooth.initial=TRUE,
                                                  browse9=FALSE,
                                                  intervention.A0=function(L0, A0) logit(1*(A0==1)),
                                                  intervention.A=function(L0, A0, L.prev, A.prev, A) logit(1*(A==1)),
                                                  browse0=FALSE, misspecify.Q=misspecify.Q))

        (est.list.test.multi.target2.M0[[m+1]] <- est.fun(copy(dt), censoring=TRUE,
                                                          targeting=2, 
                                                          smooth.initial=TRUE,
                                                          browse5=FALSE,
                                                          intervention.A0=function(L0, A0) logit(1*(A0==0)),
                                                          intervention.A=function(L0, A0, L.prev, A.prev, A) logit(1*(A==0)),
                                                          browse0=FALSE, misspecify.Q=misspecify.Q))

        (est.list.test.multi.target2.M1[[m+1]] <- est.fun(copy(dt), censoring=TRUE,
                                                          targeting=2, 
                                                          smooth.initial=TRUE,
                                                          browse9=FALSE,
                                                          intervention.A0=function(L0, A0) logit(1*(A0==1)),
                                                          intervention.A=function(L0, A0, L.prev, A.prev, A) logit(1*(A==1)),
                                                          browse0=FALSE, misspecify.Q=misspecify.Q))

        saveRDS(est.list.test.multi.M1,
                file=paste0("./simulations/output/",
                            "outlist-est-maximum-likelihood-multi-M1-2020",
                            "-K", K, ifelse(only.A0, "-A0", ""), ifelse(misspecify.Q, "-Q", ""), 
                            "-M", M, ".rds"))

        saveRDS(est.list.test.multi.M0,
                file=paste0("./simulations/output/",
                            "outlist-est-maximum-likelihood-multi-M0-2020",
                            "-K", K, ifelse(only.A0, "-A0", ""), ifelse(misspecify.Q, "-Q", ""),
                            "-M", M, ".rds"))

        saveRDS(est.list.test.multi.target2.M1,
                file=paste0("./simulations/output/",
                            "outlist-est-mininum-loss-multi-M1-2020",
                            "-K", K, ifelse(only.A0, "-A0", ""), ifelse(misspecify.Q, "-Q", ""),
                            "-M", M, ".rds"))

        saveRDS(est.list.test.multi.target2.M0,
                file=paste0("./simulations/output/",
                            "outlist-est-minimum-loss-multi-M0-2020",
                            "-K", K, ifelse(only.A0, "-A0", ""), ifelse(misspecify.Q, "-Q", ""),
                            "-M", M, ".rds"))
    }
    
    if (run.ltmle) {

        df.wide <- as.data.frame(cbind(dt[, "L0"],
                                       dN.A0=rep(1, nrow(dt)),
                                       dt[, -c("L0", "id")]))
        
        col.order <- names(df.wide)[-1]

        Anodes <- col.order[substr(col.order, 1, 1)=="A"]
        Ynodes <- col.order[substr(col.order, 1, 1)=="Y"]
        Cnodes <- col.order[substr(col.order, 1, 1)=="C"]
        Nnodes <- col.order[substr(col.order, 1, 4)=="dN.A"]
        Lnodes <- setdiff(col.order, c("L0", Anodes, Ynodes, Cnodes))

        for (ii in Cnodes) {
            df.wide[, names(df.wide)==ii] <- BinaryToCensoring(is.censored=df.wide[, names(df.wide)==ii])
        }

        if (TRUE) {
            abar0 <- (df.wide[, Nnodes] == 1) * 0 + (df.wide[, Nnodes] == 0) * df.wide[, Anodes]
            abar1 <- (df.wide[, Nnodes] == 1) * 1 + (df.wide[, Nnodes] == 0) * df.wide[, Anodes]
            #abar <- (df.wide[, Nnodes] == 1) * 0 + (df.wide[, Nnodes] == 0) * 0
    
            #--- specify intervention step 2: when dA(dt)=1 : 
            det.g.fun <- function(data, current.node, nodes) {
                if (substr(names(data)[current.node], 1, 1)=="C") {
                    if (FALSE) {
                        C.node.index <- current.node
                        observed.C <- data[, C.node.index]
                        is.deterministic <- observed.C=="censored" | is.na(observed.C)
                        #browser()
                        prob1 <- observed.C[is.deterministic]#observed.C[is.deterministic]
                        # if ( names(data)[current.node] == "C2") browser()
                        return(list(is.deterministic = is.deterministic, prob1 = prob1))
                    }
                    #return(NULL)
                } else {                
                    N.nodes.index <- grep("dN.A", names(data))
                    N.node.index <- max(N.nodes.index[N.nodes.index < current.node]) 
                    A.node.index <- current.node
                    N <- data[, N.node.index]
                    observed.A <- data[, A.node.index]
                    is.deterministic <- N == 0 | is.na(N)
                    prob1 <- observed.A[is.deterministic]
                    #if ( names(data)[current.node] == "A2") browser()
                    return(list(is.deterministic = is.deterministic, prob1 = prob1))
                }
            }

            print("running ltmle with deterministic arule")

            r0 <- suppressMessages(try(ltmle(df.wide,#[,!names(df.wide)%in%Cnodes],##[, col.order]
                                             Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes,
                                             abar=as.matrix(abar0),
                                             #abar=rep(0, length(Anodes)), 
                                             Ynodes=Ynodes, survivalOutcome=TRUE, 
                                             deterministic.g.function=det.g.fun,
                                             estimate.time=FALSE)))

            r1 <- suppressMessages(try(ltmle(df.wide,#[,!names(df.wide)%in%Cnodes],##[, col.order]
                                             Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes,
                                             abar=as.matrix(abar1),
                                             #abar=rep(0, length(Anodes)), 
                                             Ynodes=Ynodes, survivalOutcome=TRUE, 
                                             deterministic.g.function=det.g.fun,
                                             estimate.time=FALSE)))
        } else {

            r0 <- suppressMessages(try(ltmle(df.wide,#[,!names(df.wide)%in%Cnodes],##[, col.order]
                                             Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes,
                                             #abar=as.matrix(abar),
                                             abar=rep(0, length(Anodes)), 
                                             Ynodes=Ynodes, survivalOutcome=TRUE, 
                                             #deterministic.g.function=det.g.fun,
                                             estimate.time=FALSE)))

            r1 <- suppressMessages(try(ltmle(df.wide,#[,!names(df.wide)%in%Cnodes],##[, col.order]
                                             Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes,
                                             #abar=as.matrix(abar),
                                             abar=rep(1, length(Anodes)), 
                                             Ynodes=Ynodes, survivalOutcome=TRUE, 
                                             #deterministic.g.function=det.g.fun,
                                             estimate.time=FALSE)))

        }
        
        if (is(r0, "try-error")) {
            ltmle.list.0[[m+1]] <- "ERROR"
        } else {
            ltmle.list.0[[m+1]] <- list(est=r0$estimates["tmle"],
                                        sd=summary(r0)$treatment$std.dev)
        }

        if (is(r1, "try-error")) {
            ltmle.list.1[[m+1]] <- "ERROR"
        } else {
            ltmle.list.1[[m+1]] <- list(est=r1$estimates["tmle"],
                                        sd=summary(r1)$treatment$std.dev)
        }

        saveRDS(ltmle.list.0,
                file=paste0("./simulations/output/",
                            "outlist-est-ltmle-0-2020",
                            "-K", K, ifelse(only.A0, "-A0", ""), ifelse(misspecify.Q, "-Q", ""),
                            "-M", M, ".rds"))

        saveRDS(ltmle.list.1,
                file=paste0("./simulations/output/",
                            "outlist-est-ltmle-1-2020",
                            "-K", K, ifelse(only.A0, "-A0", ""), ifelse(misspecify.Q, "-Q", ""),
                            "-M", M, ".rds"))
    }    
  
}




