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
## source relevant scripts
#-------------------------------------------------------------------------------------------#

source("./R/sim-data.R")
source("./R/est-fun.R")
source("./R/fit-density-by-hazard.R")
source("./simulation/set-wd.R")
source("./simulation/coverage-fun.R")
source("./simulation/repeat-fun.R")


#-------------------------------------------------------------------------------------------#
## set parameters
#-------------------------------------------------------------------------------------------#

K <- 40#100#100#100#100#80#100
run.ltmle <- TRUE##TRUE#FALSE
run.ctmle <- FALSE#FALSE#FALSE
run.ctmle2 <- FALSE#FALSE#FALSE
compute.true.eic <- FALSE
compute.true.psi <- FALSE
misspecify.Q <- FALSE
only.A0 <- FALSE
M <- 300

#-------------------------------------------------------------------------------------------#
## true values
#-------------------------------------------------------------------------------------------#

if (compute.true.psi) {
    print(psi0.test.multi.M0 <- sim.data(1e6, seed=10011,
                                         only.A0=only.A0,
                                         intervention.A=function(L0, L.prev, A.prev, A1) cbind(logit(1)),
                                         K=K
                                         ))

    print(psi0.test.multi.M1 <- sim.data(1e6, seed=10011,
                                         only.A0=only.A0,
                                         intervention.A=function(L0, L.prev, A.prev, A1) cbind(logit(0)),
                                         K=K
                                         ))


    saveRDS(psi0.test.multi.M0,
            file=paste0("./simulation/output/",
                        "outlist-est-true-0-2020",
                        "-K", K, ifelse(misspecify.Q, "-Q", ""),
                        "-M", M, ".rds"))

    saveRDS(psi0.test.multi.M1,
            file=paste0("./simulation/output/",
                        "outlist-est-true-1-2020",
                        "-K", K, ifelse(misspecify.Q, "-Q", ""),
                        "-M", M, ".rds"))

    print(psi0.test.multi.M1 - psi0.test.multi.M0)
}

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
            file=paste0("./simulation/output/",
                        "outlist-est-true-sd-0-2020",
                        "-K", K, ifelse(misspecify.Q, "-Q", ""),
                        "-M", M, ".rds"))

    saveRDS(true.eic1,
            file=paste0("./simulation/output/",
                        "outlist-est-true-sd-1-2020",
                        "-K", K, ifelse(misspecify.Q, "-Q", ""),
                        "-M", M, ".rds"))

    ## est.list.test.multi.M0 <- list(true=c(psi0.test.multi.M0, true.eic0))
    ## est.list.test.multi.M1 <- list(true=c(psi0.test.multi.M1, true.eic1))
    ## est.list.test.multi.target2.M0 <- list(true=c(psi0.test.multi.M0, true.eic0))
    ## est.list.test.multi.target2.M1 <- list(true=c(psi0.test.multi.M1, true.eic1))

    ## ltmle.list.0 <- list(true=c(psi0.test.multi.M0, true.eic0))
    ## ltmle.list.1 <- list(true=c(psi0.test.multi.M1, true.eic1))
    
} 



#-------------------------------------------------------------------------------------------#
## repeat simulations (parallelize
#-------------------------------------------------------------------------------------------#

if (system("echo $USER",intern=TRUE)%in%c("jhl781")){ 
    no_cores <- 25
} else {
    no_cores <- detectCores() - 1
}

registerDoParallel(no_cores)

out <- foreach(m=1:M, #.combine=list,
               .multicombine = TRUE) %dopar% {
                   repeat.fun(m, K=K,
                              only.A0=only.A0, run.ltmle=run.ltmle, run.ctmle=run.ctmle, run.ctmle2=run.ctmle2,
                              misspecify.Q=misspecify.Q)
               }

stopImplicitCluster()

saveRDS(out,
        file=paste0("./simulation/output/",
                    "outlist-est",
                    ifelse(run.ltmle, "-ltmle", ""),
                    ifelse(run.ctmle, "-ctmle", ""),
                    ifelse(run.ctmle2, "-ctmle2", ""), 
                    "-2020", "-K", K, ifelse(only.A0, "-A0", ""), ifelse(misspecify.Q, "-Q", ""), 
                    "-M", M, ".rds"))




