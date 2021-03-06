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

source("./R/sim-data.R")
source("./R/est-fun2.R")
source("./R/fit-density-by-hazard.R")
source("./simulation-stochastic/coverage-fun.R")
source("./simulation-stochastic/repeat-fun.R")
source("./simulation-stochastic/compute-true.R")

#-------------------------------------------------------------------------------------------#
## set parameters
#-------------------------------------------------------------------------------------------#

K <- 15#5#5#100#100#100#100#80#100
run.ltmle <- FALSE##TRUE#FALSE
run.ctmle <- TRUE#TRUE#FALSE#FALSE
run.ctmle2 <- FALSE#FALSE#FALSE
compute.true.eic <- FALSE
compute.true.psi <- FALSE
misspecify.Q <- FALSE
M <- 600#250#600#600#250#350#500#551
n <- 1000
only.A0 <- FALSE
rescue.itt <- FALSE
estimate.Gstar <- TRUE

#-------------------------------------------------------------------------------------------#
## true values (outputs to file)
#-------------------------------------------------------------------------------------------#

compute.true(compute.true.psi=compute.true.psi, compute.true.eic=compute.true.eic)

#-------------------------------------------------------------------------------------------#
## repeat simulations (parallelize)
#-------------------------------------------------------------------------------------------#

if (system("echo $USER",intern=TRUE)%in%c("jhl781")){ 
    no_cores <- 18
} else {
    no_cores <- detectCores() - 1
}

registerDoParallel(no_cores)

out <- foreach(m=1:M, .errorhandling="pass"#, #.combine=list, .multicombine = TRUE
               ) %dopar% {
                   repeat.fun(m, K=K, n=n, estimate.Gstar=estimate.Gstar, 
                              run.ltmle=run.ltmle, run.ctmle=run.ctmle,
                              run.ctmle2=run.ctmle2,
                              rescue.itt=rescue.itt, 
                              misspecify.Q=misspecify.Q)
               }

stopImplicitCluster()

saveRDS(out,
        file=paste0("./simulation-stochastic/output/",
                    "outlist-est-rescue", ifelse(rescue.itt, "-itt", ""),
                    ifelse(run.ltmle, "-ltmle", ""),
                    ifelse(run.ctmle, "-ctmle", ""),
                    ifelse(run.ctmle2, "-ctmle2", ""),
                    ifelse(estimate.Gstar, "-estimateG", ""),
                    ifelse(n==1000, "", paste0("-n", n)), 
                    "-2020", "-K", K, ifelse(misspecify.Q, "-Q", ""), 
                    "-M", M, ".rds"))




