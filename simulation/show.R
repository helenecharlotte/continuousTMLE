#-------------------------------------------------------------------------------------------#
## parameters of interest
#-------------------------------------------------------------------------------------------#

K <- 100#100#100#100#100#80#100
run.ltmle <- FALSE##TRUE#FALSE
run.ctmle <- FALSE#FALSE#FALSE
run.ctmle2 <- TRUE#FALSE#FALSE
compute.true.eic <- FALSE
misspecify.Q <- TRUE
only.A0 <- FALSE
M <- 50

#-------------------------------------------------------------------------------------------#
## get output from simulations
#-------------------------------------------------------------------------------------------#

out <- readRDS(file=paste0("./simulation/output/",
                           "outlist-est",
                           ifelse(run.ltmle, "-ltmle", ""),
                           ifelse(run.ctmle, "-ctmle", ""),
                           ifelse(run.ctmle2, "-ctmle2", ""), 
                           "-2020", "-K", K, ifelse(only.A0, "-A0", ""), ifelse(misspecify.Q, "-Q", ""), 
                           "-M", M, ".rds"))

#-------------------------------------------------------------------------------------------#
## get true values
#-------------------------------------------------------------------------------------------#

psi0.A0 <- readRDS(file=paste0("./simulation/output/",
                               "outlist-est-true-0-2020",
                               "-K", K, ifelse(misspecify.Q, "-Q", ""),
                               "-M", M, ".rds"))

psi0.A1 <- readRDS(file=paste0("./simulation/output/",
                               "outlist-est-true-1-2020",
                               "-K", K, ifelse(misspecify.Q, "-Q", ""),
                               "-M", M, ".rds"))

true.eic.A0 <- readRDS(file=paste0("./simulation/output/",
                                   "outlist-est-true-sd-0-2020",
                                   "-K", K, ifelse(misspecify.Q, "-Q", ""),
                                   "-M", M, ".rds"))

true.eic.A1 <- readRDS(file=paste0("./simulation/output/",
                                   "outlist-est-true-sd-1-2020",
                                   "-K", K, ifelse(misspecify.Q, "-Q", ""),
                                   "-M", M, ".rds"))

#-------------------------------------------------------------------------------------------#
## extract results of interest
#-------------------------------------------------------------------------------------------#

psiA0 <- unlist(lapply(out, function(xout, A=0, which=1) {
    yout <- xout[[A+1]]
    yout[[length(yout)]][which]
}))

psiA0.init <- unlist(lapply(out, function(xout, A=0, which=1) {
    yout <- xout[[A+1]]
    yout[[1]][which]
}))

psiA1 <- unlist(lapply(out, function(xout, A=1, which=1) {
    yout <- xout[[A+1]]
    yout[[length(yout)]][which]
}))


psiA1.init <- unlist(lapply(out, function(xout, A=1, which=1) {
    yout <- xout[[A+1]]
    yout[[1]][which]
}))

sdA0 <- unlist(lapply(out, function(xout, A=0, which=3) {
    yout <- xout[[A+1]]
    yout[[length(yout)]][which]
}))

sdA1 <- unlist(lapply(out, function(xout, A=1, which=3) {
    yout <- xout[[A+1]]
    yout[[length(yout)]][which]
}))

#-------------------------------------------------------------------------------------------#
## extract results of interest
#-------------------------------------------------------------------------------------------#

message("----------------------------")
message("look at estimates from tmle:")
message("------------------")
message(paste0("init (A=0): ", round(mean(psiA0.init), 4)))
message(paste0("tmle (A=0): ", round(mean(psiA0), 4)))
message(paste0("true (A=0): ", round(psi0.A0, 4)))
message("------------------")
message(paste0("init (A=1): ", round(mean(psiA1.init), 4)))
message(paste0("tmle (A=1): ", round(mean(psiA1), 4)))
message(paste0("true (A=1): ", round(psi0.A1, 4)))

message("-----------------------------------")
message("look at coverage of tmle estimator:")
message(paste0("coverage (A=0): ", round(cov.fun(psiA0, sdA0, psi0.A0), 4)))
message(paste0("coverage (A=1): ", round(cov.fun(psiA1, sdA1, psi0.A1), 4)))

message("----------------------------")
message("look at se estimates (A=0):")
message(paste0("mean sigma : ", round(mean(sdA0), 4)))
message(paste0("mse        : ", round(mse(psiA0), 4)))

message("----------------------------")
message("look at se estimates (A=0):")
message(paste0("mean sigma : ", round(mean(sdA1), 4)))
message(paste0("mse        : ", round(mse(psiA1), 4)))

message("-------------------------------")
message("look at efficiency (A=0):")
message(paste0("max-like           : ", round(mse(psiA0)/mean(sdA0), 4)))
message(paste0("max-like, true eic : ", round(mse(psiA0)/true.eic.A0, 4)))

message("-------------------------------")
message("look at efficiency (A=1):")
message(paste0("max-like           : ", round(mse(psiA1)/mean(sdA1), 4)))
message(paste0("max-like, true eic : ", round(mse(psiA1)/true.eic.A1, 4)))
