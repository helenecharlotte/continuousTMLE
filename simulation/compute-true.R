#-------------------------------------------------------------------------------------------#
## true values of parameters
#-------------------------------------------------------------------------------------------#

compute.true <- function(compute.true.psi=FALSE, compute.true.eic=FALSE, only.A0=FALSE) {
                         
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
                            "-K", K, ifelse(only.A0, "-A0", ""),
                            "-M", M, ".rds"))

        saveRDS(psi0.test.multi.M1,
                file=paste0("./simulation/output/",
                            "outlist-est-true-1-2020",
                            "-K", K, ifelse(only.A0, "-A0", ""),
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
                                                only.A0=only.A0,
                                                smooth.initial=TRUE,
                                                browse9=FALSE, maxIter=1, 
                                                compute.true.eic=TRUE,
                                                intervention.A0=function(L0, A0) logit(1*(A0==0)),
                                                intervention.A=function(L0, A0, L.prev, A.prev, A) logit(1*(A==0)),
                                                browse0=FALSE, misspecify.Q=misspecify.Q))
    
        true.eic.1 <- suppressMessages(est.fun(copy(dt), censoring=TRUE,
                                               targeting=1,
                                               only.A0=only.A0,
                                               smooth.initial=TRUE,
                                               browse9=FALSE, maxIter=1, 
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
                            "-K", K, ifelse(only.A0, "-A0", ""),
                            "-M", M, ".rds"))

        saveRDS(true.eic1,
                file=paste0("./simulation/output/",
                            "outlist-est-true-sd-1-2020",
                            "-K", K, ifelse(only.A0, "-A0", ""),
                            "-M", M, ".rds"))
    }
} 
