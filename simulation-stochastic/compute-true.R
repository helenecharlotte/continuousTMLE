#-------------------------------------------------------------------------------------------#
## true values of parameters
#-------------------------------------------------------------------------------------------#

compute.true <- function(compute.true.psi=FALSE, compute.true.eic=FALSE) {
                         
    if (compute.true.psi) {

        #--- Causal parameter; 
    
        print(psi0.M0 <- sim.data(1e6, seed=10011,
                                  rescue=TRUE,
                                  intervention.A0=function(L0, A0) cbind(logit(1)),
                                  K=K
                                  ))

        print(psi0.M1 <- sim.data(1e6, seed=10011,
                                  rescue=TRUE,
                                  intervention.A0=function(L0, A0) cbind(logit(0)),
                                  K=K
                                  ))

        saveRDS(psi0.M0,
                file=paste0("./simulation-stochastic/output/",
                            "outlist-est-true-rescue-0-2020",
                            "-K", K, 
                            "-M", M, ".rds"))

        saveRDS(psi0.M1,
                file=paste0("./simulation-stochastic/output/",
                            "outlist-est-true-rescue-1-2020",
                            "-K", K,
                            "-M", M, ".rds"))

        print(psi0.M1 - psi0.M0)

        #--- ITT parameter; 

        print(psi0.itt.M0 <- sim.data(1e6, seed=10011,
                                      intervention.A0=function(L0, A0) cbind(logit(1)),
                                      rescue=TRUE,
                                      rescue.itt=TRUE,
                                      K=K
                                      ))

        print(psi0.itt.M1 <- sim.data(1e6, seed=10011,
                                      intervention.A0=function(L0, A0) cbind(logit(0)),
                                      rescue=TRUE,
                                      rescue.itt=TRUE,
                                      K=K
                                      ))

        saveRDS(psi0.itt.M0,
                file=paste0("./simulation-stochastic/output/",
                            "outlist-est-true-rescue-itt-0-2020",
                            "-K", K, 
                            "-M", M, ".rds"))

        saveRDS(psi0.itt.M1,
                file=paste0("./simulation-stochastic/output/",
                            "outlist-est-true-rescue-itt-1-2020",
                            "-K", K,
                            "-M", M, ".rds"))

        print(psi0.itt.M1 - psi0.itt.M0)
    }

    if (compute.true.eic) { #--- FIXME: not adapted yet.

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
                file=paste0("./simulation-stochastic/output/",
                            "outlist-est-true-sd-0-2020",
                            "-K", K, ifelse(only.A0, "-A0", ""),
                            "-M", M, ".rds"))

        saveRDS(true.eic1,
                file=paste0("./simulation-stochastic/output/",
                            "outlist-est-true-sd-1-2020",
                            "-K", K, ifelse(only.A0, "-A0", ""),
                            "-M", M, ".rds"))
    }
} 
