run.fun <- function(competing.risk=FALSE,
                    treat.effect="ate",
                    tau=1.2, 
                    M=1000,
                    n=1000,
                    verbose=FALSE,
                    one.step=FALSE,
                    iterative=FALSE,
                    get.truth=FALSE,
                    target=1, 
                    setting=1,
                    misspecify.outcome=FALSE,
                    censoring.informative=FALSE,
                    fit.outcome="cox",
                    fit.cens="cox",
                    fit.cr="cox",
                    weighted.norm=FALSE,
                    save.output=TRUE,
                    cr3=FALSE,
                    no_cores=1) {

    if (fit.cens=="hal" | fit.outcome=="hal") {
        no_cores <- 5
    }

    betaA <- -0.15
    betaL <- 1.1
    nu    <- 1.7
    eta   <- 0.7
    t0    <- 0.9
    
    if (setting==1) {
        square.effect1        <- TRUE#FALSE#FALSE#TRUE
        square.effect2        <- FALSE#TRUE#TRUE#FALSE
        interaction.Atime     <- TRUE#FALSE#FALSE#TRUE
        reversed.setting      <- TRUE#FALSE#FALSE#TRUE
    } else {
        square.effect1        <- FALSE#FALSE#TRUE
        square.effect2        <- TRUE#TRUE#FALSE
        interaction.Atime     <- FALSE#FALSE#TRUE
        reversed.setting      <- FALSE#FALSE#TRUE
    }
 
    randomize.A           <- TRUE
    no.censoring          <- FALSE

    misspecify.cens    <- FALSE

    #-------------------------------------------------------------------------------------------#
    ## a bit of reprocessing
    #-------------------------------------------------------------------------------------------#

    if (interaction.Atime) betaA <- -0.7
    change.point <- NULL
    
    if (reversed.setting) {
        betaA <- 0.5
        t0 <- 0.7
    }

    if (square.effect2) {
        if (misspecify.outcome) {
            outcome.model <- Surv(time, delta==1)~A+L1
        } else {
            outcome.model <- Surv(time, delta==1)~A+L1.squared
        }
    }

    if (square.effect1 & !interaction.Atime) {
        if (misspecify.outcome) {
            outcome.model <- Surv(time, delta==1)~A+L1
        } else {
            outcome.model <- Surv(time, delta==1)~A+L1.squared
        }
    }

    if (interaction.Atime) {
        if (misspecify.outcome) {
            change.point <- NULL
            if ((square.effect2 | square.effect1) & reversed.setting) {
                outcome.model <- Surv(time, delta==1)~A+L1
            } else {
                outcome.model <- Surv(time, delta==1)~A+L1+L2+L3
            }
        } else {
            change.point <- t0
            if (square.effect2 | square.effect1) {
                if (reversed.setting) {
                    outcome.model <- Surv(time, delta==1)~A+L1.squared+L2+L3
                } else {
                    outcome.model <- Surv(time, delta==1)~A+L1+L2.squared+L3
                }
            } else {
                outcome.model <- Surv(time, delta==1)~A+L1+L2+L3
            }
        }
    }

    if (misspecify.cens) {
        cens.model <- Surv(time, delta==0)~A+L1+L2+L3
    } else {
        if (censoring.informative) {
            if (square.effect1) {
                cens.model <- Surv(time, delta==0)~L1+L2+L3+A*L1.squared
            } else {
                cens.model <- Surv(time, delta==0)~L1+L2+L3+A*L1
            }
        } else {
            cens.model <- Surv(time, delta==0)~A+L1+L2+L3
        }
    }

    
    #-------------------------------------------------------------------------------------------#
    ## get the true value of target parameter(s)
    #-------------------------------------------------------------------------------------------#

    if (get.truth) {

        source("./R/sim.data.continuous.R")

        par(mfrow=c(1,2))

        psi0.A1 <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                            categorical=FALSE, competing.risk=competing.risk, cr.both=TRUE,
                            intervention.A=1, tau=tau, t0=t0,
                            verbose=FALSE, no.cr=ifelse(cr3, 3, 2),
                            reversed.setting=reversed.setting,
                            square.effect2=square.effect2,
                            square.effect1=square.effect1,
                            interaction.Atime=interaction.Atime)
        psi0.A0 <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                            categorical=FALSE, competing.risk=competing.risk, cr.both=TRUE,
                            intervention.A=0, tau=tau, t0=t0,
                            verbose=FALSE, no.cr=ifelse(cr3, 3, 2),
                            reversed.setting=reversed.setting,
                            square.effect2=square.effect2,
                            square.effect1=square.effect1,
                            interaction.Atime=interaction.Atime)

        psi0 <- psi0.A1[,-1, drop=FALSE] - psi0.A0[,-1, drop=FALSE]

        print(psi0.save <-
                  cbind(tau=psi0.A1[,1], psi0=psi0, psi0.A1=psi0.A1[,-1], psi0.A0=psi0.A0[,-1]))

        saveRDS(psi0.save,
                file=paste0("./simulation-contmle/output/",
                            "save-psi0",
                            ifelse(competing.risk, "-competingrisk", ""),
                            ifelse(competing.risk & cr3, "-cr3", ""), 
                            paste0("-tau", tau, collapse=""), 
                            ifelse(square.effect2, "-squareL2unif", ""),
                            ifelse(square.effect1, "-squareL1unif", ""),
                            ifelse(interaction.Atime, "-interactionAtime", ""),
                            ifelse(reversed.setting, "-reversedsetting", ""),
                            ".rds"))

        return(psi0)

    }

    #-------------------------------------------------------------------------------------------#
    ## repeat simulations (parallelize)
    #-------------------------------------------------------------------------------------------#

    registerDoParallel(no_cores)

    out <- foreach(m=1:M, .errorhandling="pass"
                   ) %dopar% {

                       dt <- sim.data(n, betaA=betaA, betaL=betaL, nu=nu, eta=eta, t0=t0,
                                      seed=m+100, no.cr=ifelse(cr3, 3, 2),
                                      competing.risk=competing.risk,
                                      categorical=FALSE, randomize.A=randomize.A,
                                      censoring.informative=censoring.informative,
                                      censoring=!no.censoring,
                                      square.effect2=square.effect2,
                                      square.effect1=square.effect1,
                                      reversed.setting=reversed.setting,
                                      interaction.Atime=interaction.Atime)

                       if (FALSE) {
                           n2 <- nrow(dt[delta==2])
                           set.seed(1)
                           id2 <- sample(dt[delta==2, id], n2/2)
                           dt[id%in%id2, delta:=3]
                           dt[, table(delta)]
                       }

                       #dt[, delta:=1*(delta==3)+3*(delta==1)+2*(delta==2)]
                       #target <- 1

                       out <- list("est"=contmle(dt, sl.method=3,
                                                 deps.size=0.1, V=10,
                                                 no.small.steps=500,
                                                 weighted.norm=weighted.norm,
                                                 one.step=one.step, 
                                                 target=target, iterative=iterative,
                                                 treat.effect=treat.effect,
                                                 pi.star.fun=function(L) 0, 
                                                 tau=tau, 
                                                 output.km=TRUE,
                                                 verbose=verbose,
                                                 estimation=list("outcome"=list(fit=fit.outcome,
                                                                                model=outcome.model,
                                                                                changepoint=change.point),
                                                                 "cens"=list(fit=fit.cens,
                                                                             model=cens.model,
                                                                             changepoint=NULL),
                                                                 "cr2"=list(fit=fit.cr,
                                                                            model=Surv(time, delta==2)~A+L1+L2+L3,
                                                                            changepoint=NULL),
                                                                 "cr3"=list(fit=fit.cr,
                                                                            model=Surv(time, delta==3)~A+L1+L2+L3,
                                                                            changepoint=NULL)
                                                                 ),
                                                 sl.models=list(
                                                     mod1=c(Surv(time, delta==1)~A+L1+L2+L3, changepoint=c(0.3, 0.7)),
                                                     mod2=c(Surv(time, delta==1)~A+L2.squared+L1*L2+L3, changepoint=NULL),
                                                     mod3=c(Surv(time, delta==1)~A+L1.squared+L1*L2+L3, changepoint=c(0.3, 0.7)),
                                                     mod4=c(Surv(time, delta==1)~A+L2.squared, changepoint=c(0.3, 0.7)),
                                                     mod5=c(Surv(time, delta==1)~A+L1.squared, changepoint=c(0.3, 0.7)),
                                                     mod6=c(Surv(time, delta==1)~A+L1.squared+L2+L3, changepoint=c(0.3, 0.7)),
                                                     mod7=c(Surv(time, delta==1)~A+L2.squared, changepoint=NULL),
                                                     mod8=c(Surv(time, delta==1)~A+L1.squared, changepoint=NULL),
                                                     mod9=c(Surv(time, delta==1)~A+L1+L2+L3, changepoint=NULL),
                                                     mod10=c(Surv(time, delta==1)~A*L1+L2+L3, changepoint=NULL),
                                                     mod11=c(Surv(time, delta==1)~A*L1.squared+L2+L3, changepoint=NULL)
                                                 )))
                       print(paste0("m=", m))
                       names(out) <- paste0("m=", m)
                       print(out)
                       
                       return(out)

                   }

    stopImplicitCluster()


    if (save.output) {
        saveRDS(out,
                file=paste0("./simulation-contmle/output/",
                            "outlist-contmle",
                            ifelse(competing.risk, "-competingrisk", ""),
                            ifelse(competing.risk & cr3, "-cr3", ""), 
                            paste0("-tau", tau, collapse=""), 
                            paste0("-effect-", ifelse(treat.effect%in%c("ate","both"), "ate", paste0("A", treat.effect))),
                            ifelse(n!=1000, paste0("-n", n), ""),
                            ifelse(one.step, "-onestep", ""),
                            ifelse(competing.risk & iterative & length(tau)==1, "-iterative", ""),
                            ifelse(is.character(weighted.norm), paste0("-", weighted.norm), ""),
                            ifelse(competing.risk, paste0("-target", target, collapse="-"), ""),
                            ifelse(square.effect2, "-squareL2unif", ""),
                            ifelse(square.effect1, "-squareL1unif", ""),
                            ifelse(interaction.Atime, "-interactionAtime", ""),
                            ifelse(reversed.setting, "-reversedsetting", ""),
                            ifelse(randomize.A, "-randomizeA", ""),
                            ifelse(censoring.informative, "-informativeCensoring", ""),
                            ifelse(no.censoring, "-noCensoring", ""),
                            ifelse(misspecify.outcome, "-misspecify-outcome", ""),
                            ifelse(misspecify.cens, "-misspecify-cens", ""),
                            paste0("-outcome-fit", fit.outcome),
                            paste0("-censoring-fit", fit.cens),
                            ifelse(competing.risk, paste0("-cr-fit", fit.cr), ""), 
                            "-M", M, ".rds"))
    } else {
        return(out)
    }


} 
