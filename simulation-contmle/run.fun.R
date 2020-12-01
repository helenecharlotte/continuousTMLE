run.fun <- function(competing.risk=FALSE,
                    treat.effect="both",
                    tau=1.2, 
                    M=1000,
                    n=1000,
                    verbose=FALSE,
                    one.step=FALSE,
                    get.truth=FALSE,
                    setting=1,
                    misspecify.outcome=FALSE,
                    censoring.informative=FALSE,
                    fit.outcome="cox",
                    fit.cens="cox",
                    fit.cr="cox",
                    save.output=TRUE,
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
                    outcome.model <- Surv(time, delta==1)~A+L1.squared
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

        print(psi0.A1 <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                                  categorical=FALSE, competing.risk=competing.risk, cr.both=TRUE,
                                  intervention.A=1, tau=tau, t0=t0,
                                  verbose=FALSE,
                                  reversed.setting=reversed.setting,
                                  square.effect2=square.effect2,
                                  square.effect1=square.effect1,
                                  interaction.Atime=interaction.Atime))
        print(psi0.A0 <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                                  categorical=FALSE, competing.risk=competing.risk, cr.both=TRUE,
                                  intervention.A=0, tau=tau, t0=t0,
                                  verbose=FALSE,
                                  reversed.setting=reversed.setting,
                                  square.effect2=square.effect2,
                                  square.effect1=square.effect1,
                                  interaction.Atime=interaction.Atime))

        psi0 <- psi0.A1 - psi0.A0

        saveRDS(list(psi0=psi0, psi0.A1=psi0.A1, psi0.A0=psi0.A0),
                file=paste0("./simulation-contmle/output/",
                            "save-psi0",
                            ifelse(competing.risk, "-competingrisk", ""), 
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
                                      seed=m+100,
                                      competing.risk=competing.risk,
                                      categorical=FALSE, randomize.A=randomize.A,
                                      censoring.informative=censoring.informative,
                                      censoring=!no.censoring,
                                      square.effect2=square.effect2,
                                      square.effect1=square.effect1,
                                      reversed.setting=reversed.setting,
                                      interaction.Atime=interaction.Atime)
                   
                       #print(paste0("m=", m))

                       #dt[, delta:=(delta==1)*2 + (delta==2)*1]

                       #print(dt[, table(time<=tau, delta, A)])

                       out <- list("est"=contmle(dt, sl.method=3, verbose=verbose,
                                                 change.point=change.point,
                                                 outcome.model=outcome.model,
                                                 cens.model=cens.model,
                                                 cr.model=Surv(time, delta==2)~A+L1+L2+L3,
                                                 cr=competing.risk,
                                                 deps.size=0.1,
                                                 no.small.steps=500,
                                                 separate.cr=FALSE,
                                                 weighted.norm=FALSE,
                                                 one.step=one.step, #deps.size=0.0001,
                                                 fit.outcome=fit.outcome,
                                                 fit.cens=fit.cens,
                                                 fit.cr=fit.cr,
                                                 treat.effect=treat.effect,
                                                 tau=tau,#seq(0.1, 0.6, length=10)[1],#c(0.3,0.5), #c(0.3,0.5,0.8),#tau,
                                                 cut.time=10, V=5,#10,
                                                 cut.time.A=10,
                                                 cut.covars=8, cut.L1.A=8,#3,
                                                 cut.L.interaction=3,
                                                 lambda.cvs=seq(0.0000001, 0.01, length=50),
                                                 #lambda.cvs.cens=seq(0.001, 0.07, length=50),
                                                 output.km=TRUE, output.hr=FALSE,
                                                 sl.models=list(mod1d=c(Surv(time, delta==1)~A+L1+L2+L3, t0=0.9),
                                                                mod2d=c(Surv(time, delta==1)~A+L2.squared+L1*L2+L3, t0=NULL),
                                                                mod3=c(Surv(time, delta==1)~A+L1.squared+L1*L2+L3, t0=0.7),
                                                                mod4=c(Surv(time, delta==1)~A+L2.squared, t0=0.7),
                                                                mod5a=c(Surv(time, delta==1)~A+L1.squared, t0=0.3),
                                                                mod6a=c(Surv(time, delta==1)~A+L1.squared, t0=0.9),
                                                                mod5b=c(Surv(time, delta==1)~A+L1.squared+L2+L3, t0=0.7),
                                                                mod6b=c(Surv(time, delta==1)~A+L2.squared, t0=NULL),
                                                                mod6c=c(Surv(time, delta==1)~A+L1.squared, t0=NULL),
                                                                mod7=c(Surv(time, delta==1)~A+L1+L2+L3, t0=NULL),
                                                                mod8c=c(Surv(time, delta==1)~A*L1+L2+L3, t0=NULL),
                                                                mod9c=c(Surv(time, delta==1)~A*L1.squared+L2+L3, t0=NULL),
                                                                mod10=c(Surv(time, delta==1)~A+L1.squared, t0=0.7)
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
                            "outlist-continuous-tmle",
                            ifelse(competing.risk, "-competingrisk", ""), 
                            paste0("-tau", tau, collapse=""), 
                            paste0("-effect-", ifelse(treat.effect=="both", "ate", paste0("A", treat.effect))),
                            ifelse(n!=1000, paste0("-n", n), ""),
                            ifelse(one.step, "-onestep", ""),
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
