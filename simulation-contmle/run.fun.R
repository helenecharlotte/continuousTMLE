run.fun <- function(competing.risk=FALSE,
                    treat.effect="ate",
                    tau=1.2, 
                    M=1000,
                    n=1000,
                    verbose=FALSE,
                    check.sup=FALSE,
                    use.observed.times=FALSE, length.times=length(tau), check.times.size=NULL,
                    one.step=FALSE,
                    push.criterion=FALSE,
                    simultaneous.ci=FALSE,
                    iterative=FALSE,
                    get.truth=FALSE,
                    target=1, 
                    setting=1,
                    print.forms=FALSE,
                    add.noise=NULL, 
                    misspecify.outcome=FALSE,
                    censoring.informative=FALSE,
                    no.effect.A=FALSE,
                    randomize.A=TRUE,
                    fit.outcome="cox",
                    fit.cens="cox",
                    fit.cr="cox",
                    cox.hal=FALSE,
                    cut.covars=8,
                    hal.screening=FALSE,
                    return.data=NULL,
                    V=10, lambda.cvs=seq(0.0000001, 0.01, length=50),
                    weighted.norm=FALSE,
                    save.output=TRUE,
                    cr3=FALSE,
                    browse=FALSE,
                    pi.star.fun=function(L) 0,
                    no_cores=1, override.nc=FALSE) {

    if (any(fit.cens=="hal") | any(fit.outcome=="hal") | any(fit.cr=="hal") & n>500) {
        if (!override.nc) no_cores <- 5
    } 

    if (!no.effect.A) betaA <- -0.15 else betaA <- 0
    betaL <- 1.1
    nu    <- 1.7
    eta   <- 0.7
    t0    <- 0.9
    
    if (setting==1 | setting==3) {
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
 
    no.censoring          <- FALSE

    misspecify.cens    <- FALSE

    #-------------------------------------------------------------------------------------------#
    ## a bit of reprocessing
    #-------------------------------------------------------------------------------------------#

    if (interaction.Atime & !no.effect.A) betaA <- -0.7
    change.point <- NULL
    
    if (reversed.setting) {
        if (!no.effect.A) betaA <- 0.5
        t0 <- 0.7
    }

    if (setting==3) {
        betaA <- 0.5
        t0 <- 0.8
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

    #--- "back-up models" for hal; 

    if (fit.cens[1]=="hal") {
        cens.model <- Surv(time, delta==0)~A#A+L1+L2+L3
    }

    if (fit.outcome[1]=="hal") {
        outcome.model <- Surv(time, delta==1)~A##A+L1+L2+L3
        change.point <- NULL
    }

    if (fit.cr[1]=="hal") {
        cr.model <- Surv(time, delta==2)~A+L1+L2+L3
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
                            print.forms=print.forms,
                            interaction.Atime=interaction.Atime)
        psi0.A0 <- sim.data(1e6, betaA=betaA, betaL=betaL, nu=nu, eta=eta,
                            categorical=FALSE, competing.risk=competing.risk, cr.both=TRUE,
                            intervention.A=0, tau=tau, t0=t0,
                            verbose=FALSE, no.cr=ifelse(cr3, 3, 2),
                            reversed.setting=reversed.setting,
                            square.effect2=square.effect2,
                            square.effect1=square.effect1,
                            print.forms=print.forms,
                            interaction.Atime=interaction.Atime)

        psi0 <- psi0.A1[,-1, drop=FALSE] - psi0.A0[,-1, drop=FALSE]

        print(psi0.save <-
                  cbind(tau=psi0.A1[,1], psi0=psi0, psi0.A1=psi0.A1[,-1], psi0.A0=psi0.A0[,-1]))

        saveRDS(psi0.save,
                file=paste0("./simulation-contmle/output/",
                            "save-psi0",
                            ifelse(setting==3, "-setting3", ""),
                            ifelse(no.effect.A, "-no-effect-A", ""),
                            ifelse(competing.risk, "-competingrisk", ""),
                            ifelse(competing.risk & cr3, "-cr3", ""), 
                            ifelse(length(tau)>10,
                                   paste0("-tau", length(tau),
                                          paste0("-tau", round(tau,3)[1:3], collapse="")),
                                   paste0("-tau", round(tau,3), collapse="")), 
                            ifelse(square.effect2, "-squareL2unif", ""),
                            ifelse(square.effect1, "-squareL1unif", ""),
                            ifelse(interaction.Atime, "-interactionAtime", ""),
                            ifelse(reversed.setting, "-reversedsetting", ""),
                            ".rds"))

        return(psi0)

    }

    if (length(return.data)>0) {
        dt <- sim.data(n, betaA=betaA, betaL=betaL, nu=nu, eta=eta, t0=t0,
                       seed=return.data+100, no.cr=ifelse(cr3, 3, 2),
                       competing.risk=competing.risk,
                       categorical=FALSE, randomize.A=randomize.A,
                       censoring.informative=censoring.informative,
                       censoring=!no.censoring,
                       square.effect2=square.effect2,
                       square.effect1=square.effect1,
                       print.forms=print.forms,
                       reversed.setting=reversed.setting,
                       interaction.Atime=interaction.Atime)

        if (length(add.noise)>0) {

            if (!is.numeric(add.noise)) add.noise <- 10            

            L.noise <- matrix(sapply(1:add.noise, function(jj) rnorm(nrow(dt), mean=jj/add.noise, sd=0.2)), nrow=nrow(dt))
            dt <- cbind(dt, L.noise)
            
        }

        return(dt)         
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
                                      print.forms=print.forms,
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
                       # dt[, delta:=1*(delta==0)+0*(delta==1)]
                       #target <- 1

                       if (browse) browser()

                       if (length(add.noise)>0) {

                           if (!is.numeric(add.noise)) add.noise <- 10            

                           L.noise <- matrix(sapply(1:add.noise, function(jj) rnorm(nrow(dt), mean=jj/add.noise, sd=0.2)), nrow=nrow(dt))
                           dt <- cbind(dt, L.noise)

                           outcome.model <- as.formula(paste0(paste0(outcome.model[2],
                                                                     outcome.model[1],
                                                                     outcome.model[3], "+",
                                                                     paste0(paste0("V", 1:10), collapse="+"))))
                           
                       }
                       
                       out <- list("est"=contmle(dt, sl.method=3,
                                                 deps.size=0.1, V=V,
                                                 no.small.steps=500,
                                                 weighted.norm=weighted.norm,
                                                 one.step=one.step,
                                                 push.criterion=push.criterion,
                                                 simultaneous.ci=simultaneous.ci,
                                                 target=target, iterative=iterative,
                                                 treat.effect=treat.effect,
                                                 pi.star.fun=pi.star.fun, 
                                                 tau=tau,
                                                 cut.covars=cut.covars,
                                                 hal.screening=hal.screening,
                                                 lambda.cvs=lambda.cvs,
                                                 output.km=TRUE,
                                                 verbose=verbose,
                                                 check.sup=check.sup,
                                                 use.observed.times=use.observed.times,
                                                 length.times=length.times,
                                                 check.times.size=check.times.size,
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
                                                 sl.models=list(mod1=list(Surv(time, delta==1)~A+L1+L2+L3),
                                                                mod2=list(Surv(time, delta==1)~A+L1.squared+L2+L3),
                                                                mod3=list(Surv(time, delta==1)~L2.squared+A+L1.squared+L2+L3),
                                                                mod4=list(Surv(time, delta==1)~A+L1.squared),
                                                                mod5=list(Surv(time, delta==1)~A*L1+L2+L3),
                                                                mod6=list(Surv(time, delta==1)~A*L1.squared+L2+L3))))
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
                            ifelse(setting==3, "-setting3", ""),
                            ifelse(no.effect.A, "-no-effect-A", ""),
                            ifelse(competing.risk, "-competingrisk", ""),
                            ifelse(competing.risk & cr3, "-cr3", ""),
                            ifelse(is.character(weighted.norm), paste0("-", weighted.norm), ""),
                            paste0("-effect-", ifelse(treat.effect%in%c("ate","both"), "ate", paste0("A", treat.effect))),
                            "-M", M, ".rds"))
        saveRDS(out,
                file=paste0("./simulation-contmle/output/",
                            "outlist-contmle",
                            ifelse(setting==3, "-setting3", ""),
                            ifelse(no.effect.A, "-no-effect-A", ""),
                            ifelse(competing.risk, "-competingrisk", ""),
                            ifelse(competing.risk & cr3, "-cr3", ""), 
                            ifelse(use.observed.times, paste0("-use-observed-times-length", length.times),
                            ifelse(length(tau)>10 | (max(tau)>1.2 & length(tau)==10),
                                   paste0("-tau", length(tau),
                                          paste0("-tau", round(tau,3)[1:3], collapse="")),
                                   paste0("-tau", round(tau,3), collapse=""))), 
                            paste0("-effect-", ifelse(treat.effect%in%c("ate","both"), "ate", paste0("A", treat.effect))),
                            ifelse(n!=1000, paste0("-n", n), ""),
                            ifelse(one.step, "-onestep", ""),
                            ifelse(!competing.risk & iterative & length(tau)>1, "-iterative", ""),
                            ifelse(competing.risk & (iterative & length(target)>1) & length(tau)==1, "-iterative", ""),
                            ifelse(is.character(weighted.norm), paste0("-", weighted.norm), ""),
                            ifelse(competing.risk, paste0("-target", target, collapse="-"), ""),
                            ifelse(square.effect2, "-squareL2unif", ""),
                            ifelse(square.effect1, "-squareL1unif", ""),
                            ifelse(interaction.Atime, "-interactionAtime", ""),
                            ifelse(reversed.setting, "-reversedsetting", ""),
                            ifelse(randomize.A, "-randomizeA", ""),
                            ifelse(censoring.informative, "-informativeCensoring", ""),
                            ifelse(no.censoring, "-noCensoring", ""),
                            ifelse(all(fit.outcome=="cox") & misspecify.outcome, "-misspecify-outcome", ""),
                            ifelse((any(fit.outcome=="hal") | any(fit.cens=="hal") | any(fit.cr=="hal")) & cox.hal, "-cox-hal", ""),
                            ifelse(all(fit.cens=="cox") & misspecify.cens, "-misspecify-cens", ""),
                            paste0("-outcome-fit", paste0(fit.outcome, collapse="")),
                            paste0("-censoring-fit", paste0(fit.cens, collapse="")),
                            ifelse(competing.risk, paste0("-cr-fit", paste0(fit.cr, collapse="")), ""), 
                            "-M", M, ".rds"))
    } else {
        return(out)
    }


} 
