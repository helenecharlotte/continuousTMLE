sim.data <- function(n, loop.max=20, endoffollowup=30,
                     betaA=0.1, betaL=0.3,
                     nu=0.5, eta=4/sqrt(2)*(1/8),
                     firstevent=TRUE,
                     censoring=TRUE,
                     seed=sample(4034244, 1),
                     interaction.AL=FALSE,
                     interaction.Atime=FALSE, t0=0.3,
                     browse=FALSE, verbose=FALSE,
                     randomize.A=FALSE,
                     censoring.informative=TRUE, censoring.high=FALSE, 
                     categorical=TRUE, intervention.A=NULL, tau=2
                     ) {
    
    set.seed(seed)

    if (censoring.high) censoring.alpha <- -1.3 else censoring.alpha <- -2.1

    if (firstevent) loop.max <- 1

    if (length(intervention.A)>0) censoring <- FALSE
    
    if (categorical) {
        Lcont <- runif(n, -3, 3)#rbinom(n, 1, 0.5)
        Lgrid <- seq(min(Lcont), max(Lcont), length=5)    
        L1 <- findInterval(Lcont, Lgrid)/10
    } else {
        #L <- runif(n, -3, 3)/3 #rbinom(n, 1, 0.5)
        L1 <- runif(n, 0, 1)
        L2 <- runif(n, 0, 1)
        L3 <- runif(n, 0, 1) 
    }

    if (interaction.AL | interaction.Atime) {
        L1 <- runif(n, 0, 1) 
        L2 <- runif(n, 0, 1)#rnorm(n, mean=1, 1)
        L3 <- runif(n, 0, 1) 
        if (length(intervention.A)>0) A <- intervention.A else if (randomize.A) A <- rbinom(n, 1, plogis(qlogis(0.5))) else if (categorical)
                                                                                                                           A <- rbinom(n, 1, plogis(0.4+0.3*L1)) else A <- rbinom(n, 1, plogis(0.4+0.3*L1-0.3*L2))
    } else {
        if (length(intervention.A)>0) A <- intervention.A else if (randomize.A) A <- rbinom(n, 1, plogis(qlogis(0.5))) else A <- rbinom(n, 1, plogis(0.4-0.1*L1))
    }


    #-- true density's dependence on covariates/treatment:
    if (interaction.AL & !interaction.Atime) {
        phiT <- function(t, A, L1, L2, L3, betaA, betaL) {
            return(exp(A*betaA+L1^2*betaL-0.7*L1^2-0.25*A*L1^2-0.25*(abs(L2*L1))))
        }
    } else if (interaction.Atime) {
        if (interaction.AL) {
            print("Hallo")
            phiT <- function(t, A, L1, L2, L3, betaA, betaL) {
                return(exp(#-0.45+#0.55*A*(t<=tau/3)-0.65*A*(t>=tau/3)+
                (t<=t0)*betaA*A+
                (t>t0)*(-0.45)*betaA*A*(2.5*L3)-
                L1*betaL-0.6*L2+0.8*L3+#-0.3*L3*L1
                + 0.1))
            }
        } else {
            print("Hallo2")
            phiT <- function(t, A, L1, L2, L3, betaA, betaL) {
                return(exp(#-0.45+#0.55*A*(t<=tau/3)-0.65*A*(t>=tau/3)+
                (t<=t0)*betaA*A+
                (t>t0)*(-0.45)*betaA*A-
                L1*betaL-0.6*L2+0.8*L3
                + 0.1))
            }
            phiT <- function(t, A, L1, L2, L3, betaA, betaL) {
                return(exp(#-0.45+#0.55*A*(t<=tau/3)-0.65*A*(t>=tau/3)+
                (t<=t0)*betaA*A+
                (t>t0)*(-0.45)*betaA*A-
                L1*betaL-1.2*L2+0.8*L3+#-0.3*L3*L1+#0.8*L3
                + 0.1))
            }
        }
    } else {
        phiT <- function(t, A, L1, L2, L3, betaA, betaL) {
            return(exp(A*betaA + L1*betaL))
        }
        #print(phiT)
    }    
    
    lambdaT <- function(t, A, L1, L2, L3, betaA, betaL, eta, nu) {
        return(phiT(t, A, L1, L2, L3, betaA, betaL)*eta*nu*t^{nu-1})
    }

    #-- censoring density's dependence on covariates/treatment:
    if (censoring.informative) {
        print("Yes")
        phiC <- function(t, A, L1, L2, L3) {
            return(exp(-0.1+censoring.alpha))#-A*0.1 - L1*0.4 - L2*0.5 + ifelse(censoring.high, 0.4, 0.7)*censoring.alpha))#7 + L3*0.35 - 2.1))
        }
    } else {
        phiC <- function(t, A, L1, L2, L3) {
            return(exp(-L1*0.2 + censoring.alpha))#7 + L3*0.35 - 2.1))
        }
        print("Here")
        phiC <- function(t, A, L1, L2, L3) {
            return(exp(-L3*1.2+0.8*L1 + 0.9 + censoring.alpha))#7 + L3*0.35 - 2.1))
        }
        phiC <- function(t, A, L1, L2, L3) {
            return(exp(1.2*L1+0.9 + censoring.alpha))#7 + L3*0.35 - 2.1))
        }
        phiC <- function(t, A, L1, L2, L3) {
            return(exp(1.5*L1-0.65*A-0.65*L3*L1+0.15 + censoring.alpha))#7 + L3*0.35 - 2.1))
        }
        phiC <- function(t, A, L1, L2, L3) {
            return(exp(0.5*(2*A-1)*L1+L1+0.15 + censoring.alpha))#7 + L3*0.35 - 2.1))
        }
        phiC <- function(t, A, L1, L2, L3) {
            return(exp(1.2*L1-0.35 + censoring.alpha))#7 + L3*0.35 - 2.1))
        }
        phiC <- function(t, A, L1, L2, L3) {
            return(exp(-L3*1.2+0.8*L1 + 0.9 + censoring.alpha))#7 + L3*0.35 - 2.1))
        }
        phiC <- function(t, A, L1, L2, L3) {
            return(exp(-L3*0.8+1.2*L1*A + 1.1 + censoring.alpha))#7 + L3*0.35 - 2.1))
        }
        #phiC <- function(t, A, L1, L2, L3) {
        #    return(exp(-1.2*L1+1.3+ censoring.alpha))#7 + L3*0.35 - 2.1))
        #}
    }
    
    lambdaC <- function(t, A, L1, L2, L3, eta, nu) {
        return(phiC(t, A, L1, L2, L3)*eta*nu*t^{nu-1})
    }

    if (censoring) {
        phi <- function(t, A, L1, L2, L3, betaA, betaL) phiT(t, A, L1, L2, L3, betaA, betaL) +
                                                            phiC(t, A, L1, L2, L3)
    } else {
        phi <- function(t, A, L1, L2, L3, betaA, betaL) phiT(t, A, L1, L2, L3, betaA, betaL)
    }

    #-- we simulate by using inverse cumulative hazard:
    if (interaction.Atime) {
        Lambda.inv <- function(u, t, A, L1, L2, L3, nu, eta) {
            return( rowSums(cbind((u <= (eta*phi(t, A, L1, L2, L3, betaA, betaL))*t0^{nu}) *
                                  (( (u + eta*phi(t, A, L1, L2, L3, betaA, betaL)*t^{nu}) /
                                     (eta*phi(t, A, L1, L2, L3, betaA, betaL)) )^{1/nu} - t),
            (u > (eta*phi(t, A, L1, L2, L3, betaA, betaL))*t0^{nu}) *
            (( (u - (eta*phi(t, A, L1, L2, L3, betaA, betaL))*t0^{nu} +
                eta*phi(t, A, L1, L2, L3,-0.45*betaA, betaL)*t0^{nu}) /
               (eta*phi(t, A, L1, L2, L3,-0.45*betaA, betaL)) )^{1/nu} - t)), na.rm=TRUE) )
        }
        Lambda.inv1 <- function(u, t, A, L1, L2, L3, nu, eta) {
            return(( (u + eta*phi(t, A, L1, L2, L3, betaA, betaL)*t^{nu}) /
                     (eta*phi(t, A, L1, L2, L3, betaA, betaL)) )^{1/nu} - t)
        }
    } else {
        Lambda.inv <- function(u, t, A, L1, L2, L3, nu, eta) {
            return(( (u + eta*phi(t, A, L1, L2, L3, betaA, betaL)*t^{nu}) /
                     (eta*phi(t, A, L1, L2, L3, betaA, betaL)) )^{1/nu} - t)
        }
    }


    #-- intialize monitoring times: 
    Tlist <- list(cbind(time=rep(0, n), delta=rep(0, n), id=1:n))

    #-- function to loop over for time-points: 
    loop.fun <- function(k, Tprev) {

        #-- simulate event time: 
        U <- -log(runif(n))
        Tout <- Lambda.inv(U, Tprev, A, L1, L2, L3, nu, eta) + Tprev
   
        #-- which event:
        if (censoring) {
            denom <- (lambdaT(Tout, A, L1, L2, L3, betaA, betaL, eta, nu) +
                      lambdaC(Tout, A, L1, L2, L3, eta, nu))
            probT <- lambdaT(Tout, A, L1, L2, L3, betaA, betaL, eta, nu) / (denom)
            probC <- lambdaC(Tout, A, L1, L2, L3, eta, nu) / (denom)
            which <- apply(cbind(probC, probT), 1, function(p) sample(0:1, size=1, prob=p))

            which <- (Tout<=endoffollowup)*which
        } else {
            which <- 1
        }
        
        #-- return event time
        return(cbind(time=Tout, delta=which, id=1:n))
    }

    #-- run simulations: 
    for (k in 1:loop.max) {
        Tlist[[k+1]] <- loop.fun(k, Tlist[[k]][,1])
    }
    
    #-- collect data:
    dt <- data.table(do.call("rbind", Tlist))

    #-- order & throw away obs after end of followup: 
    setorder(dt, id, time, delta)
    dt <- dt[time<=endoffollowup]

    #-- only first event:
    if (firstevent) {
        dt[time>0, idN:=1:.N, by="id"]
        dt <- dt[idN==1][, -"idN", with=FALSE]
        
    }
    
    if (browse) browser()
    
    #-- merge with baseline information:
    baseline <- data.table(A=A, L1=L1, L2=L2, L3=L3, id=as.numeric(1:n))
    setkey(baseline, id); setkey(dt, id)
    dt <- merge(dt, baseline, by="id")

    if (verbose) dt[, hist(time)]

    if (length(intervention.A)>0) return(mean(dt[, time<=tau])) else return(dt)
}
