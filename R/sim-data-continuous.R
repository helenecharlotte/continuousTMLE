sim.data <- function(n, loop.max=20, tau=3000,
                     betaA=0.1, betaL=0.3,
                     nu=0.5, eta=4/sqrt(2)*(1/8),
                     firstevent=TRUE,
                     censoring=TRUE, 
                     categorical=TRUE
                     ) {

    if (categorical) {
        Lcont <- runif(n, -3, 3)#rbinom(n, 1, 0.5)
        Lgrid <- seq(min(Lcont), max(Lcont), length=5)    
        L <- findInterval(Lcont, Lgrid)/10
    } else {
        L <- runif(n, -3, 3)/3 #rbinom(n, 1, 0.5)
    }

    A <- rbinom(n, 1, 0.4+0.1*L)
    
    #-- true density's dependence on covariates/treatment: 
    phiT <- function(t, A, L, betaA, betaL) {
        return(exp(A*betaA + L*betaL))
    }

    lambdaT <- function(t, A, L, betaA, betaL, eta, nu) {
        return(phiT(t, A, L, betaA, betaL)*eta*nu*t^{nu-1})
    }

    #-- censoring density's dependence on covariates/treatment: 
    phiC <- function(t, A, L) {
        return(exp(A*0.01 + L*0.2 - 1.25))
    }

    lambdaC <- function(t, A, L, eta, nu) {
        return(phiC(t, A, L)*eta*nu*t^{nu-1})
    }

    if (censoring) {
        phi <- function(t, A, L, betaA, betaL) phiT(t, A, L, betaA, betaL) + phiC(t, A, L)
    } else {
        phi <- function(t, A, L, betaA, betaL) phiT(t, A, L, betaA, betaL)
    }

    #-- we simulate by using inverse cumulative hazard:
    Lambda.inv <- function(u, t, A, L, nu, eta) {
        return(( (u + eta*phi(t, A, L, betaA, betaL)*t^{nu}) /
                 (eta*phi(t, A, L, betaA, betaL)) )^{1/nu} - t)
    }

    #-- intialize monitoring times: 
    Tlist <- list(cbind(time=rep(0, n), delta=rep(0, n), id=1:n))
    
    #-- function to loop over for time-points: 
    loop.fun <- function(k, Tprev) {

        #-- simulate event time: 
        U <- -log(runif(n))
        Tout <- Lambda.inv(U, Tprev, A, L, nu, eta) + Tprev

        #-- which event:
        if (censoring) {
            denom <- (lambdaT(Tout, A, L, betaA, betaL, eta, nu) +
                      lambdaC(Tout, A, L, eta, nu))
            probT <- lambdaT(Tout, A, L, betaA, betaL, eta, nu) / (denom)
            probC <- lambdaC(Tout, A, L, eta, nu) / (denom)
            which <- apply(cbind(probC, probT), 1, function(p) sample(0:1, size=1, prob=p))

            which <- (Tout<=tau)*which
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

    #-- order & throw away obs after tau: 
    setorder(dt, id, time, delta)
    dt <- dt[time<=tau]

    #-- only first event:
    if (firstevent) {
        dt[time>0, idN:=1:.N, by="id"]
        dt <- dt[idN==1][, -"idN", with=FALSE]
    }
    
    #-- merge with baseline information:
    baseline <- data.table(A=A, L=L, id=as.numeric(1:n))
    setkey(baseline, id); setkey(dt, id)
    dt <- merge(dt, baseline, by="id")

    return(dt)
}
