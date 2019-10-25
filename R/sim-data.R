sim.data <- function(n=50, tau=365.24, nL=3, nL0=1,  
                     loop.max=20,
                     eta = 4/sqrt(2)*(1/8), # scale parameter
                     nu = 0.5, # shape parameter
                     alpha.Y=1.55, eta.Y=0.07,
                     alpha.N=1.2, eta.N=0.05,
                     alpha.C=1.1, eta.C=0.07,
                     betaL.Y=c(0.3, -0.1, 0),
                     betaA.Y=-0.2, betaA.C=0,
                     betaL.C=c(0, 0, 0),
                     betaL0.Y=0.8,
                     betaAL0.Y=0,
                     randomize.treatment=FALSE,
                     set.treatment=NULL,
                     no.censoring=FALSE,
                     browse=FALSE, seed=sample(4034244, 1),
                     L.A=TRUE,
                     compute.target=NULL,
                     dt.out=TRUE,
                     compute.true.ic=FALSE) {

    set.seed(seed)
   
    rexpit <- function(x) rbinom(n=n, size=1, prob=plogis(x))

    if (browse) browser()
    
    #-- baseline covariate and baseline treatment:
    if (nL0>1) {
        L0 <- sapply(1:nL0, function(l) runif(n))
        colnames(L0) <- paste0("L0.", 1:nL0)
    } else {
        L0 <- runif(n, 0, 1)
    }

    
    if (length(compute.target)>0) {
        A  <- rep(compute.target, n)
    } else if (randomize.treatment) {
        A <- rexpit(0.5)
    } else if (nL0>1) {
        A  <- rexpit(0.5+0.15*L0[, 1]-0.25*L0[, 2])
    } else {
        A  <- rexpit(0.5+0.01*L0)
    }

    if (length(set.treatment)>0) {
        A <- rep(set.treatment, n)
    }    
    
    #-- initial values of time-varying covariates:
    L <- sapply(1:nL, function(l) runif(n))
    colnames(L) <- paste0("L", 1:nL)

    #------- intensities:
   
    #-- phi for outcome: 
    phiY <- function(t, A, L0, L1, L2, L3, k, etaY=eta.Y, alphaY=alpha.Y,
                     betaA=betaA.Y, betaL0=betaL0.Y, betaL=betaL.Y, betaAL0=betaAL0.Y) {
        if (nL0>1){
            add.L0 <- apply(do.call("cbind", lapply(1:length(betaL0), function(j) betaL0[j] * L0[, j])),
                            1, sum)
            return(etaY*alphaY^{(k>2)}*exp(add.L0 + betaAL0*A*L0[,4] + betaL[1]*L1 +
                                           betaL[2]*L2 + betaL[3]*L3 + betaA*A))
        } else {
            return(etaY*alphaY^{(k>2)}*exp(betaL0*L0 + betaL[1]*L1 + betaL[2]*L2 + betaL[3]*L3 + betaA*A))
        }
    }

    #-- phi for time-varying covariates: 
    phiN <- function(t, A, L0, L1, L2, L3, k, etaN=eta.N, alphaN=alpha.N,
                     betaA=0.4, betaL0=0.5, Kmax=4) {
        if (nL0>1){
            return(etaN*alphaN^{min(k, Kmax)}*exp(betaL0*L0[,1] + betaA*A))
        } else {
            return(etaN*alphaN^{min(k, Kmax)}*exp(betaL0*L0 + betaA*A))
        }
    }

    #-- phi for time-varying censoring:
    if (length(compute.target)>0 | no.censoring) {
        phiC <- function(t, A, L0, L1, L2, L3, k, etaC=eta.C, alphaC=alpha.C,
                         betaA=betaA.C, betaL0=0.1, betaL.C) {
            return(0)
        }
    } else {
        phiC <- function(t, A, L0, L1, L2, L3, k, etaC=eta.C, alphaC=alpha.C,
                         betaA=betaA.C, betaL0=0.1, betaL=betaL.C) {
            if (nL0>1){
                return(etaC*alphaC^{(k>2)}*exp(betaL0*L0[, 3] + betaA*A + betaL[1]*L1 +
                                               betaL[2]*L2 + betaL[3]*L3))
            } else{
                return(etaC*alphaC^{(k>2)}*exp(betaL0*L0 + betaA*A + betaL[1]*L1 +
                                               betaL[2]*L2 + betaL[3]*L3))
            }
        }
    }
    
    #-- collected phi: 
    phi <- function(t, A, L0, L1, L2, L3, k) {
        phiY(t, A, L0, L1, L2, L3, k) + phiN(t, A, L0, L1, L2, L3, k) + phiC(t, A, L0, L1, L2, L3, k)
    }

    #-- intensity for outcome: 
    lambdaY <- function(t, A, L0, L1, L2, L3, k, etaY=eta, nuY=nu) {
        phiY(t, A, L0, L1, L2, L3, k)*etaY*nuY*t^{nuY-1}
    }

    #-- intensity for monitoring of time-varying covariates: 
    lambdaN <- function(t, A, L0, L1, L2, L3, k, etaN=eta, nuN=nu) {
        phiN(t, A, L0, L1, L2, L3, k)*etaN*nuN*t^{nuN-1}
    }

    #-- intensity for censoring: 
    lambdaC <- function(t, A, L0, L1, L2, L3, k, etaC=eta, nuC=nu) {
        phiC(t, A, L0, L1, L2, L3, k)*etaC*nuC*t^{nuC-1}
    }

    #-- we simulate by using inverse cumulative hazard:
    Lambda.inv <- function(u, t, k, L, A, L0) {
        L1 <- L[, 1]; L2 <- L[, 2]; L3 <- L[, 3]
        return(( (u + eta*phi(t, A, L0, L1, L2, L3, k)*t^{nu}) /
                 (eta*phi(t, A, L0, L1, L2, L3, k)) )^{1/nu} - t)
    }

    #-- intialize monitoring times: 
    Tlist <- list(cbind(time=rep(0, n), delta=rep(1, n), id=1:n, L=L))

    #-- function to loop over for time-points: 
    loop.fun <- function(k, Tprev, Lprev) {

        #-- simulate event time: 
        U <- -log(runif(n))
        Tout <- Lambda.inv(U, Tprev, k, Lprev, A, L0) + Tprev

        L1 <- Lprev[, 1]; L2 <- Lprev[, 2]; L3 <- Lprev[, 3]
    
        #-- character of event:
        denom <- (lambdaY(Tout, A, L0, L1, L2, L3, k) +
                  lambdaN(Tout, A, L0, L1, L2, L3, k) +
                  lambdaC(Tout, A, L0, L1, L2, L3, k))
        probY <- lambdaY(Tout, A, L0, L1, L2, L3, k) / (denom)
        probN <- lambdaN(Tout, A, L0, L1, L2, L3, k) / (denom)
        probC <- lambdaC(Tout, A, L0, L1, L2, L3, k) / (denom)
        which <- apply(cbind(probY, probN, probC), 1, function(p) sample(0:2, size=1, prob=p))

        which <- (Tout>tau)*3 + (Tout<=tau)*which
        Tout <- (Tout>tau)*tau + (Tout<=tau)*Tout

        if (L.A) {
            Lout <- (which==1)*sapply(1:nL, function(l) rnorm(n, mean=A*(l-1)*0.02 +
                                                                     L0/l*0.05 + Lprev[, l], sd=0.1)) +
                (which!=1)*Lprev
        } else {
            Lout <- (which==1)*sapply(1:nL, function(l) rnorm(n, mean=L0/l*0.05 + Lprev[, l], sd=0.1)) +
                (which!=1)*Lprev
        }
        return(cbind(time=Tout, delta=which, id=1:n, L=Lout))
    }

    #-- run simulations: 
    for (k in 1:loop.max) {
        Tlist[[k+1]] <- loop.fun(k, Tlist[[k]][,1], Tlist[[k]][,4:(4+nL-1)])
    }

    #-- collect data:
    dt <- data.table(do.call("rbind", Tlist))

    #-- order & throw away observations after censoring/event: 
    setorder(dt, id, time)
    dt[, keep:=c(1, cumprod(delta==1)[-.N]), by="id"]
    dt <- dt[keep==1][, -"keep"]

    #-- merge with baseline information:
    setkey(dt, id)
    if (nL0>1) {
        baseline <- data.table(A=A, L0, id=as.numeric(1:n))
    } else {
        baseline <- data.table(A=A, L0=L0, id=as.numeric(1:n))
    }
    setkey(dt, id)
    dt <- merge(dt, baseline, by="id")

    if (compute.true.ic) {

        if (length(compute.target)>0) {
            compute.target <- 0
        }
        dt[, idN:=1:.N, by="id"]
        dt[, N:=.N, by="id"]

        n <- length(dt[, unique(id)])
        unique.times <- sort(unique(dt[, time]))

        dt.all <- do.call("rbind", lapply(dt[, unique(id)], function(i) {
            dt.i <- dt[id==i]
            out.i <- rbind(dt.i, data.table(time=unique.times[!(unique.times %in% dt.i[, time])]),
                           fill=TRUE)
            out.i <- out.i[order(time)]
            out.i[is.na(delta), delta:=4]
            fill.na(out.i)
            return(out.i)
        }))

        dt[, k:=c(0, idN[-.N]), by="id"]

        LambdaY <- function(t, A, L0, L1, L2, L3, k, etaY=eta, nuY=nu) {
            phiY(t, A, L0, L1, L2, L3, k)*etaY*t^{nuY}
        }

        LambdaC <- function(t, A, L0, L1, L2, L3, k, etaC=eta, nuC=nu) {
            phiC(t, A, L0, L1, L2, L3, k)*etaC*t^{nuC}
        }
        
        dt.all[, lambda.Y:=lambdaY(time, A, cbind(L0.1,L0.2,L0.3,L0.4), L1, L2, L3, k, etaY=eta, nuY=nu)]
        dt.all[, Lambda.Y:=LambdaY(time, A, cbind(L0.1,L0.2,L0.3,L0.4), L1, L2, L3, k, etaY=eta, nuY=nu)]
        dt.all[, Lambda.C:=LambdaC(time, A, cbind(L0.1,L0.2,L0.3,L0.4), L1, L2, L3, k, etaC=eta, nuC=nu)]
        dt.all[, pi.A:=plogis(0.5+0.15*L0.1-0.25*L0.2)]

        dt.all[, surv.Y:=exp(-Lambda.Y), by="id"]
        dt.all[, surv.C.tminus:=c(1, exp(-Lambda.C)[-.N]), by="id"]
        dt.all[, surv.Y.t0:=surv.Y[.N], by="id"]

        dt.all[, time.diff:=c(time[-1], time[.N])-time, by="id"]

        dt.all[, keep:=cumprod(delta != 0 & delta!= 2), by="id"]
        dt.all[, keep:=c(1, keep[-.N]), by="id"]

        dt.all[, Ht:=(keep) * (A==compute.target) *
                     1/ ( pi.A^compute.target*(1-pi.A)^compute.target * surv.C.tminus )]
        dt.all[, Ht.lambda:=-surv.Y.t0/surv.Y]
        
        out <- dt.all[, sum( (keep) * Ht * Ht.lambda *
                             ( (delta==0) - lambda.Y * time.diff
                             )), by="id"]

        ic.squared <- out[, 2][[1]]^2
        sqrt(mean(ic.squared)/n)

        dt.all[, test:=cumsum(lambda.Y*time.diff), by="id"]
        dt.all[ , c("lambda.Y", "Lambda.Y", "test")]
    }
    
    #-- output dataset:
    if (dt.out) {
        # if events at same date:
        dt[, time.round:=round(time, 6)]
        dt[, N.tmp:=.N, by=c("id", "time.round")]
        if (dt[, any(N.tmp>1)]) {
            dt[, any.Y:=any(delta==0), by="id"]
            dt[, any.C:=any(delta==2), by="id"]
            if (dt[, any(N.tmp>1 & any.C & delta %in% c(1,4))]) {
                dt <- dt[!(N.tmp>1 & any.C & delta %in% c(1,4))]
            } 
            if (dt[, any(N.tmp>1 & any.Y & delta!=0)]) {
                dt <- dt[!(N.tmp>1 & any.Y & delta!=0)]
            }
            dt <- dt[, -c("any.C", "any.Y"), with=FALSE]
        }
        dt <- dt[, -c("time.round", "N.tmp"), with=FALSE]
        return(dt)
    } else if (compute.true.ic) {
        out <- sqrt(mean(ic.squared)/n) # OBS: if want sd, divide by sqrt(n).
        return(out)
    } else {
        dt[, Y:=1*(delta[.N]==0), by="id"]
        out <- mean(dt[time==0, Y])
        return(out)
    }
}








