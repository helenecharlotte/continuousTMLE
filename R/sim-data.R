sim.data <- function(n, seed=sample(4034244, 1),
                              censoring=FALSE,
                              intervention.A=NULL, 
                              intervention.dN.A=NULL,
                              time.since=3,
                              verbose=FALSE, q0=3,
                              misspecify.Q=FALSE,
                              only.A0=FALSE,
                              form.dN.L=NULL,
                              form.dN.A=NULL,
                              form.A=NULL, 
                              form.Y=NULL,
                              form.C=NULL,
                              form.L=NULL,
                              form.A0=NULL,
                              intervention.A0=NULL,
                              rescue=FALSE, ignore.rescue=FALSE, 
                              obs.mean=FALSE,
                              K=1,
                              browse=FALSE) {


    if (length(form.dN.L)==0) form.dN.L <- function(L0, dN.L.prev, L.prev, A.prev) -0.2-0.05*K-0.025*(K>7)-0.25*dN.L.prev-0.15*L0-0.1*(A.prev==1)+0.3*L.prev
    if (length(form.dN.A)==0) form.dN.A <- function(L0, dN.A.prev, L.prev, A.prev, no.jumps.A, L.star) -0.75-0.05*K-0.42*dN.A.prev+0.15*L0+0.3*(A.prev==2)+0.4*(A.prev==1)-0.25*L.prev
    if (length(form.C)==0) form.C <- function(L0, L.prev, A.prev, A0) -3.95+(K>40)*5-0.4*K^{2/3}-0.24*(K>2 & K<=4)-0.4*(K>4 & K<=9)-(K>9)*0.4*K^{1/5}+0.2*(K>25)*K^{1/4}+0.1*L0+0.2*(A0==1)+0.9*(A0==2)+2.15*L.prev
    if (length(form.Y)==0) form.L <- function(L0, L.prev, A.prev, A0) 0.5-0.4*A0+0.15*L0-0.25*(A.prev==1)+0.4*L.prev
    if (length(form.A0)==0) form.A0 <- function(L0) cbind(-0.1+0.25*L0)
    if (length(form.A)==0) form.A <- function(L0, L.prev, A.prev, A0) cbind(-1+(1-A0)*0.6+(1-A.prev)*0.4+L.prev*0.6-0.15*(K>15)*L.prev)
    if (length(form.Y)==0) {
        if (!only.A0) {
            form.Y <- function(L0, L.prev, A.prev, A0, no.jumps.A, dN.A.prev) -1.1-
                0.33*K/3*(K>2 & K<=4)-0.25*K^{2/3}-0.25*(K>4 & K<=9)-
                (K>25 & K<45)*0.3*K^{1/5}-
                (K>75)*0.31+(K>85)*0.2-
                (K>25 & K<75)*0.5*K^{1/5}+0.6*(K>25)*K^{1/4}-0.25*A.prev+
                0.4*L.prev-0.25*A0+0.35*L.prev*A0+(K>75)*0.1*A0+(K>85)*0.01*A0
        } else {
            form.Y <- function(L0, L.prev, A.prev, A0, no.jumps.A, dN.A.prev) -1.1-
                0.33*K/3*(K>2 & K<=4)-0.25*K^{2/3}-0.25*(K>4 & K<=9)-
                (K>25 & K<45)*0.3*K^{1/5}-
                (K>75)*0.31+(K>85)*0.2-
                (K>25 & K<75)*0.5*K^{1/5}+0.6*(K>25)*K^{1/4}+L0*0.2-0.25*A0+
                (K>75)*0.1*A0+(K>85)*0.01*A0
        }
    }
    
    if (length(seed)>0) {
        set.seed(seed)
    }
    
    rexpit <- function(x) rbinom(n=n, size=1, prob=plogis(x))
    rmulti <- function(x) {
        out <- rep(0, n)
        for (p in 1:ncol(x)) {
            out <- p*rexpit(logit(x[,p]))*(out==0) + out 
        }
        out[out==0] <- ncol(x)+1
        return(out-1)
    }

    print(k.grid <- ceiling(seq(0, K, length=q0+1)[-c(1,q0+1)]))
    
    ## rmulti <- function(x) {
    ##     apply(x, 1, function(p) {
    ##         browser()
    ##         sample(0:2, prob=c(plogis(p[1])plogis(p[2]), plogis(p[1]), plogis(p[2])), size=1)
    ##     })}

    L0 <- #runif(n)
        sample(1:6, n, replace=1000)/6

    ## if (length(true.value)>0) {
    ##     A0 <- rep(true.value, n)
    ## } else
    A0 <- rmulti(plogis(form.A0(L0)))
    
    if (length(intervention.A)>0 & length(intervention.A0)==0) {
        A0 <- rmulti(plogis(intervention.A(L0, 0, 0, A0)))        
    }
    if (length(intervention.A0)>0){
        A0 <- rmulti(plogis(intervention.A0(L0, A0)))        
    }

    A.prev <- A0
    dN.A.prev <- rep(0, n)
    dN.L.prev <- rep(0, n)
    L.prev <- rep(0, n)
    Y.prev <- rep(0, n)
    C.prev <- rep(0, n)

    no.jumps.A.k <- rep(1, n)

    no.jumps.A <- 0
    L.star <- 1
    L.star.last <- 0

    if (verbose) {
        dt <- data.table(id=1:n, L0=L0, A0=A0, no.jumps.A=no.jumps.A)
    } else {
        dt <- data.table(id=1:n, L0=L0, A0=A0)
    }
    
    for (k in 1:K) {

        #-- generate outcome:
        
        Y1 <- Y.prev + rexpit(form.Y(L0, L.prev, A.prev, A0, no.jumps.A, dN.A.prev))*(1-C.prev)*(1-Y.prev)

        #-- generate covariates: 

        dN.L1 <- rexpit(form.dN.L(L0, dN.L.prev, L.prev, A.prev))
        L1 <- dN.L1*rexpit(form.L(L0, L.prev, A.prev, A0))+(1-dN.L1)*L.prev

        #-- generate treatment: 

        if (length(intervention.dN.A)>0) {
            dN.A1 <- rexpit(logit(intervention.dN.A(L0, dN.A.prev, L.prev, A.prev,
                                                    plogis(form.dN.A(L0, dN.A.prev, L.prev, A.prev,
                                                                     no.jumps.A, L.star)),
                                                    no.jumps.A, k, no.jumps.A.k,
                                                    L.star, k.grid, K)))
        } else {
            dN.A1 <- rexpit(form.dN.A(L0, dN.A.prev, L.prev, A.prev, no.jumps.A, L.star))
        }

        ## if (length(true.value)>0) {
        ##     A1 <- rep(true.value, n)
        ## } else

        A1 <- rmulti(plogis(form.A(L0, L.prev, A.prev, A0)))
        
        if (rescue & !ignore.rescue & length(intervention.A0)>0) {#(rescue & !ignore.rescue) {
            ## A1 <- rexpit(qlogis(0.5*plogis(form.A(L0, L.prev, A.prev, 1))+
            ##                     0.5*plogis(form.A(L0, L.prev, A.prev, 0))))
            A1 <- rmulti(0.5*plogis(form.A(L0, L.prev, A.prev, 1))+
                         0.5*plogis(form.A(L0, L.prev, A.prev, 0)))
        } else if (length(intervention.A)>0){
            A1 <- rmulti(plogis(intervention.A(L0, L.prev, A.prev, A1)))
        } 

        if (rescue & k>1) {
            A1 <- dN.A1*A1 + (1-dN.A1)*A.prev
        } else if (rescue) {
            A1 <- dN.A1*A1 + (1-dN.A1)*0
        } else {
            A1 <- dN.A1*A1 + (1-dN.A1)*A.prev
        }
        #-- generate censoring: 
        
        if (length(intervention.A)==0 & length(intervention.dN.A)==0 & length(intervention.A0)==0 & censoring) {
            C1 <- (1-Y1)*(C.prev+(1-C.prev)*rexpit(form.C(L0, L.prev, A.prev, A0)))
        } else {
            C1 <- rep(0, n)
        }
        
        no.jumps.A <- no.jumps.A+dN.A1
        if (k %in% (k.grid+1)) {
            no.jumps.A.k <- dN.A1
        } else {
            no.jumps.A.k <- no.jumps.A.k+dN.A1
        }
        L.star.last <- dN.A1*k + (1-dN.A1)*L.star.last
        L.star <- dN.A1+(1-dN.A1)*(L.star.last>=k+1-time.since)

        if (verbose) {
            dt.tmp <- data.table(Y1, dN.L1, L1, dN.A1, A1, C1, no.jumps.A, no.jumps.A.k)
            names(dt.tmp) <- paste0(c("Y", "dN.L", "L", "dN.A", "A", "C", "no.jumps.A", "no.jumps.A.k"), k)
        } else {
            dt.tmp <- data.table(Y1, dN.L1, L1, dN.A1, A1, C1)
            names(dt.tmp) <- paste0(c("Y", "dN.L", "L", "dN.A", "A", "C"), k)
        }
        
        dN.A.prev <- dN.A1
        A.prev <- A1
        dN.L.prev <- dN.L1
        L.prev <- L1
        Y.prev <- Y1
        C.prev <- C1
        
        dt <- cbind(dt, dt.tmp)

    }

    dt[, (paste0("Y", K+1)):=Y.prev + rexpit(form.Y(L0, L.prev, A.prev, A0, no.jumps.A, dN.A.prev))*(1-C.prev)*(1-Y.prev)]

    if (browse) browser()
    
    if (length(intervention.A)>0 | length(intervention.A0)>0 | length(intervention.dN.A)>0 | obs.mean) {
        return(mean(dt[, paste0("Y", K+1), with=FALSE][[1]]))
    } else {
        return(dt[1:nrow(dt)])
    }
}
