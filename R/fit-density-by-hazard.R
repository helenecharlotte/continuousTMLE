numextract <- function(string){ 
    as.numeric(str_extract(string, "\\-*\\d+\\.*\\d*"))
}

if (FALSE) {
    fit.density(dt, "A0", "L0")
    fit.density(dt, "A1", c("L0", "L1", "A0"), subset="dN.A1")
}

fit.density <- function(dt, Avar, covars, subset=NULL, test=FALSE) {

    ## dt.covars <- setDT(expand.grid(lapply(covars, function(var) {
    ##     dt[, unique(.SD), .SDcols=var][[1]]
    ## })))

    ## names(dt.covars) <- covars

    if (length(subset)>0) {
        dt.tmp <- dt[get(subset)==1, c(Avar, covars), with=FALSE]
    } else {
        dt.tmp <- dt[, c(Avar, covars), with=FALSE]
    }
    
    dt.tmp[, id:=1:nrow(dt.tmp)]

    avals <- sort(unique(dt.tmp[, Avar, with=FALSE][[1]]))
    size.A <- length(avals)-1
    
    for (a in avals) {
        dt.tmp[, (paste0("Avar", a)):=1*(a<=get(Avar))]
    }

    A.melt <- melt(dt.tmp, id=c("id", Avar, covars))[order(id)]

    A.melt[, a:=numextract(variable)]
    A.melt <- A.melt[, -"variable"]

    ## form <- formula(paste0("Y~a+", paste0(covars, collapse="+")))

    ## #dt.covars <- unique(dt.tmp[, covars, with=FALSE])
    
    ## hazard.A <- matrix(1, nrow=nrow(unique(dt.covars)), ncol=size.A+1)#rep(1, size.A+1)
    ## for (a1 in 0:(size.A-1)) {
    ##     A.melt[, Y:=(get(Avar)==a1)]
    ##     dt.covars[, a:=a1]
    ##     fit.A <- glm(form, data=A.melt[value==1 & a>=a1], family=binomial())
    ##     hazard.A[,a1+1] <- predict(fit.A, newdata=dt.covars, type="response")
    ## }

    if (TRUE) {
        A.melt[, atrisk:=1*(get(Avar)>=a)]
        A.melt[, Y:=1*(get(Avar)==a)]
        A.melt[, a:=factor(a)]
        
        form <- paste0("Y~", paste0(covars, "*a", collapse="+"))
        fit.A <- glm(formula(form), data=A.melt[atrisk==1], family=binomial())

        dt.covars <- setDT(expand.grid(lapply(c(covars, Avar), function(var) {
            dt[, unique(.SD), .SDcols=var][[1]]
        })))

        names(dt.covars) <- c(covars, Avar)
        dt.covars[, a:=factor(get(Avar))]

        #for (a1 in 0:(size.A-1)) {
        dt.covars[, hazard.A:=predict(fit.A, newdata=dt.covars, type="response")]
        #}

        dt.covars <- dt.covars[order(L0, get(Avar))]
        dt.covars[, haz.cumprod:=cumprod(1-hazard.A), by=covars]
        dt.covars[, est:=c(1,haz.cumprod[-.N])*hazard.A, by=covars]

        setnames(dt.covars, "est", paste0("fit.", Avar))
        return(dt.covars[, -c("a", "hazard.A", "haz.cumprod")])
        
    } else {

        out <- do.call("rbind", lapply(0:size.A, function(a) {
            if (a==0) est <- hazard.A[, 1] else est <- apply(hazard.A[, 1:a, drop=FALSE], 1, function(xx) {
                prod(1-xx)
            })*hazard.A[, a+1]
            data.table(dt.covars[, -"a", with=FALSE], Avar=a, fit=est)
        }))

        setnames(out, c("Avar", "fit"), c(Avar, paste0("fit.", Avar)))
        return(out)  
    }
}


if (FALSE) {

    library(data.table)
    library(stringr)
        
    n <- 1000
    L0 <- rbinom(n, 1, 0.5)
    L1 <- rbinom(n, 1, 0.5)
    #L0 <- runif(n, -3, 3)
        
    rexpit <- function(x) rbinom(n=n, size=1, prob=plogis(x))
    rmulti <- function(x) {
        out <- rep(0, n)
        for (p in 1:ncol(x)) {
            out <- p*rexpit(qlogis(x[,p]))*(out==0) + out 
        }
        out[out==0] <- ncol(x)+1
        return(out-1)
    }

    A <- rmulti(plogis(cbind(L0*0.2-0.2, 0.4*L1-0.4)))

    dt <- data.table(L0, L1, A)

    
    (fit <- fit.density(dt, "A", c("L0", "L1"), subset=NULL))
    #fit.density(dt, "A", "L0", subset=NULL)[order(L0)]    

    dt[L0==1, table(A)/nrow(dt[L0==1])]
    dt[L0==0 & L1==0, table(A)/nrow(dt[L0==0 & L1==0])]
    fit[L0==0 & L1==0]
    dt[L0==0 & L1==1, table(A)/nrow(dt[L0==0 & L1==1])]
    fit[L0==0 & L1==1]
    dt[L0==1 & L1==1, table(A)/nrow(dt[L0==1 & L1==1])]
    fit[L0==1 & L1==1]


    if (FALSE) {

        n <- 1000
        L0 <- rbinom(n, 1, 0.5)
        L1 <- rbinom(n, 1, 0.5)
        L0 <- runif(n, -3, 3)
    
        A <- rmulti(plogis(cbind(L0*0.2-0.2, 0.4*L1-0.4)))
        dt <- data.table(L0, L1, A)
    }

    
}



if (FALSE) { #-- continuous covariates

    library(data.table)
    library(stringr)

       
    n <- 1e4
    L0 <- rbinom(n, 1, 0.5)
    L1 <- rbinom(n, 1, 0.5)
    L0 <- runif(n, -3, 3)
   
    
    rexpit <- function(x) rbinom(n=n, size=1, prob=plogis(x))
    rmulti <- function(x) {
        out <- rep(0, n)
        for (p in 1:ncol(x)) {
            out <- p*rexpit(qlogis(x[,p]))*(out==0) + out 
        }
        out[out==0] <- ncol(x)+1
        return(out-1)
    }

    A <- rmulti(plogis(cbind(L0*0.2-0.2, 0.4*L1-0.4)))

    dt <- data.table(L0, L1, A)

    dt[, A.true:=do.call("rbind", lapply(1:nrow(dt), function(ii) {
        L0 <- dt[ii, L0]
        L1 <- dt[ii, L1]
        a <- dt[ii, A]
        return((as.numeric(table(A <- rmulti(plogis(cbind(L0*0.2-0.2, 0.4*L1-0.4)))))/n)[a+1])
    }))]

    (fit <- fit.density(dt, "A", c("L0", "L1"), subset=NULL))
    #fit.density(dt, "A", "L0", subset=NULL)[order(L0)]

    merge(dt, fit, by=c("L0", "L1", "A"))
    

    dt[L0==1, table(A)/nrow(dt[L0==1])]
    dt[L0==0 & L1==0, table(A)/nrow(dt[L0==0 & L1==0])]
    fit[L0==0 & L1==0]
    dt[L0==0 & L1==1, table(A)/nrow(dt[L0==0 & L1==1])]
    fit[L0==0 & L1==1]
    dt[L0==1 & L1==1, table(A)/nrow(dt[L0==1 & L1==1])]
    fit[L0==1 & L1==1]


    if (FALSE) {

        n <- 1000
        L0 <- rbinom(n, 1, 0.5)
        L1 <- rbinom(n, 1, 0.5)
        L0 <- runif(n, -3, 3)
    
        A <- rmulti(plogis(cbind(L0*0.2-0.2, 0.4*L1-0.4)))
        dt <- data.table(L0, L1, A)
    }

    
}


if (FALSE) { #-- continuous A??

    library(data.table)
    library(stringr)
       
    n <- 1e3
    L0 <- rbinom(n, 1, 0.5)
    L1 <- rbinom(n, 1, 0.5)
    L0 <- runif(n, -3, 3)

    Acont <- rnorm(n, mean=L0, sd=0.5)

    #agrid <- seq(min(Acont), max(Acont), length=10)
    agrid <- unique(sort(c(min(Acont)-0.1, sample(Acont, size=40), max(Acont)+0.1)))

    diff <- diff(agrid)

    A <- findInterval(Acont, agrid)

    rexpit <- function(x) rbinom(n=n, size=1, prob=plogis(x))
    rmulti <- function(x) {
        out <- rep(0, n)
        for (p in 1:ncol(x)) {
            out <- p*rexpit(qlogis(x[,p]))*(out==0) + out 
        }
        out[out==0] <- ncol(x)+1
        return(out-1)
    }

    #A <- rmulti(plogis(cbind(L0*0.2-0.2, 0.4*L1-0.4)))
    
    dt <- data.table(L0, L1, A)

    if (FALSE) {
        A.true <- do.call("rbind", lapply(1:nrow(dt), function(ii) {
            L0 <- dt[ii, L0]
            A <- rnorm(n, mean=L0, sd=0.5)
            return(A)
        }))

        hist(A.true[2,])
        abline(v=L0[2])
    
    }



    (fit <- fit.density(dt, "A", c("L0"), subset=NULL))
    #fit.density(dt, "A", "L0", subset=NULL)[order(L0)]

    fit[, Aval:=agrid[A]]

    fit[, fit.A.final:=fit.A/diff[A]]
    fit[, test:=diff[A]]


    par(mfrow=c(1,2))
    plot(fit[L0==0, Aval], fit[L0==0, fit.A.final])
    abline(v=0, col="red")
    lines(fit[L0==0, Aval], dnorm(fit[L0==0, Aval], mean=0, sd=0.5), type="l")
    plot(fit[L0==1, Aval], fit[L0==1, fit.A.final])
    abline(v=1, col="red")
    lines(fit[L0==1, Aval], dnorm(fit[L0==1, Aval], mean=1, sd=0.5), type="l")

    par(mfrow=c(2,4))
    plot(fit[L0==fit[5000,L0], agrid[A]], fit[L0==fit[5000,L0], fit.A.final])
    abline(v=fit[5000,L0], col="red")
    lines(fit[L0==fit[5000,L0], agrid[A]], dnorm(fit[L0==fit[5000,L0], agrid[A]], mean=fit[5000,L0], sd=0.5), type="l")
    plot(fit[L0==fit[7500,L0], agrid[A]], fit[L0==fit[7500,L0], fit.A.final])
    abline(v=fit[7500,L0], col="red")
    lines(fit[L0==fit[7500,L0], agrid[A]], dnorm(fit[L0==fit[7500,L0], agrid[A]], mean=fit[7500,L0], sd=0.5), type="l")
    plot(fit[L0==fit[12000,L0], agrid[A]], fit[L0==fit[12000,L0], fit.A.final])
    abline(v=fit[12000,L0], col="red")
    lines(fit[L0==fit[12000,L0], agrid[A]], dnorm(fit[L0==fit[12000,L0], agrid[A]], mean=fit[12000,L0], sd=0.5), type="l")
    plot(fit[L0==fit[25000,L0], agrid[A]], fit[L0==fit[25000,L0], fit.A.final])
    abline(v=fit[25000,L0], col="red")
    lines(fit[L0==fit[25000,L0], agrid[A]], dnorm(fit[L0==fit[25000,L0], agrid[A]], mean=fit[25000,L0], sd=0.5), type="l")

    #--------------------
    ##--- try with cophh? 
    library(survival)

    dt[, D:=1]
    fit.cox <- coxph(Surv(A, D) ~ L0, data=dt)

    bhaz <- setDT(basehaz(fit.cox))#[,1]
    setnames(bhaz, c("time", "hazard"), c("A", "bhaz"))
    bhaz[, dbhaz:=c(diff(bhaz), 0)]
    #dbhaz <- diff(bhaz[, "hazard"])

    #dt <- merge(bhaz, dt, by="A")

    #dt[, haz:=NULL]
    #dt[, dbhaz:=c(dbhaz, 0)]
    #dt[, fit.density:=exp(-bhaz*predict(fit.cox, type="risk"))*predict(fit.cox, type="risk")]

    fit2 <- merge(fit, bhaz, by="A")[order(L0,A)]
    fit2[, fit.cox:=as.numeric(exp(-bhaz*predict(fit.cox, type="risk", newdata=fit2))*
               dbhaz*predict(fit.cox, type="risk", newdata=fit2))]

    #par(mfrow=c(1,4))
    plot(fit2[L0==fit2[5000,L0], agrid[A]], fit2[L0==fit2[5000,L0], fit.cox])
    abline(v=fit2[5000,L0], col="red")
    lines(fit2[L0==fit2[5000,L0], agrid[A]], dnorm(fit2[L0==fit2[5000,L0], agrid[A]], mean=fit2[5000,L0], sd=0.5), type="l")
    plot(fit2[L0==fit2[7500,L0], agrid[A]], fit2[L0==fit2[7500,L0], fit.cox])
    abline(v=fit2[7500,L0], col="red")
    lines(fit2[L0==fit2[7500,L0], agrid[A]], dnorm(fit2[L0==fit2[7500,L0], agrid[A]], mean=fit2[7500,L0], sd=0.5), type="l")
    plot(fit2[L0==fit2[12000,L0], agrid[A]], fit2[L0==fit2[12000,L0], fit.cox])
    abline(v=fit2[12000,L0], col="red")
    lines(fit2[L0==fit2[12000,L0], agrid[A]], dnorm(fit2[L0==fit2[12000,L0], agrid[A]], mean=fit2[12000,L0], sd=0.5), type="l")
    plot(fit2[L0==fit2[25000,L0], agrid[A]], fit2[L0==fit2[25000,L0], fit.cox])
    abline(v=fit2[25000,L0], col="red")
    lines(fit2[L0==fit2[25000,L0], agrid[A]], dnorm(fit2[L0==fit2[25000,L0], agrid[A]], mean=fit2[25000,L0], sd=0.5), type="l")


    #--------------------------
    ##--- try now with poisson?

    # ........

    gridsize <- length(agrid)
    
    tmp <- data.table(id=rep(1:n, each=gridsize),
                      #W=rep(W, each=gridsize),
                      Aint1=rep(1:gridsize, length=n*gridsize),
                      Agrid=rep(agrid, length=n*gridsize))

    dt[, id:=1:.N]

    dt1 <- merge(dt, tmp, by="id")[Aint1<=A+1]
    dt1[A<=Agrid, Agrid:=A]

    dt1[, Adiff:=c(0, diff(Agrid)), by=id]

    dt1[, Aval:=agrid[A]]

    dt1[, RT:=sum(Adiff), by=c("A", "L0")]
    dt1[, D:=sum(A==Aint1+1), by=c("A", "L0")]

    dt.pois <- unique(dt1[, c("A", "L0", "RT", "D"), with=FALSE])

    dt.pois[, table(D)]
    dt.pois[, range(RT)]

    par(mfrow=c(1,2)) #-- noget her ser i hvert fald meget rigtigt ud! 
    plot(dt.pois[L0==0, Aval], dt.pois[L0==0, D])
    plot(dt.pois[L0==1, Aval], dt.pois[L0==1, D])
    par(mfrow=c(1,1))

    plot(dt.pois[, A], dt.pois[, D])

    
    if (FALSE) {
        merge(dt, fit, by=c("L0", "L1", "A"))
    

        dt[L0==1, table(A)/nrow(dt[L0==1])]
        dt[L0==0 & L1==0, table(A)/nrow(dt[L0==0 & L1==0])]
        fit[L0==0 & L1==0]
        dt[L0==0 & L1==1, table(A)/nrow(dt[L0==0 & L1==1])]
        fit[L0==0 & L1==1]
        dt[L0==1 & L1==1, table(A)/nrow(dt[L0==1 & L1==1])]
        fit[L0==1 & L1==1]
    }

    

    

    if (FALSE) {

        n <- 1000
        L0 <- rbinom(n, 1, 0.5)
        L1 <- rbinom(n, 1, 0.5)
        L0 <- runif(n, -3, 3)
    
        A <- rmulti(plogis(cbind(L0*0.2-0.2, 0.4*L1-0.4)))
        dt <- data.table(L0, L1, A)
    }

    
}


