cox.hal.sl <- function(mat2, dt, X=NULL, time.var="time", A.name="A",
                       delta.outcome=1, cols.obs=cols.obs,
                       cut.covars=5, cut.time=10, cut.time.A=4,
                       cut.L.A=5, cut.L.interaction=5,
                       covars=c("L1", "L2", "L3"),
                       lambda.cv=NULL,
                       penalize.time=TRUE, adjust.penalization=TRUE,
                       lambda.cvs=seq(0, 0.003, length=51)[-1], 
                       V=10, verbose=TRUE
                       ) {

    set.seed(19192)

    n <- nrow(dt)
    unique.times <- sort(unique(dt[, get(time.var)]))
    cv.split <- matrix(sample(1:n, size=n), ncol=V)

    outlist.hal <- list()

    partial.loss.fun <- function(cox.fit, risk.set) {
        tmp <- copy(mat2)
        Xnew <- X[, cols.obs]
        tmp[, fit.lp:=predict(cox.fit, type="link", newx=Xnew)]
        tmp[, risk:=0]
        tmp[id%in%risk.set, risk:=1]
        tmp <- tmp[rev(order(get(time.var)))]
        tmp[, term2:=cumsum(risk*exp(fit.lp))]
        tmp[term2==0, term2:=1]
        return(sum(tmp[risk==0, (delta.obs==delta.outcome)*
                                (fit.lp - log(term2))]))
        #test.dt <- mat2[!id%in%risk.set]
        #risk.dt <- mat2[id%in%risk.set]#[rev(order(time))]
        #Xnew <- X[mat2$id %in% risk.set,][, cols.obs]#[rev(order(risk.dt$time)),]

        
        
        #risk.dt[, fit.lp:=predict(cox.fit, type="link", newx=Xnew)]
        #risk.dt <- risk.dt[rev(order(time))]
        #risk.dt[, term2:=c(1, cumsum(exp(fit.lp))[-.N])]
        #return(sum(risk.dt[(delta.obs==delta.outcome)>0, (delta.obs==delta.outcome)*
        #                                                 (fit.lp - log(term2))]))
    }

    for (vv in 1:V) {
    
        test.set <- cv.split[,vv]#sample(1:n, floor(n/10))
        train.set <- dt[, id][!dt[, id] %in% test.set]

        mat2.train <- mat2[id %in% train.set]

        Y <- mat2.train[, Surv(time.obs, delta.obs==delta.outcome)]
        X.obs <- X[mat2$id %in% train.set,][, cols.obs]

        penalty.factor <- rep(1, ncol(X.obs))
        penalty.factor[1] <- 0
        
        #----------------------------
        #---- 2. compute cve for hal
        #----------------------------

        outlist.hal[[vv]] <- lapply(lambda.cvs, function(lambda.cv) {
            train.fit <- glmnet(x=as.matrix(X.obs), y=Y, family="cox", maxit=1000,
                                penalty.factor=penalty.factor, 
                                lambda=lambda.cv)
            if (sum(abs(coef(train.fit)[,1]))==0) {
                return(-Inf)
            } else {
                return(partial.loss.fun(train.fit, train.set))
            }
        })
    }

    #-- remove if all failed to run in given split? 
    outlist.hal <- lapply(outlist.hal, function(out) {
        if (all(unique(unlist(out))==Inf)) return(NULL) else return(out)
        message("obs: HAL cv failed in at least one of the sample splits")
    })

    cve.hal <- unlist(lapply(1:length(lambda.cvs), function(mm) {
        sum(-unlist(lapply(outlist.hal, function(out) out[[mm]])))
    }))

    if (verbose) print(cbind(lambda.cvs,cve.hal))

    if (all(unique(cve.hal)==Inf)) {
        return(NULL)
    } else {
        return(lambda.cvs[cve.hal==min(cve.hal)])
    }
}
