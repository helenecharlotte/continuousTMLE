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

    partial.loss.fun <- function(cox.fit, test.set, risk.set, lambda.cv) {
        tmp <- copy(mat2)
        Xnew <- X[, cols.obs]
        tmp[, fit.lp:=predict(cox.fit, type="link", newx=Xnew, s=lambda.cv)]
        tmp[, risk:=0]
        tmp[id%in%risk.set, risk:=1]
        # tmp[!id%in%risk.set, risk:=1]
        tmp <- tmp[rev(order(get(time.var)))]
        tmp[, term2:=cumsum(risk*exp(fit.lp))]
        tmp[term2==0, term2:=1]
        return(sum(tmp[id%in%test.set, (delta.obs==delta.outcome)*
                                       (fit.lp - log(term2))]))        #(fit.lp - log(term2))]))
        # return(sum(tmp[risk==1, (delta.obs==delta.outcome)*
        #                        (fit.lp - (term2))]))#(fit.lp - log(term2))]))
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

        if (FALSE) for (lambda.cv in lambda.cvs)  {
                       train.fit <- glmnet(x=as.matrix(X.obs), y=Y, family="cox", maxit=1000,
                                           penalty.factor=penalty.factor, 
                                           lambda=lambda.cv)
                       mat2[id %in% train.set, (paste0("fit.lp.",lambda.cv)):=predict(train.fit, type="link", newx=X.obs, s=lambda.cv)]
                   }

        if (TRUE) {
            outlist.hal[[vv]] <- lapply(lambda.cvs, function(lambda.cv) {
                train.fit <- glmnet(x=as.matrix(X.obs), y=Y, family="cox", maxit=1000,
                                    penalty.factor=penalty.factor, 
                                    lambda=lambda.cv)
                #if (verbose) print(paste0("lambda=", lambda.cv))
                #if (verbose) print(coef(train.fit))
                mat2[id %in% train.set, (paste0("fit.lp.", lambda.cv)):=predict(train.fit, type="link", newx=X.obs, s=lambda.cv)]
                if (sum(abs(coef(train.fit,s=lambda.cv)[,1]))==0) {
                    return(-Inf)
                } else {
                    return(partial.loss.fun(train.fit, test.set, train.set, lambda.cv=lambda.cv))
                    #return(partial.loss.fun(train.fit, test.set, 1:n, lambda.cv=lambda.cv)-
                    #       partial.loss.fun(train.fit, test.set, train.set, lambda.cv=lambda.cv))
                }
            })
        }
    }

    if (FALSE) {

        Y <- mat2[, Surv(time.obs, delta.obs==delta.outcome)]
        X.obs <- X[, cols.obs]
        
        for (lambda.cv in lambda.cvs) {
            fit <- glmnet(x=as.matrix(X.obs), y=Y, family="cox", maxit=1000,
                          penalty.factor=penalty.factor, 
                          lambda=lambda.cv)
            mat2[, (paste0("fit.all.lp.exp.", lambda.cv)):=exp(predict(fit, type="link", newx=X.obs, s=lambda.cv))]
        }

        basehaz <- glmnet_basesurv(mat2[, time.obs],
                                   mat2[, delta.obs==delta.outcome],
                                   predict(fit, X.obs, type="link"), centered=TRUE)
        
        mat2[, c("fit.cox1", "fit.all.lp.exp.0.002"), with=FALSE]
        
        for (lambda.cv in lambda.cvs) {
            mat2[, (paste0("fit.lp.exp",lambda.cv)):=exp(get(paste0("fit.lp.",lambda.cv)))]
            print(mat2[, sum((delta.obs==delta.outcome)*get(paste0("fit.lp.",lambda.cv)))])
        }

        for (lambda.cv in lambda.cvs) {
            mat2[, (paste0("fit.lp.exp",lambda.cv)):=exp(get(paste0("fit.lp.",lambda.cv)))]
            mat2[, fit.lp:=get(paste0("fit.lp.",lambda.cv))]
            mat2 <- mat2[rev(order(get(time.var)))]
            mat2[, term2:=cumsum(exp(fit.lp))]
            print(mat2[, sum((delta.obs==delta.outcome)*log(term2))])
        }
        
        test.error <- lapply(lambda.cvs, function(lambda.cv) {
            mat2[, fit.lp:=get(paste0("fit.lp.",lambda.cv))]
            mat2 <- mat2[rev(order(get(time.var)))]
            mat2[, term2:=cumsum(exp(fit.lp))]
            return(mat2[, -sum((delta.obs==delta.outcome)*(fit.lp - log(term2)))])
        })

        (list(lambda.cv=lambda.cvs[unlist(test.error)==min(unlist(test.error))][1],
              cve=min(unlist(test.error))))
    }
    

    #-- remove if all failed to run in given split? 
    outlist.hal <- lapply(outlist.hal, function(out) {
        if (all(unique(unlist(out))==Inf)) return(NULL) else return(out)
        #message("obs: HAL CV failed in at least one of the sample splits")
    })

    cve.hal <- unlist(lapply(1:length(lambda.cvs), function(mm) {
        sum(-unlist(lapply(outlist.hal, function(out) out[[mm]])))
    }))

    if (verbose) print(cbind(lambda.cvs,cve.hal))

    if (all(unique(cve.hal)==Inf)) {
        return(NULL)
    } else {
        return(list(lambda.cv=lambda.cvs[cve.hal==min(cve.hal)][1],
                    cve=min(cve.hal)))
        #return(list(lambda.cv=lambda.cvs[unlist(test.error)==min(unlist(test.error))][1],
        #            cve=min(unlist(test.error))))
    }
}
