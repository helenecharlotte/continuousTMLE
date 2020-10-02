poisson.hal.sl <- function(mat, dt, outcome.model, change.point, X=NULL, delta.outcome=1, cols.obs=cols.obs,
                           cut.covars=5, cut.time=10, cut.time.A=4,
                           cut.L1.A=5, cut.L.interaction=5,
                           covars=c("L1", "L2", "L3"),
                           poisson.cens=FALSE, lambda.cv=NULL,
                           penalize.time=TRUE, adjust.penalization=TRUE,
                           lambda.cvs=seq(0, 0.003, length=51)[-1], 
                           V=10
                           ) {

    n <- nrow(dt)
    unique.times <- sort(unique(dt[, time]))
    unique.T <- sort(unique(dt[delta==1, time]))
    unique.Tdiff <- unique.T - c(0, unique.T[-length(unique.T)])
    cv.split <- matrix(sample(1:n, size=n), ncol=V)

    outlist.cox <- list()
    outlist.hal <- list()

    log.like.loss.fun <- function(dN, dLambda, Lambda) {
        if (any(dN==0 & dLambda==0)) {
            dLambda[dN==0 & dLambda==0] <- 1
        }
        return(-sum(log(dLambda)*dN - exp(-Lambda)))
    }

    lebesgue.log.like.loss.fun <- function(dN, lambda, Lambda) {
        if (any(dN==0 & lambda==0)) {
            lambda[dN==0 & lambda==0] <- 1
        }
        return(-sum(log(lambda)*dN - exp(-Lambda)))
    }

    for (vv in 1:V) {
    
        test.set <- cv.split[,vv]#sample(1:n, floor(n/10))
        train.set <- dt[, id][!dt[, id] %in% test.set]

        #----------------------------
        #---- 1. compute cve for cox
        #----------------------------

        if (!length(lambda.cvs)>0) {

            dt.train <- dt[id %in% train.set]
            dt.test <- dt[id %in% test.set]

            if (length(change.point)>1) { #-- if there is a change-point:

                t0 <- change.point

                dt.train[, time.indicator:=(time<=t0)]

                dt.train2 <- rbind(dt.train, dt.train)[order(id)]
                dt.train2[, period:=1:.N, by="id"]

                dt.train2[period==1, `:=`(tstart=0, tstop=(time<=t0)*time+(time>t0)*t0)]
                dt.train2[period==2, `:=`(tstart=t0, tstop=time)]
                dt.train2[period==1 & !time.indicator, delta:=0]

                mod1 <- as.character(outcome.model)
                mod2 <- paste0(gsub(substr(mod1[2], which(strsplit(mod1[2], "")[[1]]=="(")+1,
                                           which(strsplit(mod1[2], "")[[1]]==",")-1), "tstart, tstop", mod1[2]),
                               "~", 
                               gsub("\\+A", "", gsub(" ", "", paste0("I((period==1)&(", A.name, "==1))", #"*L3",
                                                                     " + I((period==2)&(", A.name, "==1))", #"*L3",
                                                                     " + ",
                                                                     mod1[3]))))

                fit.cox.vv <- coxph(formula(mod2), data=dt.train2[!time.indicator | period==1])

                dt.test[, time.indicator:=(time<=t0)]
                    
                dt.test2 <- rbind(dt.test, dt.test)[order(id)]
                dt.test2[, period:=1:.N, by="id"]

                dt.test2[period==1, `:=`(tstart=0, tstop=(time<=t0)*time+(time>t0)*t0)]
                dt.test2[period==2, `:=`(tstart=t0, tstop=time)]
                dt.test2[period==1 & !time.indicator, delta:=0]

                dt.test2 <- dt.test2[!(period==2 & time<=t0)]
                dt.test2[period==1, time:=(time<=t0)*time+(time>t0)*t0]
                dt.test2[period==1 & time==t0, delta:=0]
                
                dt.test2[, fit.lambda.cox.vv:=predict(fit.cox.vv, newdata=dt.test2,
                                                      type="risk")]
                    
                dt.test22 <- dt.test2[rev(order(time))]
                dt.test22[, denom:=cumsum(fit.lambda.cox.vv), by="period"]
                dt.test22[, dHaz:=delta/c(1, denom[-.N])]

                setorder(dt.test22, time)
                dt.test22[, chaz:=cumsum(dHaz)]
                    
                dt.test2 <- merge(dt.test2, dt.test22[, c("id", "time", "chaz", "dHaz"),
                                                      with=FALSE], by=c("id", "time"))

                dt.test2[, Lambda:=(chaz*fit.lambda.cox.vv), by=c("id")]
                dt.test2[, dLambda:=dHaz*fit.lambda.cox.vv, by=c("id")]
                dt.test2[, dN:=delta]
                dt.test2[dN==0, dLambda:=1]
                    
                dt.test2[, idN:=1:.N, by="id"]
                dt.test2[, N:=.N, by="id"]

                outlist.cox[[vv]] <- log.like.loss.fun(dN=dt.test2[idN==N]$dN,
                                                       dLambda=dt.test2[idN==N]$dLambda,
                                                       Lambda=dt.test2[idN==N]$Lambda)
            } else { #-- if there is no change-point:
                fit.cox.vv <- coxph(as.formula(deparse(outcome.model)),
                                    data=dt.train)
                dt.test[, fit.lambda.cox.vv:=predict(fit.cox.vv, newdata=dt.test,
                                                     type="risk")]

                dt.test2 <- dt.test[rev(order(time))]
                dt.test2[, denom:=cumsum(fit.lambda.cox.vv)]
                dt.test2[, dHaz:=delta/c(1, denom[-.N])]#[, dHaz:=delta/denom]

                setorder(dt.test2, time)
                dt.test2[, chaz:=cumsum(dHaz)]
                    
                dt.test <- merge(dt.test, dt.test2[, c("id", "time", "chaz", "dHaz"),
                                                   with=FALSE], by=c("id", "time"))
                    
                dt.test[, Lambda:=(chaz*fit.lambda.cox.vv), by=c("id")]
                dt.test[, dLambda:=dHaz*fit.lambda.cox.vv, by=c("id")]
                dt.test[, dN:=delta]
                dt.test[dN==0, dLambda:=1]
                    
                dt.test[, idN:=1:.N, by="id"]
                dt.test[, N:=.N, by="id"]

                outlist.cox[[vv]] <- log.like.loss.fun(dN=dt.test[idN==N]$dN,
                                                       dLambda=dt.test[idN==N]$dLambda,
                                                       Lambda=dt.test[idN==N]$Lambda)

            }
        }

        #----------------------------
        #---- 2. compute cve for hal
        #----------------------------

        print(paste0("v=", vv, "/", V))

        mat.train <- mat[id %in% train.set]

        mat.train[, RT:=sum(tdiff*(time<=time.obs)), by=c("x", "A")]
        mat.train[, D:=sum(event*(time<=time.obs)), by=c("x", "A")]

        if (length(mat[, unique(A)])>1) reduce.A <- TRUE else reduce.A <- FALSE

        tmp.train <- unique(mat.train[time<=time.obs & (!reduce.A | A.obs==A), c("RT", "D", "x"), with=FALSE])

        X.obs <- unique.matrix(X[mat$id %in% train.set,][mat.train$time<=mat.train$time.obs,
                                                         cols.obs])[tmp.train$RT>0,]

        Y <- tmp.train[RT>0, D] 
        offset <- tmp.train[RT>0, log(RT)]

        if (penalize.time) {
            penalty.factor <- c(0, rep(1, ncol(X.obs)-1))
        } else {
            #print("no penalization of coefficients for time indicators (main effects)")
            penalty.factor <- rep(1, ncol(X.obs))
            penalty.factor[1:cut.time] <- 0
        }

        if (length(lambda.cvs)>0) {
            outlist.hal[[vv]] <- lapply(lambda.cvs, function(lambda.cv) {
                fit.vv <- glmnet(x=X.obs, y=Y, lambda=lambda.cv,#lambda.cv,
                                 family="poisson",
                                 offset=offset,
                                 #penalty.factor=c(0, rep(1, ncol(X.obs) - 2)),
                                 penalty.factor=penalty.factor,
                                 maxit=1000)
                mat[id %in% test.set, fit.lambda.vv:=exp(predict(fit.vv, X[mat$id %in% test.set, cols.obs],
                                                                 newoffset=0))]

                mat[id %in% test.set, fit.pois.dLambda.vv:=fit.lambda.vv*tdiff]
                mat[id %in% test.set, fit.pois.Lambda.vv:=cumsum(fit.pois.dLambda.vv), by=c("id", "A")]

                return(lebesgue.log.like.loss.fun(dN=mat[id %in% test.set & time==time.obs, delta.obs],
                                                  lambda=mat[id %in% test.set & time==time.obs, fit.lambda.vv],
                                                  Lambda=mat[id %in% test.set & time==time.obs, fit.pois.Lambda.vv]))
            })
        } else {

            if (length(lambda.cv)>0) {
                fit.vv <- glmnet(x=X.obs, y=Y, lambda=0,#lambda.cv,
                                 family="poisson",
                                 offset=offset,
                                 #penalty.factor=c(0, rep(1, ncol(X.obs) - 2)),
                                 penalty.factor=penalty.factor,
                                 maxit=1000)
            } else {
                print("CV for penalization")
                fit.vv <- cv.glmnet(x=X.obs, y=Y, 
                                    family="poisson",
                                    offset=offset,
                                    #penalty.factor=c(0, rep(1, ncol(X.obs) - 2)),
                                    penalty.factor=penalty.factor,
                                    maxit=1000)
                if (adjust.penalization) {
                    #-- this was an extra ad hoc step to not penalize "too much"
                    fit.vv <- glmnet(x=X.obs, y=Y, 
                                     family="poisson",
                                     offset=offset, lambda=fit.vv$lambda.1se*0.2,
                                     #penalty.factor=c(0, rep(1, ncol(X.obs) - 2)),
                                     penalty.factor=penalty.factor,
                                     maxit=1000)
                }
                coef(fit.vv)
            }
        
            mat[id %in% test.set, fit.lambda.vv:=exp(predict(fit.vv, X[mat$id %in% test.set, cols.obs],
                                                             newoffset=0))]

            mat[id %in% test.set, fit.pois.dLambda.vv:=fit.lambda.vv*tdiff]
            mat[id %in% test.set, fit.pois.Lambda.vv:=cumsum(fit.pois.dLambda.vv), by=c("id", "A")]

            outlist.hal[[vv]] <- lebesgue.log.like.loss.fun(dN=mat[id %in% test.set & time==time.obs, delta.obs],
                                                            lambda=mat[id %in% test.set & time==time.obs, fit.lambda.vv],
                                                            Lambda=mat[id %in% test.set & time==time.obs, fit.pois.Lambda.vv])

            #outlist.hal[[vv]] <- log.like.loss.fun(dN=mat[id %in% test.set & time==time.obs, delta.obs],
            #                                       dLambda=mat[id %in% test.set & time==time.obs, fit.pois.dLambda.vv],
            #                                       Lambda=mat[id %in% test.set & time==time.obs, fit.pois.dLambda.vv])
        }

    }

    cve.cox <- sum(unlist(outlist.cox))
    if (length(lambda.cvs)>0) {
        cve.hal <- unlist(lapply(1:length(lambda.cvs), function(mm) {
            sum(unlist(lapply(outlist.hal, function(out) out[[mm]])))
        }))
        plot(lambda.cvs,cve.hal)
        return(lambda.cvs[cve.hal==min(cve.hal)])
    } else {
        cve.hal <- sum(unlist(outlist.hal))
        return(cve.cox>cox.hal)
    }
}
