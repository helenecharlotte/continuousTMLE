cox.sl <- function(dt, tau=1.2, V=5, A.name="A", method=2, only.cox.sl=FALSE,
                   covars=c("L1", "L2", "L3"), 
                   outcome.models=list(mod1=c(Surv(time, delta==1)~A+L1+L2+L3, t0=0.9),
                                       mod2=c(Surv(time, delta==1)~A*L1.squared+L1*L2+L3, t0=NULL),
                                       mod3=c(Surv(time, delta==1)~A+L1.squared, t0=NULL),
                                       mod4=c(Surv(time, delta==1)~A+L1+L2+L3, t0=1.2),
                                       mod5=c(Surv(time, delta==1)~A+L1+L2+L3, t0=0.3),
                                       mod6=c(Surv(time, delta==1)~A+L1+L2+L3, t0=NULL),
                                       mod7=c(Surv(time, delta==1)~A+L3.squared, t0=NULL),
                                       mod8=c(Surv(time, delta==1)~A+L1+L2, t0=0.9),
                                       mod9=c(Surv(time, delta==1)~A, t0=NULL)
                                       )
                   ) {

    n <- nrow(dt)
    unique.times <- sort(unique(dt[, time]))
    cv.split <- matrix(sample(1:n, size=n), ncol=V)

    for (covar in covars)
        dt[, (paste0(covar, ".squared")):=(get(covar))^2]
    
    if (any(method==1)) {

        outlist1 <- list()
        
        #-- use nelson-aalen for baseline hazard;
        fit.baseline <- coxph(Surv(time, delta == 1)~1,
                              data=dt)
        bhaz.cox <- rbind(data.table(time=0, hazard=0),
                          data.table(time=tau, hazard=0),
                          merge(data.table(time=unique.times),
                                setDT(basehaz(fit.baseline, centered=TRUE)),
                                by="time", all.x=TRUE))
        bhaz.cox[, "dhaz0":=c(0, diff(hazard))]
        bhaz.cox[is.na(dhaz0), "dhaz0":=0]
        setnames(bhaz.cox, "hazard", "chaz0")
        bhaz.cox[, "chaz0":=cumsum(dhaz0)]
        setorder(bhaz.cox, time)
        bhaz.cox[, chaz0:=cumsum(dhaz0)]
        
        dev.loss.fun <- function(delta, Lambda) {
            M <- delta - Lambda
            #return(M^2)
            d <- sign(M)*sqrt(- 2*(M+delta*log(delta - M)))
            return(sum(na.omit(d)^2))
        }

        L2.loss.fun <- function(dN, dLambda) {
            print("use L2 loss")
            return(sum((dN - dLambda)^2))
        }

        for (vv in 1:V) {
    
            test.set <- cv.split[,vv]#sample(1:n, floor(n/10))
            train.set <- dt[, id][!dt[, id] %in% test.set]

            dt.train <- dt[id %in% train.set]
            dt.test <- dt[id %in% test.set]

            outlist1[[vv]] <- lapply(outcome.models, function(outcome.model) {
                if (length(outcome.model)>1) { #-- if there is a change-point:

                    #print("estimate time-varying hazard")

                    t0 <- outcome.model[[2]]

                    dt.train[, time.indicator:=(time<=t0)]

                    dt.train2 <- rbind(dt.train, dt.train)[order(id)]
                    dt.train2[, period:=1:.N, by="id"]

                    #dt.train2 <- dt.train2[!time.indicator | period==1]
                    dt.train2[period==1, `:=`(tstart=0, tstop=(time<=t0)*time+(time>t0)*t0)]
                    dt.train2[period==2, `:=`(tstart=t0, tstop=time)]
                    dt.train2[period==1 & !time.indicator, delta:=0]

                    mod1 <- as.character(outcome.model[[1]])
                    mod2 <- paste0(gsub(substr(mod1[2], which(strsplit(mod1[2], "")[[1]]=="(")+1,
                                               which(strsplit(mod1[2], "")[[1]]==",")-1), "tstart, tstop", mod1[2]),
                                   "~", 
                                   gsub("\\+A", "", gsub(" ", "", paste0("I((period==1)&(", A.name, "==1))", #"*L3",
                                                                         " + I((period==2)&(", A.name, "==1))", #"*L3",
                                                                         " + ",
                                                                         mod1[3]))))


                    fit.cox.vv <- coxph(formula(mod2), data=dt.train2[!time.indicator | period==1])
                    #if (t0==0.9 & vv==10)
                   # print(t0)
                   # print(fit.cox.vv)

                    dt.test[, time.indicator:=(time<=t0)]
                    
                    dt.test2 <- rbind(dt.test, dt.test)[order(id)]
                    dt.test2[, period:=1:.N, by="id"]

                    dt.test2[period==1, `:=`(tstart=0, tstop=(time<=t0)*time+(time>t0)*t0)]
                    dt.test2[period==2, `:=`(tstart=t0, tstop=time)]
                    dt.test2[period==1 & !time.indicator, delta:=0]

                    dt.test2 <- dt.test2[!(period==2 & time<=t0)]
                    dt.test2[period==1, time:=(time<=t0)*time+(time>t0)*t0]
                    #dt.test2[period==2, time:=(time<=tau)*time+(time>tau)*tau]
                    dt.test2[period==1 & time==t0, delta:=0]
                    #dt.test2[period==2 & time==tau, delta:=0]
                
                    dt.test2[, fit.lambda.cox.vv:=predict(fit.cox.vv, newdata=dt.test2,
                                                          type="risk")]
                    
                   # browser()

                    unique.times.vv <- data.table(time=dt.test[delta==1, unique(time)])

                    mat.vv <- do.call("rbind", lapply(1:nrow(dt.test), function(i) {
                        tmp <- cbind(id=dt.test[i, id], unique.times.vv)
                        tmp[, time.obs:=dt.test[id==dt.test[i, id], time]]
                        tmp[, delta.obs:=dt.test[id==dt.test[i, id], delta]]
                        tmp[, period:=(time<=t0)*1+(time>t0)*2]
                        tmp <- merge(tmp, dt.test[id==dt.test[i, id], -c("time", "delta"), with=FALSE],
                                     by="id")
                        return(tmp[time<=time.obs])
                    }))

                    mat.vv[, fit.lambda.cox.vv:=predict(fit.cox.vv, newdata=mat.vv,
                                                        type="risk")]


                    dt.test2[, fit.lambda.cox.vv:=predict(fit.cox.vv, newdata=dt.test2,
                                                          type="risk")]

                    dt.test22 <- dt.test2[rev(order(time))]
                    dt.test22[, denom:=cumsum(fit.lambda.cox.vv), by="period"]
                    dt.test22[, dHaz:=delta/denom]
                    
                    mat.vv <- merge(mat.vv, dt.test22[time!=t0 & delta>0, c("time", "dHaz"), with=FALSE],
                                    by="time")

                    mat.vv[, dN:=delta.obs*(time==time.obs)]
                    mat.vv[, dLambda:=dHaz*fit.lambda.cox.vv]
                    
                    return(L2.loss.fun(dN=mat.vv$dN, dLambda=mat.vv$dLambda))

                    
                    bhaz.cox <- unique(rbind(bhaz.cox,
                                             data.table(time=t0, chaz0=0, dhaz0=0))[order(time)],
                                       by="time")
                    bhaz.cox[, chaz0:=cumsum(dhaz0)]
                    
                    dt.test2 <- merge(dt.test2, bhaz.cox[, c("time", "dhaz0", "chaz0"), with=FALSE],
                                      by="time")[order(id)]

                    dt.test2[, chaz0.inc:=chaz0-c(0,chaz0[-.N]), by=c("id")]

                    dt.test2[, Lambda.cox.vv:=cumsum(chaz0.inc*fit.lambda.cox.vv), by=c("id")]
                    dt.test2[, Lambda.cox.vv:=predict(fit.cox.vv, newdata=dt.test2,
                                                      type="expected")]

                    dt.test22 <- dt.test2[rev(order(time))]
                    dt.test22[, denom:=cumsum(fit.lambda.cox.vv), by="period"]
                    dt.test22[, dHaz:=delta/denom]
                    dHaz.dt <- unique(dt.test22[, c("time", "dHaz", "period"), with=FALSE][order(time)])
                    dHaz.dt[, cHaz:=cumsum(dHaz), by="period"]
                    dt.test2 <-  merge(dt.test2, dHaz.dt[, -c("dHaz", "period")], by="time")
                    
                    dt.test2[, Lambda.cox.vv:=cumsum(cHaz*fit.lambda.cox.vv), by=c("id")]
                    
                    dt.test2[, N:=.N, by="id"]
                    
                    dt.test2 <- dt.test2[period==N][, -"N", with=FALSE]

                   # if ((dev.loss.fun(delta=dt.test2$delta, Lambda=dt.test2$Lambda.cox.vv))==Inf) browser()
                    
                    return(dev.loss.fun(delta=dt.test2$delta, Lambda=dt.test2$Lambda.cox.vv))
                    
                } else { #-- if there is no change-point:
                    fit.cox.vv <- coxph(as.formula(deparse(outcome.model[[1]])), #outcome.model,
                                        data=dt.train)
                    #print(outcome.model)
                 #   print(fit.cox.vv)
                    dt.test[, fit.lambda.cox.vv:=predict(fit.cox.vv, newdata=dt.test,
                                                         type="risk")]
                    #dt.test[, time:=(time<=tau)*time+(time>tau)*tau]
                    #dt.test[time==tau, delta:=0]

 unique.times.vv <- data.table(time=dt.test[delta==1, unique(time)])

                    mat.vv <- do.call("rbind", lapply(1:nrow(dt.test), function(i) {
                        tmp <- cbind(id=dt.test[i, id], unique.times.vv)
                        tmp[, time.obs:=dt.test[id==dt.test[i, id], time]]
                        tmp[, delta.obs:=dt.test[id==dt.test[i, id], delta]]
                        tmp[, period:=(time<=t0)*1+(time>t0)*2]
                        tmp <- merge(tmp, dt.test[id==dt.test[i, id], -c("time", "delta"), with=FALSE],
                                     by="id")
                        return(tmp[time<=time.obs])
                    }))

                    mat.vv[, fit.lambda.cox.vv:=predict(fit.cox.vv, newdata=mat.vv,
                                                        type="risk")]

                    dt.test[, fit.lambda.cox.vv:=predict(fit.cox.vv, newdata=dt.test,
                                                         type="risk")]

                    dt.test2 <- dt.test[rev(order(time))]
                    dt.test2[, denom:=cumsum(fit.lambda.cox.vv)]
                    dt.test2[, dHaz:=delta/denom]
                    
                    mat.vv <- merge(mat.vv, dt.test2[time!=t0 & delta>0, c("time", "dHaz"), with=FALSE],
                                    by="time")

                    mat.vv[, dN:=delta.obs*(time==time.obs)]
                    mat.vv[, dLambda:=dHaz*fit.lambda.cox.vv]
                    
                    return(L2.loss.fun(dN=mat.vv$dN, dLambda=mat.vv$dLambda))


                    
                    dt.test <- merge(dt.test, bhaz.cox[, c("time", "dhaz0", "chaz0"), with=FALSE],
                                     by="time")

                    dt.test[, Lambda.cox.vv:=chaz0*fit.lambda.cox.vv]

                    dt.test <- dt.test[rev(order(time))]
                    dt.test[, denom:=cumsum(fit.lambda.cox.vv)]
                    dt.test[, dHaz:=delta/denom]
                    dHaz.dt <- unique(dt.test[, c("time", "dHaz"), with=FALSE][order(time)])
                    dHaz.dt[, cHaz:=cumsum(dHaz)]
                    dt.test <-  merge(dt.test, dHaz.dt[, -c("dHaz")], by="time")
                    
                    dt.test[, Lambda.cox.vv:=cumsum(cHaz*fit.lambda.cox.vv), by=c("id")]
                
                    return(dev.loss.fun(delta=dt.test$delta, Lambda=dt.test$Lambda.cox.vv))
                }
            
            })

            #print(outlist1[[vv]])
        }
    }

    if (any(method==2)) {
    
        partial.loss.fun <- function(cox.fit, risk.set, t0=NULL) {
            risk.dt <- dt[id%in%risk.set][rev(order(time))]
            if (length(t0)>0) {
                risk.dt2 <- risk.dt[time>t0]
                risk.dt2[, time:=t0]
                risk.dt <- rbind(risk.dt2, risk.dt)[order(id)]
                risk.dt[, period:=1:.N, by="id"]
                risk.dt[, max.period:=.N, by="id"] 
                risk.dt <- risk.dt[rev(order(time))]

                risk.dt[, fit.lp:=predict(cox.fit, type="lp", newdata=risk.dt)]
                risk.dt[, term2:=c(1, cumsum(exp(fit.lp))[-.N]), by="period"]
                return(sum(risk.dt[delta>0, delta*(period==max.period)*(fit.lp - log(term2))]))
            } else {
                risk.dt[, fit.lp:=predict(cox.fit, type="lp", newdata=risk.dt)]
                risk.dt[, term2:=c(1, cumsum(exp(fit.lp))[-.N])]
                return(sum(risk.dt[delta>0, delta*(fit.lp - log(term2))]))
            }
        }
   
        outlist2 <- list()
    
        for (vv in 1:V) {
    
            test.set <- cv.split[,vv]#sample(1:n, floor(n/10))
            train.set <- dt[, id][!dt[, id] %in% test.set]

            outlist2[[vv]] <- lapply(outcome.models, function(outcome.model) {
                tryCatch(
                    if (length(outcome.model)==1) {
                        train.fit <- coxph(outcome.model[[1]], data=dt[id%in%train.set])
                        return(partial.loss.fun(train.fit, 1:n)-partial.loss.fun(train.fit, train.set))
                    } else {

                        t0 <- outcome.model[[2]]

                        dt.train <- dt[id%in%train.set]
                        dt.train[, time.indicator:=(time<=t0)]

                        dt2 <- rbind(dt.train, dt.train)[order(id)]
                        dt2[, period:=1:.N, by="id"]

                        dt2[period==1, `:=`(tstart=0, tstop=(time<=t0)*time+(time>t0)*t0)]
                        dt2[period==2, `:=`(tstart=t0, tstop=time)]
                        dt2[period==1 & !time.indicator, delta:=0]

                        mod1 <- as.character(outcome.model[[1]])
                        mod2 <- paste0(gsub(substr(mod1[2], which(strsplit(mod1[2], "")[[1]]=="(")+1,
                                                   which(strsplit(mod1[2], "")[[1]]==",")-1), "tstart, tstop", mod1[2]),
                                       "~", 
                                       gsub("\\+A", "", gsub(" ", "", paste0("I((period==1)&(", A.name, "==1))",
                                                                             " + I((period==2)&(", A.name, "==1))", " + ",
                                                                             mod1[3]))))

                        train.fit <- coxph(formula(mod2), data=dt2[!time.indicator | period==1])
                        #if (t0==0.9 & vv==10) print(train.fit)
                        return(partial.loss.fun(train.fit, 1:n, t0=outcome.model[[2]])-
                               partial.loss.fun(train.fit, train.set, t0= outcome.model[[2]]))
                    }
                  , error=function(e) Inf)
            })
        }
    }

    if (any(method==1)) {
        cve1 <- unlist(lapply(1:length(outcome.models), function(mm) {
            sum(unlist(lapply(outlist1, function(out) out[[mm]])))
        }))
        print(paste0("model picked by method1: ", (names(outcome.models)[cve1==min(cve1[abs(cve1)<Inf])])))
    }

    if (any(method==2)) {
        cve2 <- unlist(lapply(1:length(outcome.models), function(mm) {
            sum(-unlist(lapply(outlist2, function(out) out[[mm]])))
        }))
        print(paste0("model picked by method2: ", (names(outcome.models)[cve2==min(cve2[abs(cve2)<Inf])])))
    }

    if (length(method)==2) {
        print(cbind(cve1, cve2))
    }

    if (only.cox.sl) return(c(method1=(names(outcome.models)[cve1==min(cve1[abs(cve1)<Inf])]),
                              method2=(names(outcome.models)[cve2==min(cve2[abs(cve2)<Inf])])))

    if (any(method==1)) {
        return(names(outcome.models)[cve1==min(cve1[abs(cve1)<Inf])])
    } else {
        return(names(outcome.models)[cve2==min(cve2[abs(cve2)<Inf])])
    }
}






