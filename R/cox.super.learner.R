cox.sl <- function(dt, V=5, A.name="A", method=2, only.cox.sl=FALSE,
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

    if (any(method>=1)) {
    
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
                if (any(class(cox.fit)=="coxnet")) {
                    Xnew <- as.matrix(model.matrix(outcome.model, data=risk.dt))
                    risk.dt[, fit.lp:=predict(cox.fit, type="link", newx=Xnew)]
                } else {
                    risk.dt[, fit.lp:=predict(cox.fit, type="lp", newdata=risk.dt)]
                }
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
                    if (length(outcome.model)==1 | names(outcome.model)[1]=="coxnet") {
                        if (any(names(outcome.model)=="coxnet")) {
                            X <- model.matrix(outcome.model[[1]], data=dt[id%in%train.set])
                            y <- dt[id%in%train.set, Surv(time, delta==1)]
                            train.fit <- glmnet(x=X, y=y, family="cox", maxit=1000,
                                                lambda=outcome.model[[2]])
                        } else {
                            train.fit <- coxph(outcome.model[[1]], data=dt[id%in%train.set])
                        }
                        if (method==2) {
                            return(partial.loss.fun(train.fit, 1:n)-partial.loss.fun(train.fit, train.set))
                        } else if (method==3) {
                            return(partial.loss.fun(train.fit, train.set))
                        }
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
                        if (method==2) {
                            return(partial.loss.fun(train.fit, 1:n, t0=outcome.model[[2]])-
                                   partial.loss.fun(train.fit, train.set, t0=outcome.model[[2]]))
                        } else if (method==3) {
                            return(partial.loss.fun(train.fit, train.set, t0=outcome.model[[2]]))
                        }
                    }
                  , error=function(e) Inf)
            })
        }
    }

    if (any(method>=1)) {
        cve2 <- unlist(lapply(1:length(outcome.models), function(mm) {
            sum(-unlist(lapply(outlist2, function(out) out[[mm]])))
        }))
    }
    
    if (only.cox.sl) {
        print(cbind(cve2))
        return(c(#method1=(names(outcome.models)[cve1==min(cve1[abs(cve1)<Inf])]),
            method2=(names(outcome.models)[cve2==min(cve2[abs(cve2)<Inf])])))
    }

    return(names(outcome.models)[cve2==min(cve2[abs(cve2)<Inf])])
}






