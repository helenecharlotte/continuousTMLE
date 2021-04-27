cox.sl <- function(dt, V=5, A.name="A", method=3, only.cox.sl=FALSE, time.var="time", delta.var="delta",
                   verbose=FALSE,
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

    set.seed(19192)

    n <- nrow(dt)
    unique.times <- sort(unique(dt[, get(time.var)]))
    cv.split <- matrix(sample(1:n, size=n), ncol=V)

    if (any(method>=1)) {
    
        partial.loss.fun <- function(cox.fit, risk.set, t0=NULL) {
            risk.dt <- dt[id%in%risk.set][rev(order(time))]
            tmp <- copy(dt)
            tmp[, risk:=0]
            tmp[id%in%risk.set, risk:=1]
            if (length(t0)>0) {
                tmp2 <- tmp[get(time.var)>t0]
                tmp2[, (time.var):=t0]
                tmp <- rbind(tmp2, tmp)[order(id)]
                tmp[, period:=1:.N, by="id"]
                tmp[, max.period:=.N, by="id"]
                tmp <- tmp[rev(order(get(time.var)))]
                tmp[, fit.lp:=predict(cox.fit, type="lp", newdata=tmp)]
                tmp[, term2:=cumsum(risk*exp(fit.lp)), by="period"]
                tmp[term2==0, term2:=1]
                #return(sum(tmp[risk==0, get(delta.var)*(period==max.period)*(fit.lp - log(term2))]))
                return(sum(tmp[risk==0, get(delta.var)*(period==max.period)*(fit.lp - log(term2))]))
            } else {
                if (any(class(cox.fit)=="coxnet")) {
                    Xnew <- as.matrix(model.matrix(outcome.model, data=tmp))
                    tmp[, fit.lp:=predict(cox.fit, type="link", newx=Xnew)]
                } else {
                    tmp[, fit.lp:=predict(cox.fit, type="lp", newdata=tmp)]
                }
                tmp <- tmp[rev(order(get(time.var)))]
                tmp[, term2:=cumsum(risk*exp(fit.lp))]
                tmp[term2==0, term2:=1]
                return(sum(tmp[risk==0, get(delta.var)*(fit.lp - log(term2))]))
                #return(sum(tmp[risk==0, get(delta.var)*(fit.lp - log(term2))]))
            }
        }
   
        outlist2 <- list()
        
        for (vv in 1:V) {

            test.set <- cv.split[,vv]#sample(1:n, floor(n/10))
            train.set <- dt[, id][!dt[, id] %in% test.set]

            outlist2[[vv]] <- lapply(outcome.models, function(outcome.model) {
                if (length(names(outcome.model))==0) names(outcome.model) <- "no.name"
                tryCatch(
                    if (length(outcome.model)==1 | names(outcome.model)[1]=="coxnet") {
                        if (any(names(outcome.model)=="coxnet")) {
                            X <- model.matrix(outcome.model[[1]], data=dt[id%in%train.set])
                            y <- dt[id%in%train.set, Surv(get(time.var), get(delta.var)==1)]
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
                        dt2[period==1 & !time.indicator, (delta.var):=0]

                        mod1 <- as.character(outcome.model[[1]])
                        # mod2 <- paste0(gsub(substr(mod1[2], which(strsplit(mod1[2], "")[[1]]=="(")+1,
                        #                            which(strsplit(mod1[2], "")[[1]]==",")-1), "tstart, tstop", mod1[2]),
                        #                "~", 
                        #                gsub(paste0("\\+", A.name), "", gsub(" ", "", paste0("I((period==1)&(", A.name, "==1))",
                        #                                                                     " + I((period==2)&(", A.name, "==1))", " + ",
                        #                                                                     mod1[3]))))
                        mod2 <- paste0(gsub(substr(mod1[2], which(strsplit(mod1[2], "")[[1]]=="(")+1,
                                                   which(strsplit(mod1[2], "")[[1]]==",")-1), paste0("tstart", ", tstop"), mod1[2]),
                                       "~", 
                                       gsub(paste0("\\+", A.name, "\\+"), "", gsub(paste0("\\+", A.name, " "), "", gsub(" ", "", paste0("I((period","==1)&(", A.name, "==1))",
                                                                                                                                        " + I((period", "==2)&(", A.name, "==1))", " + ",
                                                                                                                                        paste0("+",mod1[3]))))))

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

    if (verbose) print(cve2)
    
    if (only.cox.sl) {
        print(cbind(names(outcome.models), cve2))
        return(c(#method1=(names(outcome.models)[cve1==min(cve1[abs(cve1)<Inf])]),
            method2=(names(outcome.models)[cve2==min(cve2[abs(cve2)<Inf])])))
    }

    return(list(pick=names(outcome.models)[cve2==min(cve2[abs(cve2)<Inf])][1],
                cve=min(cve2[abs(cve2)<Inf])))
}






