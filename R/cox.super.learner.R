cox.sl <- function(dt, tau=1.2, V=10, A.name="A",
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

    #-- try make one split (V=1);

    dt[, L1.squared:=L1^2]
    dt[, L2.squared:=L2^2]
    dt[, L3.squared:=L3^2]
    
    n <- nrow(dt)

    cv.split <- matrix(sample(1:n, size=n), ncol=V)

    outlist <- list()
    #outlist2 <- list()
    #outcome.model <-  Surv(time, delta==1)~A+L1.squared
    #outcome.model <- Surv(time, delta==1)~A*L1.squared+A+L1.squared+L1*L2+L3
    #outcome.model <-  Surv(time, delta==1)~A+L1+L2+L3
    
    for (vv in 1:V) {
    
        test.set <- cv.split[,vv]#sample(1:n, floor(n/10))
        train.set <- dt[, id][!dt[, id] %in% test.set]

        outlist[[vv]] <- lapply(outcome.models, function(outcome.model) {
            tryCatch(
                if (length(outcome.model)==1) {
                    train.fit <- coxph(outcome.model[[1]], data=dt[id%in%train.set])
                    return(partial.loss.fun(train.fit, 1:n)-partial.loss.fun(train.fit, train.set))
                } else {

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
                    return(partial.loss.fun(train.fit, 1:n, t0=outcome.model[[2]])-
                           partial.loss.fun(train.fit, train.set, t0= outcome.model[[2]]))
                }
            , error=function(e) Inf)
        })

        #outlist2[[vv]] <- partial.loss.fun(train.fit, test.set)    
    }

    #sum(-unlist(outlist))
    #sum(-unlist(outlist2))

    cve <- unlist(lapply(1:length(outcome.models), function(mm) {
        sum(-unlist(lapply(outlist, function(out) out[[mm]])))
    }))

    return(names(outcome.models)[cve==min(cve[abs(cve)<Inf])])
}


































if (FALSE) {

    
    #-- get nelson-aalen from full data: 
    
    fit.baseline <- coxph(Surv(time, delta == 1)~1,
                          data=dt)

    bhaz.cox <- rbind(data.table(time=0, hazard=0),
                      merge(data.table(time=unique.times),
                            setDT(basehaz(fit.baseline, centered=centered)),
                            by="time", all.x=TRUE))
    bhaz.cox[, "dhaz":=c(0, diff(hazard))]
    bhaz.cox[is.na(dhaz), "dhaz":=0]
    setnames(bhaz.cox, "hazard", "chaz")
    bhaz.cox[, "chaz":=cumsum(dhaz)]

    #smooth.dhaz <- smooth.spline(bhaz.cox$time, bhaz.cox$dhaz)

    mat <- merge(mat, bhaz.cox, by="time")


    
    #--- cox model 1: time-varying hr

    time.model.fun <- function(dt.train, mat.test, outcome.model, t0=NULL,
                               mod.period1="", mod.period2="", name="1") {

        if (length(t0)>0) { #-- if there is a change-point:

            dt.train[, time.indicator:=(time<=t0)]

            dt2 <- rbind(dt.train, dt.train)[order(id)]
            dt2[, period:=1:.N, by="id"]

            dt2[period==1, `:=`(tstart=0, tstop=(time<=t0)*time+(time>t0)*t0)]
            dt2[period==2, `:=`(tstart=t0, tstop=time)]
            dt2[period==1 & !time.indicator, delta:=0]

            mod1 <- as.character(outcome.model)
            mod2 <- paste0(gsub(substr(mod1[2], which(strsplit(mod1[2], "")[[1]]=="(")+1,
                                       which(strsplit(mod1[2], "")[[1]]==",")-1), "tstart, tstop", mod1[2]),
                           "~", 
                           gsub("\\+A", "", gsub(" ", "", paste0("I((period==1)&(", A.name, "==1))", mod.period1,
                                                                 " + I((period==2)&(", A.name, "==1))", mod.period2, " + ",
                                                                 mod1[3]))))

            fit.cox <- coxph(formula(mod2), data=dt2[!time.indicator | period==1])
        } else { #-- if there is no change-point:
            fit.cox <- coxph(as.formula(deparse(outcome.model)), #outcome.model,
                             data=dt.train)
        }

        print(fit.cox)

        if (FALSE) {


            bhaz.cox <- rbind(data.table(time=0, hazard=0),
                              merge(data.table(time=unique.times),
                                    setDT(basehaz(fit.cox, centered=centered)),
                                    by="time", all.x=TRUE))
            bhaz.cox[, "dhaz":=c(0, diff(hazard))]
            bhaz.cox[is.na(dhaz), "dhaz":=0]
            setnames(bhaz.cox, "hazard", "chaz")
            bhaz.cox[, "chaz":=cumsum(dhaz)]

            if (length(t0)>0) {
                bhaz.cox <- rbind(bhaz.cox,
                                  data.table(time=t0, dhaz=0, chaz=0),
                                  data.table(time=tau, dhaz=0, chaz=0))[order(time)]
                bhaz.cox[, period:=(time<=t0)*1+(time>t0)*2]
                bhaz.cox[, chaz:=cumsum(dhaz), by="period"]
            }
            mat.test <- merge(mat.test, bhaz.cox, by="time")

            setnames(mat.test, "dhaz", paste0("dhaz", name))
            setnames(mat.test, "chaz", paste0("chaz", name))
            
            mat.test[is.na("dhaz"), dhaz:=0]
        
            if (length(t0)>0) {
                mat.test[, period:=(time<=t0)*1+(time>t0)*2]
                mat.test[, chaz:=cumsum(dhaz), by=c("id", "period")]
            } else {
                mat.test[, chaz:=cumsum(dhaz), by=c("id")]
            }

        }

        if (length(t0)>0) {
            #bhaz.cox <- rbind(bhaz.cox,
            #                  data.table(time=t0, dhaz=0, chaz=0),
            #                  data.table(time=tau, dhaz=0, chaz=0))[order(time)]
            mat.test[, period:=(time<=t0)*1+(time>t0)*2]
            mat.test[, chaz:=cumsum(dhaz), by=c("id", "period")]
        }

       # dt2[, (paste0("fit", name, "2")):=predict(fit.cox, newdata=dt2, type="risk")]

        mat.test[, (paste0("fit", name)):=predict(fit.cox, newdata=mat.test, type="risk")]
        
       # merge(mat.test, dt2[, c("id", "period", "fit12"), with=FALSE], by=c("id", "period"))
        


        #mat.test <- merge(mat.test, dt2[, c("id", "period", paste0("fit", name)), with=FALSE],
        #             by=c("id", "period"))

        #dt2[, (paste0("fit", name)):=predict(fit.cox, newdata=dt2, type="risk")]

        return(mat.test)
    }

    
    compute.loss.fun <- function(mat, fit, dfit) {
        return(mat[, - sum( (time<=tau) * (time<=time.obs) *
                            ( log(get(fit)) *
                             (delta==1 & time==time.obs) -
                              get(fit)*get(dfit) )), by="id"])
    }


    #-- try make one split (V=1);

    n <- nrow(dt)

    cv.split <- matrix(sample(1:n, size=n), ncol=V)

    outlist <- list()

    for (vv in 1:V) { 

        test.set <- cv.split[,vv]#sample(1:n, floor(n/10))
        train.set <- dt[, id][!dt[, id] %in% test.set]

        if (FALSE) {
            cox2.test <- time.model.fun(dt, mat, outcome.model, t0=0.9)
            cox6.test <- time.model.fun(dt, mat,
                                        Surv(time, delta==1)~A+L1+L2+L3)

            mean(compute.loss.fun(cox2.test, "fit1", "dhaz")[, "V1"][[1]])
            mean(compute.loss.fun(cox6.test, "fit1", "dhaz")[, "V1"][[1]])
        }
        
        cox1 <- time.model.fun(dt[id %in% train.set], mat[id %in% test.set], outcome.model, t0=0.3)
        cox2 <- time.model.fun(dt[id %in% train.set], mat[id %in% test.set], outcome.model, t0=0.9)
        cox3 <- time.model.fun(dt[id %in% train.set], mat[id %in% test.set], outcome.model, t0=1.1)
        cox4 <- time.model.fun(dt[id %in% train.set], mat[id %in% test.set], outcome.model, t0=1.5)
        cox5 <- time.model.fun(dt[id %in% train.set], mat[id %in% test.set], outcome.model)
        cox6 <- time.model.fun(dt[id %in% train.set], mat[id %in% test.set],
                               Surv(time, delta==1)~A+L1+L2+L3)
        cox7 <- time.model.fun(dt[id %in% train.set], mat[id %in% test.set],
                               Surv(time, delta==1)~A*L1.squared+A+L1.squared+L1*L2+L3)
        cox8 <- time.model.fun(dt[id %in% train.set], mat[id %in% test.set],
                               Surv(time, delta==1)~A)
        cox9 <- time.model.fun(dt[id %in% train.set], mat[id %in% test.set],
                               Surv(time, delta==1)~A+L3+L2, t0=1.1)
        cox10 <- time.model.fun(dt[id %in% train.set], mat[id %in% test.set],
                                Surv(time, delta==1)~A+L1+L3+L2, t0=1.1)
        cox11 <- time.model.fun(dt[id %in% train.set], mat[id %in% test.set],
                                Surv(time, delta==1)~A+L1.squared)
        cox12 <- time.model.fun(dt[id %in% train.set], mat[id %in% test.set],
                                Surv(time, delta==1)~A+L1)
        
        outlist[[vv]] <- c(cox1=mean(compute.loss.fun(cox1, "fit1", "dhaz")[, "V1"][[1]]),
                           cox2=mean(compute.loss.fun(cox2, "fit1", "dhaz")[, "V1"][[1]]),
                           cox3=mean(compute.loss.fun(cox3, "fit1", "dhaz")[, "V1"][[1]]),
                           cox4=mean(compute.loss.fun(cox4, "fit1", "dhaz")[, "V1"][[1]]),
                           cox5=mean(compute.loss.fun(cox5, "fit1", "dhaz")[, "V1"][[1]]),
                           cox6=mean(compute.loss.fun(cox6, "fit1", "dhaz")[, "V1"][[1]]),
                           cox7=mean(compute.loss.fun(cox7, "fit1", "dhaz")[, "V1"][[1]]),
                           cox8=mean(compute.loss.fun(cox8, "fit1", "dhaz")[, "V1"][[1]]),
                           cox9=mean(compute.loss.fun(cox9, "fit1", "dhaz")[, "V1"][[1]]),
                           cox10=mean(compute.loss.fun(cox10, "fit1", "dhaz")[, "V1"][[1]]),
                           cox11=mean(compute.loss.fun(cox11, "fit1", "dhaz")[, "V1"][[1]]),
                           cox12=mean(compute.loss.fun(cox12, "fit1", "dhaz")[, "V1"][[1]]))
    }

    return(c(cox1=mean(unlist(lapply(outlist, function(x) x["cox1"]))),
             cox2=mean(unlist(lapply(outlist, function(x) x["cox2"]))),
             cox3=mean(unlist(lapply(outlist, function(x) x["cox3"]))),
             cox4=mean(unlist(lapply(outlist, function(x) x["cox4"]))),
             cox5=mean(unlist(lapply(outlist, function(x) x["cox5"]))),
             cox6=mean(unlist(lapply(outlist, function(x) x["cox6"]))),
             cox7=mean(unlist(lapply(outlist, function(x) x["cox7"]))),
             cox8=mean(unlist(lapply(outlist, function(x) x["cox8"]))),
             cox9=mean(unlist(lapply(outlist, function(x) x["cox9"]))),
             cox10=mean(unlist(lapply(outlist, function(x) x["cox10"]))),
             cox11=mean(unlist(lapply(outlist, function(x) x["cox11"]))),
             cox12=mean(unlist(lapply(outlist, function(x) x["cox12"])))))  


}


