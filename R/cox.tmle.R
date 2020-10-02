cox.tmle <- function(dt,
                     outcome.model=Surv(time, delta==1)~A*L1.squared+A+L1.squared+L1*L2+L3, # RIGHT NOW MUST BE CALLED "time" AND "delta"
                     cr.model=Surv(time, delta==2)~A+L1, competing.risk=NULL, 
                     change.point=NULL, mod.period1="*L3", mod.period2="*L3",
                     tau=c(1.2),
                     km.cens=FALSE, SL.method=2,
                     cens.model=Surv(time, delta==0)~L1+L2+L3+A*L1,
                     treat.model=A~L1+L2+L3,
                     treat.effect=c("1", "0", "both", "stochastic"),
                     pi.star.fun=function(A, L) 0.2, 
                     output.km=FALSE, output.hr=FALSE, centered=TRUE,
                     poisson.initial=FALSE, lambda.cv=NULL,
                     lambda.cvs=seq(0, 0.008, length=51)[-1], 
                     cut.covars=5, cut.time=10, cut.time.A=4,
                     covars=c("L1", "L2", "L3"),
                     one.step=FALSE, deps.size=0.001, no.small.steps=100, # FIXME: test again!
                     penalize.time=TRUE, adjust.penalization=TRUE,
                     poisson.cens=FALSE, cut.L1.A=5, cut.L.interaction=5, 
                     maxIter=5, verbose=TRUE, browse=FALSE, browse2=FALSE,
                     browse3=FALSE, browse4=FALSE,
                     only.cox=FALSE, only.cox.sl=FALSE,
                     SL=FALSE, SL.cens=FALSE, SL.poisson=FALSE,
                     SL.outcome.models=list(mod1=c(Surv(time, delta==1)~A+L1+L2+L3, t0=0.9),
                                            mod1a=c(Surv(time, delta==1)~A+L1+L2+L3, t0=0.2765),
                                            mod2=c(Surv(time, delta==1)~A*L1.squared+L1*L2+L3, t0=NULL),
                                            mod3=c(Surv(time, delta==1)~A+L1.squared, t0=NULL),
                                            #mod4=c(Surv(time, delta==1)~A+L1+L2+L3, t0=1.2),
                                            mod5=c(Surv(time, delta==1)~A+L1+L2+L3, t0=0.3),
                                            mod6=c(Surv(time, delta==1)~A+L1+L2+L3, t0=NULL),
                                            mod7=c(Surv(time, delta==1)~A+L3.squared, t0=NULL),
                                            mod8=c(Surv(time, delta==1)~A+L1+L2, t0=0.9),
                                            mod8a=c(Surv(time, delta==1)~A+L1+L2, t0=0.2765),
                                            mod9=c(Surv(time, delta==1)~A, t0=NULL),
                                            mod10=c(Surv(time, delta==1)~L1+L2+L3+A*L1, t0=NULL),
                                            mod11=c(Surv(time, delta==1)~1, t0=NULL),
                                            mod12=c(Surv(time, delta==1)~A+L1+L2+L3, t0=0.192018),
                                            #mod13=c(Surv(time, delta==1)~A+L1+L2+L3+I((period==2)&(L3>=0.5)), t0=0.192018),
                                            #mod13=c(Surv(time, delta==1)~A+L1+L2+L3+I((period==2)&(L3>=0.5)), t0=0.192018*1.2),
                                            #mod14=c(Surv(time, delta==1)~A+L1+L2+L3+I((period==2)&(L3>=0.5)), t0=0.9),
                                            mod15=c(Surv(time, delta==1)~A*L3+L1+L2, t0=0.2765),
                                            mod16=c(Surv(time, delta==1)~A*L3+L1+L2, t0=0.9),
                                            mod17=c(Surv(time, delta==1)~A*L3+L1+L2, t0=NULL)
                                            )) {
    
    #-- 0 -- some initializations:

    if (length(dt[, unique(delta)])>2 & length(competing.risk)==0) competing.risk <- TRUE
    if (length(competing.risk)==0) competing.risk <- FALSE 

    if (length(grep(".squared", as.character(outcome.model)[3]))>0) {
        names.squared <- unique(gsub(".squared", "",
                                     grep(".squared", unlist(strsplit(gsub("\\+", " ",
                                                                           as.character(outcome.model)[3]), " ")),
                                          value=TRUE)))
        for (col in names.squared)
            dt[, (paste0(col, ".squared")):=get(col)^2]
    }

    if (length(grep(".root", as.character(outcome.model)[3]))>0) {
        names.root <- unique(gsub(".root", "",
                                     grep(".root", unlist(strsplit(gsub("\\+", " ",
                                                                           as.character(outcome.model)[3]), " ")),
                                          value=TRUE)))
        for (col in names.root)
            dt[, (paste0(col, ".root")):=sqrt(get(col))]
    }

    
    if (length(grep(".sin", as.character(outcome.model)[3]))>0) {
        names.sin <- unique(gsub(".sin", "",
                                 grep(".sin", unlist(strsplit(gsub("\\+", " ",
                                                                   as.character(outcome.model)[3]), " ")),
                                      value=TRUE)))
        for (col in names.sin)
            dt[, (paste0(col, ".sin")):=sin(get(col)*8)]
    }

    #-- get number of subjects:
    n <- length(dt[, unique(id)])
    
    #-- get treatment colname:
    A.name <- as.character(treat.model)[2]

    #-- get unique times in dataset
    unique.times <- sort(unique(dt[, time]))

    #-- which parameters are we interested in?
    if (treat.effect[1]=="1") a <- 1 else if (treat.effect[1]=="0") a <- 0 else a <- c(1, 0)

    #-- 1 -- estimate censoring distribution:

    if (km.cens) {
        km.mod <- paste0(gsub("Surv", "Hist", as.character(cens.model)[2]), "~1")
        km.fit <- summary(fit.km <- prodlim(formula(km.mod), data=dt),
                          times=unique.times, asMatrix=TRUE)$table
        cens.km <- data.table(km.fit[, c("time", "surv")])
    }

    if (SL.cens) { #-- cox-SL

        print("use SL for censoring")

        SL.cens.outcome.models <- lapply(SL.outcome.models, function(mod) {
            mod1 <- as.character(mod[[1]])
            return(c(as.formula(paste0(gsub("1", "0", mod1[2]), mod1[1], mod1[3])), t0=mod[2][[1]]))
        })

        dtcens <- copy(dt)
        dtcens[, delta:=1-delta]

        SL.pick <- cox.sl(dtcens, tau=tau, A.name=A.name,
                          outcome.models=SL.cens.outcome.models)

        SL.model <- SL.cens.outcome.models[[SL.pick]]

        cens.model <- SL.model[[1]]

        print(paste0("model picked for censoring: ", cens.model))

        if (length(grep(".squared", as.character(cens.model)[3]))>0) {
            names.squared <- unique(gsub(".squared", "",
                                         grep(".squared", unlist(strsplit(gsub("\\+", " ",
                                                                               as.character(cens.model)[3]), " ")),
                                              value=TRUE)))
            for (col in names.squared)
                dt[, (paste0(col, ".squared")):=get(col)^2]
        }


    }

    if (length(cens.model)>0) {
        cens.cox <- coxph(as.formula(deparse(cens.model)), data=dt)
    }

    #-- 2 -- estimate treatment propensity: 

    prob.A <- predict(glm(as.formula(deparse(treat.model)), data=dt), type="response")

    #-- 3 -- estimate outcome distribution:

    if (SL) { #-- cox-SL

        print("use SL")

        if (only.cox.sl) SL.method <- c(1,2)
        SL.pick <- cox.sl(dt, tau=tau, A.name=A.name, only.cox.sl=only.cox.sl,
                          method=SL.method, covars=covars,
                          outcome.models=SL.outcome.models)
        if (only.cox.sl) return(SL.pick)
        
        SL.model <- SL.outcome.models[[SL.pick]]

        if (length(SL.model)>1) {
            print(paste0("model picked: ", outcome.model <- SL.model[[1]]))
            print(paste0("model picked: ", change.point <- SL.model[[2]]))
        } else {
            print(paste0("model picked: ", outcome.model <- SL.model[[1]]))
        }

    }

    #-- Apply either specified model or the one picked by SL
    
    if (length(change.point)>0) { #-- if there is a change-point:

        print("estimate time-varying hazard")
        
        dt[, time.indicator:=(time<=t0)]

        dt2 <- rbind(dt, dt)[order(id)]
        dt2[, period:=1:.N, by="id"]

        #dt2 <- dt2[!time.indicator | period==1]
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
                         data=dt)
    }

    if (competing.risk) {
        fit.cr.cox <- coxph(as.formula(deparse(cr.model)), #outcome.model,
                            data=dt)
    }
    
    if (verbose) print(fit.cox)

    if (only.cox) return()
    
    #-- 4 -- use kaplan-meier for estimation:

    if (output.km) {
        km.mod <- paste0(gsub("Surv", "Hist", as.character(outcome.model)[2]),
                         "~", A.name)
        km.fit <- summary(fit.km <- prodlim(formula(km.mod), data=dt),
                          times=tau, asMatrix=TRUE)$table
        if (length(a)==1) {
            km.est <- 1-as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", a),]["surv"])
            km.se <- as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", a),]["se.surv"])
        } else {
            km.est <- 1-as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "1"),]["surv"])-
                (1-as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "0"),]["surv"]))
            km.se <- sqrt((as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "1"),]["se.surv"]))^2+
                          (as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", "0"),]["se.surv"]))^2)
        }
    }

    #-- 5 -- estimate crude HR:

    if (output.hr) {
        hr.mod <- paste0(as.character(outcome.model)[2],
                         "~", A.name)
        hr <- coxph(formula(hr.mod), data=dt)
        if (verbose) print(hr)
        hr.pval <- summary(hr)$coefficients[A.name, 5]
        hr <- coef(hr)[A.name]
    }

    #-- 6 -- get baseline hazard:

    bhaz.cox <- rbind(data.table(time=0, hazard=0),
                      merge(data.table(time=unique.times),
                            setDT(basehaz(fit.cox, centered=centered)),
                            by="time", all.x=TRUE))
    bhaz.cox[, dhaz:=c(0, diff(hazard))]
    setnames(bhaz.cox, "hazard", "chaz")

    if (competing.risk) {
        bhaz.cox <- merge(bhaz.cox, rbind(data.table(time=0, hazard=0),
                                          setDT(basehaz(fit.cr.cox, centered=centered))),
                          by="time", all.x=TRUE)
        setnames(bhaz.cox, "hazard", "cr.chaz")
        bhaz.cox[, cr.dhaz:=c(0, diff(cr.chaz))]
    }
    
    #-- 6b -- add censoring baseline hazard:

    if (length(cens.model)>0) {
        bhaz.cox <- merge(bhaz.cox, rbind(data.table(time=0, hazard=0),
                                          setDT(basehaz(cens.cox, centered=centered))),
                          by="time", all.x=TRUE)
        setnames(bhaz.cox, "hazard", "cens.chaz")
        bhaz.cox[, cens.dhaz:=c(0, diff(cens.chaz))]

        if (length(change.point)>0) {
            bhaz.cox <- rbind(bhaz.cox,
                              data.table(time=change.point, dhaz=0, cens.dhaz=0, chaz=0, cens.chaz=0),
                              data.table(time=tau, dhaz=0, cens.dhaz=0, chaz=0, cens.chaz=0))[order(time)]
            bhaz.cox[, period:=(time<=change.point)*1+(time>change.point)*2]
            bhaz.cox[, chaz:=cumsum(dhaz), by="period"]
            bhaz.cox[, cens.chaz:=cumsum(cens.dhaz)]
        }

        #-- 6c -- get censoring survival one time-point back: 

        bhaz.cox[, cens.chaz.1:=c(0, cens.chaz[-.N])]
    }
    
    #-- 7 -- dublicate bhaz.cox; for each treatment option:

    if (poisson.initial) {
        mat.cox <- do.call("rbind", lapply(a, function(aa) data.table(A=aa, bhaz.cox)))
    } else {
        mat.cox <- do.call("rbind", lapply(a, function(aa) data.table(A=aa, bhaz.cox)[time<=tau]))
    }

    #-- 8 -- add subject-specific information:

    dt.a <- do.call("rbind", lapply(a, function(aa) {
        dt.tmp <- copy(dt)
        dt.tmp[, A:=aa]
    }))

    if (length(change.point)>0) {
        dt2.a <- do.call("rbind", lapply(a, function(aa) {
            dt.tmp <- copy(dt2)
            dt.tmp[, A:=aa]
        }))
        dt2.a[, fit.cox:=predict(fit.cox, newdata=dt2.a, type="risk")]
        if (competing.risk) 
            dt2.a[, fit.cr.cox:=predict(fit.cr.cox, newdata=dt2.a, type="risk")]
    } else {
        dt.a[, fit.cox:=predict(fit.cox, newdata=dt.a, type="risk")]
        if (competing.risk)
            dt.a[, fit.cr.cox:=predict(fit.cr.cox, newdata=dt.a, type="risk")]
    }

    if (length(cens.model)>0) {
        dt.a[, fit.cens.a:=predict(cens.cox, newdata=dt.a, type="risk")]
    }
    
    if (km.cens) {
        setnames(cens.km, "surv", "surv.cens")
        dt.a <- merge(dt.a, cens.km)
    }

    mat <- do.call("rbind", lapply(1:n, function(i) {
        tmp <- cbind(id=i, mat.cox)
        tmp[, time.obs:=dt[id==i, time]]
        tmp[, delta.obs:=dt[id==i, delta]]
        tmp[, A.obs:=dt[id==i, get(A.name)]]
        if (treat.effect[1]=="stochastic") {
            tmp <- merge(tmp, dt[, c("id", covars), with=FALSE], by="id")
            tmp[, pi.star:=pi.star.fun(get(A.name), cbind(L1,L2,L3))]
            tmp[, pi.star:=pi.star*get(A.name)+(1-pi.star)*(1-get(A.name))]
            tmp[, Ht:=-(pi.star / # treatment and censoring weights
                      ((prob.A[i]^get(A.name) * (1-prob.A[i])^(1-get(A.name)))))]
        } else {
            tmp[, Ht:=-((get(A.name)==1) - (get(A.name)==0)) / # treatment and censoring weights
                      ((prob.A[i]^get(A.name) * (1-prob.A[i])^(1-get(A.name))))]
        }
        if (length(cens.model)>0) {
            if (km.cens) {
                tmp <- merge(tmp, dt.a[, c("surv.cens", "time", A.name), with=FALSE],
                             by=c("time", A.name))
                tmp[, surv.km.C1:=c(1, surv.cens[-.N])]
            }  
            tmp <- merge(tmp, dt.a[id==i, c(A.name, "fit.cens.a"), with=FALSE], by=A.name)
            tmp[, surv.C1:=exp(-fit.cens.a*cens.chaz.1)]
            if (km.cens) {
                tmp[, Ht:=Ht/surv.km.C1]
            } else {
                tmp[, Ht:=Ht/exp(-fit.cens.a*cens.chaz.1)]
            }
        } 
        if (length(change.point)>0) {
            tmp[, period:=(time<=change.point)*1+(time>change.point)*2]
            if (competing.risk) {
                tmp <- merge(tmp, dt2.a[id==i, c("period", A.name, "fit.cox", "fit.cr.cox"), with=FALSE], by=c("period", A.name))
            } else {
                tmp <- merge(tmp, dt2.a[id==i, c("period", A.name, "fit.cox"), with=FALSE], by=c("period", A.name))
            }
        } else {
            if (competing.risk) {
                tmp <- merge(tmp, dt.a[id==i, c(A.name, "fit.cox", "fit.cr.cox"), with=FALSE], by=A.name)
            } else {
                tmp <- merge(tmp, dt.a[id==i, c(A.name, "fit.cox"), with=FALSE], by=A.name)
            }
        }
    }))
    
    #-- 9 -- compute clever covariates:

    if (competing.risk) {
        mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)-cumsum(cr.dhaz*fit.cr.cox)), by=c("id", "A")]
        mat[, surv.t1:=c(0, surv.t[-.N]), by=c("id", "A")]
        mat[, surv.tau:=surv.t[time==max(time[time<=tau])], by=c("id", "A")]
        mat[, F1.t:=cumsum(surv.t1*dhaz*fit.cox), by=c("id", "A")]
        mat[, F1.tau:=F1.t[time==max(time[time<=tau])], by=c("id", "A")]
    } else{
        mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
        mat[, surv.tau:=surv.t[time==max(time[time<=tau])], by=c("id", "A")]
    }
    
    if (browse) browser()

    #-- 10 -- poisson used for initial:

    if (poisson.initial) {
        mat <- poisson.hal(mat=mat, dt=dt, outcome.model=outcome.model,
                           change.point=change.point, X=NULL, delta.outcome=1, verbose=verbose,
                           cut.covars=cut.covars, cut.time=cut.time, browse=FALSE,
                           cut.time.A=cut.time.A,
                           cut.L1.A=cut.L1.A, cut.L.interaction=cut.L.interaction,
                           covars=covars,
                           poisson.cens=poisson.cens,
                           SL.poisson=SL.poisson, 
                           lambda.cv=lambda.cv, lambda.cvs=lambda.cvs,
                           penalize.time=penalize.time, adjust.penalization=adjust.penalization)

        #-- 10a -- poisson used for censoring distribution:

        if (poisson.cens) {
            mat <- poisson.hal(mat=mat[["mat"]], dt=dt, outcome.model=outcome.model,
                               change.point=change.point, X=mat[["X"]], delta.outcome=0, verbose=verbose,
                               cut.covars=cut.covars, cut.time=cut.time,
                               cut.time.A=cut.time.A, browse=FALSE,
                               cut.L1.A=cut.L1.A, cut.L.interaction=cut.L.interaction,
                               covars=covars,
                               poisson.cens=poisson.cens,
                               SL.poisson=SL.poisson, 
                               lambda.cv=lambda.cv, lambda.cvs=lambda.cvs,
                               penalize.time=penalize.time, adjust.penalization=adjust.penalization)
        }
    }

    if (competing.risk) {
        mat[surv.t>0, Ht.lambda:=1 - (F1.tau - F1.t) / surv.t]
        mat[surv.t==0, Ht.lambda:=1]
        mat[surv.t>0, Ht.lambda2:=-(F1.tau - F1.t) / surv.t]
        mat[surv.t==0, Ht.lambda2:=1]
    } else {
        mat[surv.t>0, Ht.lambda:=surv.tau / surv.t]
        mat[surv.t==0, Ht.lambda:=1]
    }

    #-- 11 -- initial fit:
    
    #init.fit <- mean(mat[get(A.name)==a, 1-surv.tau[1], by="id"][,2][[1]])

    if (treat.effect[1]=="stochastic") {
        if (competing.risk) {
            init.fit <- mean(rowSums(sapply(a, function(aa)
            (mat[get(A.name)==aa, pi.star*F1.tau[1], by="id"][,2][[1]]))))
        } else {
            init.fit <- mean(rowSums(sapply(a, function(aa)
            (mat[get(A.name)==aa, pi.star*(1-surv.tau[1]), by="id"][,2][[1]]))))
        }
    } else {
        if (competing.risk) {
            init.fit <- mean(rowSums(sapply(a, function(aa)
            (2*(aa==a[1])-1)*(mat[get(A.name)==aa, F1.tau[1], by="id"][,2][[1]]))))
        } else {
            init.fit <- mean(rowSums(sapply(a, function(aa)
            (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-surv.tau[1], by="id"][,2][[1]]))))
        }
    }

    if (treat.effect[1]=="stochastic") {
        if (competing.risk) {
            eval.ic <- function(mat) {
                out <- mat[, sum( (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda *
                                  ( (delta.obs==1 & time==time.obs) -
                                    dhaz * fit.cox ) +
                                  (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda2 *
                                  ( (delta.obs==2 & time==time.obs) -
                                    cr.dhaz * fit.cr.cox )), by="id"]
                ic.squared <- (out[, 2][[1]] +
                               rowSums(sapply(a, function(aa)
                               (mat[get(A.name)==aa, pi.star[1]*F1.tau[1], by="id"][,2][[1]]))) -
                               (length(a)==2 & !(treat.effect[1]=="stochastic"))*init.fit-
                               (length(a)==1 | (treat.effect[1]=="stochastic"))*(1-init.fit))^2
                return(sqrt(mean(ic.squared)/n))
            }
        } else {
            eval.ic <- function(mat) {
                out <- mat[, sum( pi.star * (time<=time.obs) * Ht * Ht.lambda *
                                  ( (delta.obs==1 & time==time.obs) -
                                    dhaz * fit.cox )), by="id"]
                ic.squared <- (out[, 2][[1]] +
                               rowSums(sapply(a, function(aa)
                               (mat[get(A.name)==aa, pi.star[1]*(1-surv.tau[1]), by="id"][,2][[1]]))) -
                               (length(a)==2 & !(treat.effect[1]=="stochastic"))*init.fit-
                               (length(a)==1 | (treat.effect[1]=="stochastic"))*(1-init.fit))^2
                return(sqrt(mean(ic.squared)/n))
            }
        }
    } else {
        if (competing.risk) {
            eval.ic <- function(mat) {
                out <- mat[, sum( (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda *
                                  ( (delta.obs==1 & time==time.obs) -
                                    dhaz * fit.cox ) +
                                  (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda2 *
                                  ( (delta.obs==2 & time==time.obs) -
                                    cr.dhaz * fit.cr.cox )), by="id"]
                ic.squared <- (out[, 2][[1]] +
                               rowSums(sapply(a, function(aa)
                               (2*(aa==a[1])-1)*(mat[get(A.name)==aa, F1.tau[1], by="id"][,2][[1]]))) -
                               (length(a)==2)*init.fit-
                               (length(a)==1)*(1-init.fit))^2
                return(sqrt(mean(ic.squared)/n))
            }
        } else {
            eval.ic <- function(mat) {
                out <- mat[, sum( (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda *
                                  ( (delta.obs==1 & time==time.obs) -
                                    dhaz * fit.cox )), by="id"]
                ic.squared <- (out[, 2][[1]] +
                               rowSums(sapply(a, function(aa)
                               (2*(aa==a[1])-1)*(mat[get(A.name)==aa, surv.tau[1], by="id"][,2][[1]]))) -
                               (length(a)==2)*init.fit-
                               (length(a)==1)*(1-init.fit))^2
                return(sqrt(mean(ic.squared)/n))
            }
        }
    }

    init.ic <- eval.ic(mat)

    #-- 11a -- add to list to be outputted:
    tmle.list <- list(c(init=init.fit, sd.eic=init.ic))
    if (output.km) tmle.list[[1]] <- c(tmle.list[[1]], km.est=km.est, km.se=km.se)
    if (output.hr) tmle.list[[1]] <- c(tmle.list[[1]], hr=hr, hr.pval=hr.pval)

    #-- 12 -- tmle:

    if (one.step) { # simple one-step! 

        eval.equation <- function(mat, eps, j=1) {
            if (j==1) {
                out <- mat[get(A.name)==A.obs, sum( (time<=time.obs) * Ht * Ht.lambda *
                                                    ( (delta.obs==1 & time==time.obs) -
                                                      exp( eps * Ht * Ht.lambda ) * dhaz * fit.cox )),
                           by="id"]
            } else {
                out <- mat[get(A.name)==A.obs, sum( (time<=time.obs) * Ht * Ht.lambda2 *
                                                    ( (delta.obs==2 & time==time.obs) -
                                                      exp( eps * Ht * Ht.lambda2 ) * cr.dhaz * fit.cr.cox )),
                           by="id"]
            }
            return(mean(out[, 2][[1]])) 
        }

        
        #--- decide on direction of small eps increments:
        if (competing.risk) {
            if (abs(eval.equation(mat, -0.01)+eval.equation(mat, -0.01, j=2))<
                abs(eval.equation(mat, 0.01)+eval.equation(mat, 0.01, j=2))) {
                deps <- -deps.size
            } else {
                deps <- deps.size
            }
        } else if (abs(eval.equation(mat, -0.01))<abs(eval.equation(mat, 0.01))) {
            deps <- -deps.size
        } else {
            deps <- deps.size
        }

        #--- initial sign of eic equation:
        if (competing.risk) {
            sign.eic <- sign(eval.equation(mat, 0)+eval.equation(mat, 0, j=2))#-0.01))
        } else {
            sign.eic <- sign(eval.equation(mat, 0))#-0.01))
        }
        #sign.eic <- sign(eval.equation(mat, -0.01))
        
        #--- one-step: track eic equation with small step deps: 
        for (step in 1:no.small.steps) {

            if (competing.risk) {
                mat[, fit.cox:=fit.cox*exp(deps*Ht.lambda*Ht)]
                mat[, fit.cr.cox:=fit.cr.cox*exp(deps*Ht.lambda2*Ht)]
                mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)-cumsum(cr.dhaz*fit.cr.cox)), by=c("id", "A")]
                mat[fit.cox==Inf, surv.t:=0]
                mat[, surv.tau:=surv.t[.N], by=c("id", "A")]
                mat[, F1.t:=cumsum(surv.t1*dhaz*fit.cox), by=c("id", "A")]
                mat[, F1.tau:=F1.t[time==max(time[time<=tau])], by=c("id", "A")]
            } else {
                mat[, fit.cox:=fit.cox*exp(deps*Ht.lambda*Ht)]
                mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
                mat[, surv.tau:=surv.t[.N], by=c("id", "A")]
            }

            print(eic.value <- eval.equation(mat, deps))

            if  (sign.eic*eic.value<=0) {
                break
            }            
        }

        eval.equation(mat, 0)

        #-- 12a -- evaluate target parameter:
        if (treat.effect[1]=="stochastic") {
            if (competing.risk) {
                tmle.fit <- mean(rowSums(sapply(a, function(aa)
                (mat[get(A.name)==aa, pi.star*F1.tau[1], by="id"][,2][[1]]))))
            } else {
                tmle.fit <- mean(rowSums(sapply(a, function(aa)
                (mat[get(A.name)==aa, pi.star*(1-surv.tau[1]), by="id"][,2][[1]]))))
            }        
        } else {
            if (competing.risk) {
                tmle.fit <- mean(rowSums(sapply(a, function(aa)
                (2*(aa==a[1])-1)*(mat[get(A.name)==aa, F1.tau[1], by="id"][,2][[1]]))))
            } else {
                tmle.fit <- mean(rowSums(sapply(a, function(aa)
                (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-surv.tau[1], by="id"][,2][[1]]))))
            }
        }
        
        #-- 12b -- update clever covariate:
        if (competing.risk) {
            mat[surv.t>0, Ht.lambda:=1-(F1.tau - F1.t) / surv.t]
            mat[surv.t==0, Ht.lambda:=1]
            mat[surv.t>0, Ht.lambda2:=-(F1.tau - F1.t) / surv.t]
            mat[surv.t==0, Ht.lambda2:=1]
        } else {
            mat[surv.t>0, Ht.lambda:=surv.tau/surv.t]
            mat[surv.t==0, Ht.lambda:=1]
        }

        #-- 12c -- compute sd:
        if (FALSE) { # sep 24: was this not redundant??
            eval.ic <- function(mat) {
                out <- mat[, sum( (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda *
                                  ( (delta.obs==1 & time==time.obs) -
                                    dhaz * fit.cox )), by="id"]
                ic.squared <- (out[, 2][[1]] +
                               rowSums(sapply(a, function(aa)
                               (2*(aa==a[1])-1)*(mat[get(A.name)==aa, surv.tau[1], by="id"][,2][[1]]))) -
                               (length(a)==2)*tmle.fit-
                               (length(a)==1)*(1-tmle.fit))^2
                return(sqrt(mean(ic.squared)/n))
            }
        }

        tmle.list[[length(tmle.list)+1]] <- c(tmle.fit=tmle.fit, sd.eic=eval.ic(mat))
            
    }
    
    if (!one.step) { # iterative tmle 
    
        for (iter in 1:maxIter) {

            #-- 12a -- estimate eps:
            eval.equation <- function(mat, eps, j=1) {
                if (j==1) {
                    out <- mat[get(A.name)==A.obs, sum( (time<=time.obs) * Ht * Ht.lambda *
                                                        ( (delta.obs==1 & time==time.obs) -
                                                          exp( eps * Ht * Ht.lambda ) * dhaz * fit.cox )),
                               by="id"]
                } else {
                    out <- mat[get(A.name)==A.obs, sum( (time<=time.obs) * Ht * Ht.lambda2 *
                                                        ( (delta.obs==2 & time==time.obs) -
                                                          exp( eps * Ht * Ht.lambda2 ) * cr.dhaz * fit.cr.cox )),
                               by="id"]
                }
                return(mean(out[, 2][[1]])) 
            }

            print(paste0("iter=", iter, ", estimate eps: ",
                         round(eps.hat <- nleqslv(0.01, function(eps) eval.equation(mat, eps))$x, 4)))

            if (competing.risk) {
                print(paste0("cr: ", "estimate eps: ",
                             round(eps.hat2 <- nleqslv(0.01, function(eps) eval.equation(mat, eps, j=2))$x, 4)))
            }
            
            #-- 12b -- update hazard(s):
            
            if (competing.risk) {
                mat[, fit.cox:=fit.cox*exp(eps.hat*Ht.lambda*Ht)]
                mat[, fit.cr.cox:=fit.cr.cox*exp(eps.hat2*Ht.lambda2*Ht)]
                mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)-cumsum(cr.dhaz*fit.cr.cox)), by=c("id", "A")]
                mat[fit.cox==Inf, surv.t:=0]
                mat[, surv.tau:=surv.t[.N], by=c("id", "A")]
                mat[, F1.t:=cumsum(surv.t1*dhaz*fit.cox), by=c("id", "A")]
                mat[, F1.tau:=F1.t[time==max(time[time<=tau])], by=c("id", "A")]
            } else {
                mat[, fit.cox:=fit.cox*exp(eps.hat*Ht.lambda*Ht)]
                mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
                mat[fit.cox==Inf, surv.t:=0]
                mat[, surv.tau:=surv.t[.N], by=c("id", "A")]
            }

            #-- 12c -- evaluate target parameter:
            if (treat.effect[1]=="stochastic") {
                if (competing.risk) {
                    tmle.fit <- mean(rowSums(sapply(a, function(aa)
                    (mat[get(A.name)==aa, pi.star*F1.tau[1], by="id"][,2][[1]]))))
                } else {
                    tmle.fit <- mean(rowSums(sapply(a, function(aa)
                    (mat[get(A.name)==aa, pi.star*(1-surv.tau[1]), by="id"][,2][[1]]))))
                }        
            } else {
                if (competing.risk) {
                    tmle.fit <- mean(rowSums(sapply(a, function(aa)
                    (2*(aa==a[1])-1)*(mat[get(A.name)==aa, F1.tau[1], by="id"][,2][[1]]))))
                } else {
                    tmle.fit <- mean(rowSums(sapply(a, function(aa)
                    (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-surv.tau[1], by="id"][,2][[1]]))))
                }
            }
        
            #-- 12b -- update clever covariate:
            if (competing.risk) {
                mat[surv.t>0, Ht.lambda:=1-(F1.tau - F1.t) / surv.t]
                mat[surv.t==0, Ht.lambda:=1]
                mat[surv.t>0, Ht.lambda2:=-(F1.tau - F1.t) / surv.t]
                mat[surv.t==0, Ht.lambda2:=1]
            } else {
                mat[surv.t>0, Ht.lambda:=surv.tau/surv.t]
                mat[surv.t==0, Ht.lambda:=1]
            }

            #-- 12e -- compute sd:
            if (FALSE) { # sep 24: was this not redundant??
                eval.ic <- function(mat) {
                    out <- mat[, sum( (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda *
                                      ( (delta.obs==1 & time==time.obs) -
                                        dhaz * fit.cox )), by="id"]
                    ic.squared <- (out[, 2][[1]] +
                                   rowSums(sapply(a, function(aa)
                                   (2*(aa==a[1])-1)*(mat[get(A.name)==aa, surv.tau[1], by="id"][,2][[1]]))) -
                                   (length(a)==2)*tmle.fit-
                                   (length(a)==1)*(1-tmle.fit))^2
                    return(sqrt(mean(ic.squared)/n))
                }
            }

            #if (treat.effect[1]=="stochastic") {
            #    tmle.list[[iter+1]] <- c(tmle.fit=tmle.fit, sd.eic=eval.ic(mat))
            #} else {
                tmle.list[[iter+1]] <- c(tmle.fit=tmle.fit, sd.eic=eval.ic(mat))
            #}

            if (competing.risk) {
                eval.iter <- abs(eval.equation(mat, 0)+eval.equation(mat, 0, j=2))
            } else {
                eval.iter <- abs(eval.equation(mat, 0))
            }
            
            if (eval.iter<=eval.ic(mat)/(sqrt(n)*log(n))) {
                break
            }
        }
    }

    return(tmle.list)    
}












if (FALSE) {
    mat.big.melt[, Ht.2:=Ht.lambda*Ht]
                
    test <- dcast(mat.big.melt, id+tau~k, value.var="Ht.2")
    test2 <- dcast(mat.big.melt, id+tau~k, value.var="dhaz.fit")

    # unique(mat.big.melt[, c("id", "tau", "dhaz.fit"), with=FALSE])
            
    for (jj in names(test))
        test[is.na(get(jj)), (jj):=0]
    for (jj in names(test2))
        test2[is.na(get(jj)), (jj):=0]
                 
    exp(as.matrix(test[, -c("id", "tau")]) %*% as.matrix(delta, nrow=1))
    test3 <- merge(test[, c("id", "tau"), with=FALSE][, exp:=drop(exp(as.matrix(test[, -c("id", "tau")]) %*% as.matrix(delta, nrow=1)))],
                   test2)
    for (jj in names(test3[, -c("id", "tau", "exp")]))
        test3[, (jj):=get(jj)*exp]

    test4 <- melt(test3[, -"exp"], id.vars=c("id", "tau"), variable.name="k", value.name="dhaz.fit")
    test4[, k:=as.numeric(as.character(k))]

    mat.big.melt <- merge(mat.big.melt, test4[k<=tau], by=c("id", "tau", "k"))
    setnames(mat.big.melt, c("dhaz.fit.x", "dhaz.fit.y"), c("dhaz.fit.old", "dhaz.fit"))
}
            
