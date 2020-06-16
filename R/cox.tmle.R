cox.tmle <- function(dt, outcome.model=Surv(time, delta==1)~A*L1.squared+A+L1.squared+L1*L2+L3, # FIXME: RIGHT NOW MUST BE CALLED "time" AND "delta"
                     change.point=NULL, mod.period1="", mod.period2="*L3",
                     tau=c(1.2),
                     cens.model=Surv(time, delta==0)~L1+L2+L3+A*L1,
                     treat.model=A~L1+L2+L3,
                     treat.effect=c("1", "0", "both"), 
                     output.km=FALSE, output.hr=FALSE,
                     maxIter=5, verbose=TRUE, browse=FALSE) {

    #-- 0 -- some initializations:

    if (length(grep(".squared", as.character(outcome.model)[3]))>0) {
        names.squared <- unique(gsub(".squared", "",
                                     grep(".squared", unlist(strsplit(gsub("\\+", " ",
                                                                           as.character(outcome.model)[3]), " ")),
                                          value=TRUE)))
        for (col in names.squared)
            dt[, (paste0(col, ".squared")):=get(col)^2]
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

    cens.cox <- coxph(cens.model, data=dt)

    #-- 2 -- estimate treatment propensity: 

    prob.A <- predict(glm(treat.model, data=dt), type="response")

    #-- 3 -- estimate outcome distribution:    

    if (length(change.point)>0) { #-- if there is a change-point:

        print("estimate time-varying hazard")
        
        dt[, time.indicator:=(time<=t0)]

        dt2 <- rbind(dt, dt)[order(id)]
        dt2[, period:=1:.N, by="id"]

        dt2 <- dt2[!time.indicator | period==1]
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

        fit.cox <- coxph(formula(mod2), data=dt2)
    } else { #-- if there is no change-point:
        fit.cox <- coxph(outcome.model, data=dt)
    }
    if (verbose) print(fit.cox)

    #-- 4 -- use kaplan-meier for estimation:

    if (output.km) {
        km.mod <- paste0(gsub("Surv", "Hist", as.character(outcome.model)[2]),
                         "~", A.name)
        km.fit <- summary(fit.km <- prodlim(formula(km.mod), data=dt),
                          times=tau, asMatrix=TRUE)$table
        km.est <- 1-as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", a),]["surv"])
        km.se <- as.numeric(km.fit[km.fit[,1]==paste0(A.name, "=", a),]["se.surv"])
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

    #-- 6b -- add censoring baseline hazard:

    bhaz.cox <- merge(bhaz.cox, rbind(data.table(time=0, hazard=0),
                                      setDT(basehaz(cens.cox, centered=centered))),
                      by="time", all.x=TRUE)
    setnames(bhaz.cox, "hazard", "cens.chaz")
    bhaz.cox[, cens.dhaz:=c(0, diff(cens.chaz))]

    if (FALSE) {#(length(change.point)>0) { # CHECK: MAY NOT NEED THIS.
        bhaz.cox <- rbind(bhaz.cox,
                          data.table(time=change.point, dhaz=0, cens.dhaz=0, chaz=0, cens.chaz=0),
                          data.table(time=tau, dhaz=0, cens.dhaz=0, chaz=0, cens.chaz=0))[order(time)]
        bhaz.cox[, period:=(time<=change.point)*1+(time>change.point)*2]
        bhaz.cox[, chaz:=cumsum(dhaz), by="period"]
        bhaz.cox[, cens.chaz:=cumsum(cens.dhaz)]
    }

    #-- 6c -- get censoring survival one time-point back: 

    bhaz.cox[, cens.chaz.1:=c(0, cens.chaz[-.N])]

    #-- 7 -- dublicate bhaz.cox; for each treatment option:

    mat.cox <- do.call("rbind", lapply(a, function(aa) data.table(A=aa, bhaz.cox)[time<=tau]))

    #-- 8 -- add subject-specific information:
    if (browse) browser()
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
    } else {
        dt.a[, fit.cox:=predict(fit.cox, newdata=dt.a, type="risk")]
    }
    dt.a[, fit.cens.a:=predict(cens.cox, newdata=dt.a, type="risk")]

    mat <- do.call("rbind", lapply(1:n, function(i) {
        tmp <- cbind(id=i, mat.cox)
        tmp[, time.obs:=dt[id==i, time]]
        tmp[, delta.obs:=dt[id==i, delta]]
        tmp[, A.obs:=dt[id==i, get(A.name)]]
        tmp <- merge(tmp, dt.a[id==i, c(A.name, "fit.cens.a"), with=FALSE], by=A.name)
        tmp[, Ht:=-((get(A.name)==1) - (get(A.name)==0)) / # treatment and censoring weights
                  ((prob.A[i]^get(A.name) * (1-prob.A[i])^(1-get(A.name)))*
                   exp(-fit.cens.a*cens.chaz.1))]
        if (length(change.point)>0) {
            tmp[, period:=(time<=change.point)*1+(time>change.point)*2]
            tmp <- merge(tmp, dt2.a[id==i, c("period", A.name, "fit.cox"), with=FALSE], by=c("period", A.name))
        } else {
            tmp <- merge(tmp, dt.a[id==i, c(A.name, "fit.cox"), with=FALSE], by=A.name)
        }
    }))

    #-- 9 -- compute clever covariates:
    
    mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
    mat[, surv.tau:=surv.t[.N], by=c("id", "A")]
    mat[, Ht.lambda:=surv.tau / surv.t]

    #-- 10 -- initial fit:
    
    #init.fit <- mean(mat[get(A.name)==a, 1-surv.tau[1], by="id"][,2][[1]])
    init.fit <- mean(rowSums(sapply(a, function(aa)
    (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-surv.tau[1], by="id"][,2][[1]]))))

    eval.ic <- function(mat) {
        out <- mat[, sum( (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda *
                          ( (delta.obs==1 & time==time.obs) -
                            dhaz * fit.cox )), by="id"]
        
        #ic.squared <- (out[, 2][[1]] + (mat[get(A.name)==a, surv.tau[1], by="id"][,2][[1]]) - (1-init.fit))^2
        ic.squared <- (out[, 2][[1]] +
                       rowSums(sapply(a, function(aa)
                       (2*(aa==a[1])-1)*(mat[get(A.name)==aa, surv.tau[1], by="id"][,2][[1]]))) -
                       (length(a)==2)*init.fit-
                       (length(a)==1)*(1-init.fit))^2
        return(sqrt(mean(ic.squared)/n))
    }

    init.ic <- eval.ic(mat)

    #-- 10a -- add to list to be outputted:
    tmle.list <- list(c(init=init.fit, sd.eic=init.ic))
    if (output.km) tmle.list[[1]] <- c(tmle.list[[1]], km.est=km.est, km.se=km.se)
    if (output.hr) tmle.list[[1]] <- c(tmle.list[[1]], hr=hr, hr.pval=hr.pval)

    #-- 12 -- tmle:

    for (iter in 1:maxIter) {

        #-- 12a -- estimate eps:
        eval.equation <- function(mat, eps) {
            out <- mat[get(A.name)==A.obs, sum( (time<=time.obs) * Ht * Ht.lambda *
                                                ( (delta.obs==1 & time==time.obs) -
                                                  exp( eps * Ht * Ht.lambda ) * dhaz * fit.cox )),
                       by="id"]
            return(mean(out[, 2][[1]])) 
        }
        
        print(paste0("m=", m, ", iter=", iter, ", estimate eps: ",
                     round(eps.hat <- nleqslv(0.01, function(eps) eval.equation(mat, eps))$x, 4)))
    
        #-- 12b -- update hazard:
        mat[, fit.cox:=fit.cox*exp(eps.hat*Ht.lambda*Ht)]
        mat[, surv.t:=exp(-cumsum(dhaz*fit.cox)), by=c("id", "A")]
        mat[, surv.tau:=surv.t[.N], by=c("id", "A")]

        #-- 12c -- evaluate target parameter:x
        #tmle.fit <- mean(mat[get(A.name)==a, 1-surv.tau[1], by="id"][,2][[1]])
        tmle.fit <- mean(rowSums(sapply(a, function(aa)
        (2*(aa==a[1])-1)*(mat[get(A.name)==aa, 1-surv.tau[1], by="id"][,2][[1]]))))

        #-- 12d -- update clever covariate:
        mat[surv.t>0, Ht.lambda:=surv.tau/surv.t]
        mat[surv.t==0, Ht.lambda:=0]

        #-- 12e -- compute sd:
        eval.ic <- function(mat) {
            out <- mat[, sum( (get(A.name)==A.obs) * (time<=time.obs) * Ht * Ht.lambda *
                              ( (delta.obs==1 & time==time.obs) -
                                dhaz * fit.cox )), by="id"]
            ic.squared <- (out[, 2][[1]] +
                           #(mat[A==a, surv.tau[1], by="id"][,2][[1]]) -
                           rowSums(sapply(a, function(aa)
                           (2*(aa==a[1])-1)*(mat[get(A.name)==aa, surv.tau[1], by="id"][,2][[1]]))) -
                           (length(a)==2)*tmle.fit-
                           (length(a)==1)*(1-tmle.fit))^2
            return(sqrt(mean(ic.squared)/n))
        }

        tmle.list[[iter+1]] <- c(tmle.fit=tmle.fit, sd.eic=eval.ic(mat))
            
        if  (abs(eval.equation(mat, 0))<=eval.ic(mat)/(sqrt(n)*log(n))) {
            break
        }
    }

    return(tmle.list)    
}
