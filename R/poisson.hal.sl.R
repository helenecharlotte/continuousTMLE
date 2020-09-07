poisson.hal.sl <- function(mat, dt, X=NULL, delta.outcome=1, cols.obs=cols.obs,
                           cut.covars=5, cut.time=10, cut.time.A=4,
                           cut.L1.A=5, cut.L.interaction=5,
                           covars=c("L1", "L2", "L3"),
                           poisson.cens=FALSE, lambda.cv=NULL,
                           penalize.time=TRUE, adjust.penalization=TRUE,
                           V=10
                           ) {

    n <- nrow(dt)
    unique.times <- sort(unique(dt[, time]))
    cv.split <- matrix(sample(1:n, size=n), ncol=V)

    for (vv in 1:V) {
    
        test.set <- cv.split[,vv]#sample(1:n, floor(n/10))
        train.set <- dt[, id][!dt[, id] %in% test.set]

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
            print("no penalization of coefficients for time indicators (main effects)")
            penalty.factor <- rep(1, ncol(X.obs))
            penalty.factor[1:cut.time] <- 0
        }

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


        #--- add cox:

        dt.train <- dt[id %in% train.set] 
        
                
        if (length(change.point)>0) { #-- if there is a change-point:

            print("estimate time-varying hazard")

            dt.train[, time.indicator:=(time<=t0)]

            dt.train2 <- rbind(dt.train, dt.train)[order(id)]
            dt.train2[, period:=1:.N, by="id"]

            #dt.train2 <- dt.train2[!time.indicator | period==1]
            dt.train2[period==1, `:=`(tstart=0, tstop=(time<=t0)*time+(time>t0)*t0)]
            dt.train2[period==2, `:=`(tstart=t0, tstop=time)]
            dt.train2[period==1 & !time.indicator, delta:=0]

            mod1 <- as.character(outcome.model)
            mod2 <- paste0(gsub(substr(mod1[2], which(strsplit(mod1[2], "")[[1]]=="(")+1,
                                       which(strsplit(mod1[2], "")[[1]]==",")-1), "tstart, tstop", mod1[2]),
                           "~", 
                           gsub("\\+A", "", gsub(" ", "", paste0("I((period==1)&(", A.name, "==1))", mod.period1,
                                                                 " + I((period==2)&(", A.name, "==1))", mod.period2, " + ",
                                                                 mod1[3]))))

            fit.cox.vv <- coxph(formula(mod2), data=dt.train2[!time.indicator | period==1])

            mat[id %in% test.set, fit.lambda.cox.vv:=predict(fit.cox.vv, newdata=mat[id %in% test.set],
                                                             type="risk")]
                   
                    
        } else { #-- if there is no change-point:
            fit.cox.vv <- coxph(as.formula(deparse(outcome.model)), #outcome.model,
                                data=dt.train)

            mat[, L1.squared:=L1^2]
            mat[, L2.squared:=L2^2]
            mat[, L3.squared:=L3^2]

            mat[id %in% test.set, fit.lambda.cox.vv:=predict(fit.cox.vv, newdata=mat[id %in% test.set],
                                                             type="risk")]
        }


    }

    #-- use nelson-aalen for baseline hazard;
    if (TRUE) {
        fit.baseline <- coxph(Surv(time, delta == 1)~1,
                              data=dt)
        bhaz.cox <- rbind(data.table(time=0, hazard=0),
                          merge(data.table(time=unique.times),
                                setDT(basehaz(fit.baseline, centered=centered)),
                                by="time", all.x=TRUE))
        bhaz.cox[, "dhaz0":=c(0, diff(hazard))]
        bhaz.cox[is.na(dhaz0), "dhaz0":=0]
        setnames(bhaz.cox, "hazard", "chaz0")
        bhaz.cox[, "chaz0":=cumsum(dhaz0)]
    } else {
        setnames(bhaz.cox, "dhaz", "dhaz0")
    }

    mat1 <- merge(mat, bhaz.cox[, c("time", "dhaz0"), with=FALSE], by="time")
                
    mat1[, fit.pois.dLambda.vv:=fit.lambda.vv*tdiff]
    mat1[, surv.Lambda.pois.vv:=cumsum(fit.pois.dLambda.vv), by=c("id", "A")]
    mat1[, surv.Lambda.cox.vv:=cumsum(dhaz0*fit.cox), by=c("id", "A")]
    mat1[, surv.t.pois.vv:=exp(-cumsum(fit.pois.dLambda.vv)), by=c("id", "A")]


    mat2 <- data.table(id=mat1[time==time.obs, id],
                       delta=mat1[time==time.obs, delta.obs],
                       M.cox=mat1[time==time.obs, delta.obs - surv.Lambda.cox.vv],
                       M.pois=mat1[time==time.obs, delta.obs - surv.Lambda.pois.vv])
    mat2[, d.cox:=sign(M.cox)*sqrt(-2*(M.cox+delta*log(delta-M.cox)))]
    mat2[, d.pois:=sign(M.pois)*sqrt(-2*(M.pois+delta*log(delta-M.pois)))]
    mat2[is.na(d.cox), d.cox:=0]

    return(mat2[, sum(d.cox^2)]>mat2[, sum(d.pois^2)])
}
