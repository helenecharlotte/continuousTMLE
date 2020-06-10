log.linear.targeting <- function(dt.hal, iteration=1, verbose=FALSE) {

    if (iteration>1) {
        dt.hal[, clever.covar:=get(paste0("fit.Surv.tau.", iteration-1))/
                     get(paste0("fit.Surv.", iteration-1))]
        dt.hal[, tmle.covar:=clever.weight*clever.covar]
        #-- estimate fluctation;         
        fit.tmle <- glm(as.formula(paste0("D ~ ", "tmle.covar-1",
                                          "+offset(log(RT*fit.lambda.", iteration-1, "))")),
                        family=poisson(link=log),
                        data=dt.hal[RT>0])
    } else {
        dt.hal[, clever.weight:=-( (A==1)/fit.pi - (A==0)/(1-fit.pi) ) * (time <= tau) / fit.Surv.cens]
        dt.hal[, clever.covar:=fit.Surv.tau/fit.Surv]
        dt.hal[, tmle.covar:=clever.weight*clever.covar]
        #-- estimate fluctation;         
        fit.tmle <- glm(as.formula(paste0("D ~ ", "tmle.covar-1",
                                          "+offset(log(RT*fit.lambda))")),
                        family=poisson(link=log),
                        data=dt.hal[RT>0])
    }

    if (verbose) print(fit.tmle)

    dt.hal[RT>0, (paste0("fit.lambda.", iteration)):=
                     exp(predict(fit.tmle, dt.hal[RT>0],
                                 newoffset=dt.hal[RT>0, log(RT*clever.weight)]))/RT]

    #-- updated cumulative hazard
    dt.hal[, (paste0("fit.Lambda.", iteration)):=cumsum(diff*get((paste0("fit.lambda.", iteration)))), by=covars]

    dt.hal[, (paste0("fit.Surv.", iteration)):=exp(-get(paste0("fit.Lambda.", iteration)))]

    dt.hal[, (paste0("fit.Surv.tau.", iteration)):=get(paste0("fit.Surv.", iteration))[.N], by=covars]

} 
