hal.screening <- function(covars, dt, V=5, cut.one.way=18, method.risk="test", cv.glmnet=FALSE,
                          seed=13349, grouped=TRUE, use.min=TRUE, order=1,
                          cut.time.treatment=NULL,
                          delta.var="delta", delta.value=1,
                          treatment=NULL, time.var=NULL, cut.time=5,
                          mat=NULL, browse=FALSE,
                          cut.two.way=5) {

    if (order==1) {
        fit <- fit.hal(covars=covars, dt=dt, V=V, cut.one.way=cut.one.way, method.risk=method.risk,
                       cut.time.treatment=cut.time.treatment,
                       delta.var=delta.var, delta.value=delta.value,
                       mat=mat, browse=browse,
                       cv.glmnet=cv.glmnet, treatment=treatment, time.var=time.var, cut.time=cut.time,
                       seed=seed, grouped=grouped, use.min=use.min)
        picked.covars <- covars[sapply(covars, function(covar) length(grep(covar, coef(fit$fit)@Dimnames[[1]][coef(fit$fit)@i+1]))>0)]
        if (length(picked.covars)==0) {
            for (jj in 1:10) {
                fit <- fit.hal(covars=covars, dt=dt, V=V, cut.one.way=cut.one.way*2, method.risk=method.risk,
                               cut.time.treatment=cut.time.treatment,
                               delta.var=delta.var, delta.value=delta.value,
                               mat=mat, browse=browse, lambda.cvs=fit$cve$min$lambda.cv/3,
                               cv.glmnet=cv.glmnet, treatment=treatment, time.var=time.var, cut.time=cut.time,
                               seed=seed, grouped=grouped, use.min=use.min)
                picked.covars <- covars[sapply(covars, function(covar) length(grep(covar, coef(fit$fit)@Dimnames[[1]][coef(fit$fit)@i+1]))>0)]
                if (length(picked.covars)>0) break
            }
        }
        return(picked.covars)
    } else {
        two.way <- expand.grid(c(treatment, covars), c(treatment, covars))
        fit <- fit.hal(covars=covars, dt=dt, V=V, cut.one.way=cut.one.way, method.risk=method.risk, cv.glmnet=cv.glmnet,
                       seed=seed, grouped=grouped, use.min=use.min, treatment=treatment,
                       cut.time.treatment=cut.time.treatment,
                       delta.var=delta.var, delta.value=delta.value,
                       mat=mat, browse=browse,
                       time.var=time.var, cut.time=cut.time,
                       two.way=two.way, cut.two.way=cut.two.way)$fit
        two.way.out <- two.way[apply(two.way, 1, function(row2) {
            if (row2[1]!=row2[2]) {
                (length(grep(row2[1], grep(paste0(":", row2[2]), coef(fit)@Dimnames[[1]][coef(fit)@i+1], value=TRUE)))>0)
            } else return(FALSE)
        }),, drop=FALSE]
        if (is.na(two.way.out[1,1])) two.way.out <- cbind(var1="", var2="")
        return(two.way.out)
    }
    
}
