lebesgue.loss.fun <- function(train.fit, dt, risk.set=NULL, test.set, time.var=NULL, X=NULL, lambda.cv=NULL,
                              delta.var="delta", delta.value=1) {

    tmp <- copy(dt)

    tmp[id %in% test.set, fit.lambda:=exp(predict(train.fit, X[tmp$id %in% test.set,],
                                                  newoffset=0, s=lambda.cv))]
    tmp[id %in% test.set, fit.dLambda:=fit.lambda*risk.time]
    tmp[id %in% test.set, fit.Lambda:=cumsum(fit.dLambda), by="id"]

    tmp[id %in% test.set, dN:=1*(time.obs==grid.time & get(delta.var)==delta.value)]

    return(tmp[id %in% test.set & time.obs==grid.time, -sum(log(fit.lambda)*dN - fit.Lambda)])
}
