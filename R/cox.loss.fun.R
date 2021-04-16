cox.loss.fun <- function(train.fit, dt, risk.set, test.set, X=NULL, lambda.cv=NULL,
                         delta.var="delta", delta.value=1, change.point=NULL) {

    tmp <- copy(dt)
   
    if (ifelse(length(change.point)>0, change.point>0, FALSE)) {
        tmp <- rbind(tmp, tmp)[order(id)]
        tmp[, time.period:=1:.N, by="id"]
    } 
    
    if (any(class(train.fit)%in%c("coxnet", "cv.glmnet", "glmnet"))) {
        tmp[, fit.lp:=predict(train.fit, type="link", newx=X, s=lambda.cv)]
    } else {
        tmp[, fit.lp:=predict(train.fit, type="lp", newdata=tmp)]
    }

    tmp[, risk:=0]
    tmp[id%in%risk.set, risk:=1]

    tmp <- tmp[rev(order(time))]
   
    if (ifelse(length(change.point)>0, change.point>0, FALSE)) {
        tmp[, term2:=cumsum(risk*exp(fit.lp)), by="time.period"]
    } else {
        tmp[, term2:=cumsum(risk*exp(fit.lp))]
    }
    
    tmp[term2==0, term2:=1]

    if (ifelse(length(change.point)>0, change.point>0, FALSE)) {
        return(-sum(tmp[id%in%test.set &
                        ((time.period==1 & time<=change.point) | (time.period==2 & time>change.point)),
        (get(delta.var)==delta.value)*(fit.lp - log(term2))]))
    } else {
        return(-sum(tmp[id%in%test.set, (get(delta.var)==delta.value)*(fit.lp - log(term2))]))
    }
    
}
