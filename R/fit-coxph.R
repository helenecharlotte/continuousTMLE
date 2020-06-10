fit.coxph <- function(dt, covars=c("A", "L"), which.delta=1, tau=1) {

    fit.cox <- coxph(Surv(time, delta) ~ A+L, data=dt)    

    dt.1 <- copy(dt)
    dt.0 <- copy(dt)
    dt.1[, A:=1]
    dt.1[, time:=tau]
    dt.0[, A:=0]
    dt.0[, time:=tau]

    dt[, surv.cox.A1:=exp(-predict(fit.cox, newdata=dt.1, type="expected"))]
    dt[, surv.cox.A0:=exp(-predict(fit.cox, newdata=dt.0, type="expected"))]

    return(dt[, mean(surv.cox.A0-surv.cox.A1)])
    
}
