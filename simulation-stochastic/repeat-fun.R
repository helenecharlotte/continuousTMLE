repeat.fun <- function(m, K, rescue.itt=FALSE, extra=FALSE, estimate.Gstar=FALSE,
                       run.ltmle=FALSE, run.ctmle=FALSE, run.ctmle2=TRUE,
                       misspecify.Q=FALSE, 
                       n=1000) {

    dt <- sim.data(n, seed=10011+m, censoring=TRUE, rescue=TRUE,
                   browse=FALSE,
                   K=K)

    if (run.ctmle) {

        if (!rescue.itt) {
        
            (est.list.M0 <- est.fun(copy(dt), censoring=TRUE,
                                    targeting=1, rescue=TRUE, extra=extra,
                                    estimate.Gstar=estimate.Gstar,
                                    smooth.initial=TRUE,
                                    browse3=FALSE,
                                    stochastic.A=TRUE,
                                    intervention.A0=function(L0, A0) logit(1*(A0==0)),
                                    browse0=FALSE, misspecify.Q=misspecify.Q))

            (est.list.M1 <- est.fun(copy(dt), censoring=TRUE,
                                    targeting=1, rescue=TRUE, extra=extra,
                                    estimate.Gstar=estimate.Gstar,
                                    smooth.initial=TRUE,
                                    browse9=FALSE,
                                    stochastic.A=TRUE,
                                    intervention.A0=function(L0, A0) logit(1*(A0==1)),
                                    browse0=FALSE, misspecify.Q=misspecify.Q))

        } else {

            (est.list.M0 <- est.fun(copy(dt), censoring=TRUE,
                                    targeting=1, rescue=TRUE,
                                    smooth.initial=TRUE,
                                    browse3=FALSE,
                                    intervention.A=function(L0, A0, L.prev, A.prev, A) cbind(logit(1)),
                                    test=1, 
                                    intervention.A0=function(L0, A0) logit(1*(A0==0)),
                                    browse0=FALSE, misspecify.Q=misspecify.Q))

            (est.list.M1 <- est.fun(copy(dt), censoring=TRUE,
                                    targeting=1, rescue=TRUE, 
                                    smooth.initial=TRUE,
                                    browse5=FALSE,
                                    intervention.A=function(L0, A0, L.prev, A.prev, A) cbind(logit(1)),
                                    test=1, 
                                    intervention.A0=function(L0, A0) logit(1*(A0==1)),
                                    browse0=FALSE, misspecify.Q=misspecify.Q))

        }

        return(list(est.list.M1, est.list.M0))

    }

    if (run.ctmle2) {

        if (!rescue.itt) {
            
            (est.list.M0 <- est.fun(copy(dt), censoring=TRUE,
                                    targeting=2, rescue=TRUE, 
                                    smooth.initial=TRUE,
                                    browse5=FALSE,
                                    stochastic.A=TRUE,
                                    intervention.A0=function(L0, A0) logit(1*(A0==0)),
                                    browse0=FALSE, misspecify.Q=misspecify.Q))

            (est.list.M1 <- est.fun(copy(dt), censoring=TRUE,
                                    targeting=2, rescue=TRUE, 
                                    smooth.initial=TRUE,
                                    browse9=FALSE,
                                    stochastic.A=TRUE,
                                    intervention.A0=function(L0, A0) logit(1*(A0==1)),
                                    browse0=FALSE, misspecify.Q=misspecify.Q))
        } else {

            (est.list.M0 <- est.fun(copy(dt), censoring=TRUE,
                                    targeting=2, rescue=TRUE, 
                                    smooth.initial=TRUE,
                                    browse3=FALSE,
                                    intervention.A=function(L0, A0, L.prev, A.prev, A) cbind(logit(1)),
                                    test=1, 
                                    intervention.A0=function(L0, A0) logit(1*(A0==0)),
                                    browse0=FALSE, misspecify.Q=misspecify.Q))

            (est.list.M1 <- est.fun(copy(dt), censoring=TRUE,
                                    targeting=2, rescue=TRUE, 
                                    smooth.initial=TRUE,
                                    browse5=FALSE,
                                    intervention.A=function(L0, A0, L.prev, A.prev, A) cbind(logit(1)),
                                    test=1, 
                                    intervention.A0=function(L0, A0) logit(1*(A0==1)),
                                    browse0=FALSE, misspecify.Q=misspecify.Q))

        }
       
        return(list(est.list.M1, est.list.M0))
    }
    
    if (run.ltmle) {

        df.wide <- as.data.frame(cbind(dt[, "L0"],
                                       dN.A0=rep(1, nrow(dt)),
                                       dt[, -c("L0", "id")]))
        
        col.order <- names(df.wide)[-1]

        Anodes <- col.order[substr(col.order, 1, 1)=="A"]
        Ynodes <- col.order[substr(col.order, 1, 1)=="Y"]
        Cnodes <- col.order[substr(col.order, 1, 1)=="C"]
        Nnodes <- col.order[substr(col.order, 1, 4)=="dN.A"]
        Lnodes <- setdiff(col.order, c("L0", Anodes, Ynodes, Cnodes))

        for (ii in Cnodes) {
            df.wide[, names(df.wide)==ii] <- BinaryToCensoring(is.censored=df.wide[, names(df.wide)==ii])
        }

        if (TRUE) {
            abar0 <- (df.wide[, Nnodes] == 1) * 0 + (df.wide[, Nnodes] == 0) * df.wide[, Anodes]
            abar1 <- (df.wide[, Nnodes] == 1) * 1 + (df.wide[, Nnodes] == 0) * df.wide[, Anodes]
            #abar <- (df.wide[, Nnodes] == 1) * 0 + (df.wide[, Nnodes] == 0) * 0
    
            #--- specify intervention step 2: when dA(dt)=1 : 
            det.g.fun <- function(data, current.node, nodes) {
                if (substr(names(data)[current.node], 1, 1)=="C") {
                    if (FALSE) {
                        C.node.index <- current.node
                        observed.C <- data[, C.node.index]
                        is.deterministic <- observed.C=="censored" | is.na(observed.C)
                        #browser()
                        prob1 <- observed.C[is.deterministic]#observed.C[is.deterministic]
                        # if ( names(data)[current.node] == "C2") browser()
                        return(list(is.deterministic = is.deterministic, prob1 = prob1))
                    }
                    #return(NULL)
                } else {                
                    N.nodes.index <- grep("dN.A", names(data))
                    N.node.index <- max(N.nodes.index[N.nodes.index < current.node]) 
                    A.node.index <- current.node
                    N <- data[, N.node.index]
                    observed.A <- data[, A.node.index]
                    is.deterministic <- N == 0 | is.na(N)
                    prob1 <- observed.A[is.deterministic]
                    #if ( names(data)[current.node] == "A2") browser()
                    return(list(is.deterministic = is.deterministic, prob1 = prob1))
                }
            }

            print("running ltmle with deterministic arule")

            r0 <- suppressMessages(try(ltmle(df.wide,#[,!names(df.wide)%in%Cnodes],##[, col.order]
                                             Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes,
                                             abar=as.matrix(abar0),
                                             #abar=rep(0, length(Anodes)), 
                                             Ynodes=Ynodes, survivalOutcome=TRUE, 
                                             deterministic.g.function=det.g.fun,
                                             estimate.time=FALSE)))

            r1 <- suppressMessages(try(ltmle(df.wide,#[,!names(df.wide)%in%Cnodes],##[, col.order]
                                             Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes,
                                             abar=as.matrix(abar1),
                                             #abar=rep(0, length(Anodes)), 
                                             Ynodes=Ynodes, survivalOutcome=TRUE, 
                                             deterministic.g.function=det.g.fun,
                                             estimate.time=FALSE)))
        } else {

            r0 <- suppressMessages(try(ltmle(df.wide,#[,!names(df.wide)%in%Cnodes],##[, col.order]
                                             Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes,
                                             #abar=as.matrix(abar),
                                             abar=rep(0, length(Anodes)), 
                                             Ynodes=Ynodes, survivalOutcome=TRUE, 
                                             #deterministic.g.function=det.g.fun,
                                             estimate.time=FALSE)))

            r1 <- suppressMessages(try(ltmle(df.wide,#[,!names(df.wide)%in%Cnodes],##[, col.order]
                                             Anodes=Anodes, Lnodes=Lnodes, Cnodes=Cnodes,
                                             #abar=as.matrix(abar),
                                             abar=rep(1, length(Anodes)), 
                                             Ynodes=Ynodes, survivalOutcome=TRUE, 
                                             #deterministic.g.function=det.g.fun,
                                             estimate.time=FALSE)))

        }
        
        if (is(r0, "try-error")) {
            ltmle.list.0 <- "ERROR"
        } else {
            ltmle.list.0 <- list(est=r0$estimates["tmle"],
                                 sd=summary(r0)$treatment$std.dev)
        }

        if (is(r1, "try-error")) {
            ltmle.list.1 <- "ERROR"
        } else {
            ltmle.list.1 <- list(est=r1$estimates["tmle"],
                                        sd=summary(r1)$treatment$std.dev)
        }

        return(list( ltmle.list.1,  ltmle.list.0))

    }

}
