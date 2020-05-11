cov.fun <- function(x1, ltmle=FALSE) {
    return(mean(unlist(lapply(x1[-1], function(x, psi0=x1[[1]][[1]]) {
        if (is.list(x)) {
            if (!is.na(x[[length(x)]][[1]])) {
                x <- lapply(x, function(x1) if (!ltmle & length(x1)==1) return(NULL) else return(x1))
                if (!ltmle & is.numeric(x[[length(x)]][1])) {
                    return(1*(x[[length(x)]][1]-1.96*x[[length(x)]][3] <= psi0 &
                              psi0 <= x[[length(x)]][1]+1.96*x[[length(x)]][3]))
                } else if (ltmle) {
                    return(1*(x$est-1.96*x$sd <= psi0 &
                              psi0 <= x$est+1.96*x$sd))
                }
            }
        }
    }))))
}

cov.fun2 <- function(x1, x0, ltmle=FALSE) {

    true.diff <- x1[[1]][[1]] - x0[[1]][[1]]

    if (ltmle) {
        est.x1 <- unlist(lapply(x1[-1], function(x) {x$est}))
        est.x0 <- unlist(lapply(x0[-1], function(x) {x$est}))
        sd.x1 <- unlist(lapply(x1[-1], function(x) {x$sd}))
        sd.x0 <- unlist(lapply(x0[-1], function(x) {x$sd}))
    } else { #if (is.list(x1) & is.list(x0)) {
        #if (!is.na(x0[[length(x0)]][[1]]) & !is.na(x1[[length(x1)]][[1]])) {
        est.x1 <- unlist(lapply(x1[-1], function(x) {x[[length(x)]][1]}))
        est.x0 <- unlist(lapply(x0[-1], function(x) {x[[length(x)]][1]}))
        sd.x1 <- unlist(lapply(x1[-1], function(x) {x[[length(x)]][3]}))
        sd.x0 <- unlist(lapply(x0[-1], function(x) {x[[length(x)]][3]}))       
        #}
    }
    
    diff <- est.x1 - est.x0
    sd.diff <- sqrt(sd.x1^2 + sd.x0^2)

    return(mean(diff - 1.96*sd.diff <= true.diff &
                diff + 1.96*sd.diff >= true.diff))
}
