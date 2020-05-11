cov.fun <- function(psi.hat, sd, psi0) {
    return(mean(psi.hat - 1.96*sd <= psi0 &
                psi.hat + 1.96*sd >= psi0))
}
