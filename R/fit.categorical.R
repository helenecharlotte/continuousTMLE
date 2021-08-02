### fit.categorical.R --- 
#----------------------------------------------------------------------
## Author: Helene C. W. Rytgaard, Thomas Alexander Gerds
## Created: Jul 28 2021 (11:58) 
## Version: 
## Last-Updated: Jul 29 2021 (18:31) 
##           By: Thomas Alexander Gerds
##     Update #: 21
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @description Pseudo hazard regression model for ordinal categoricals.
#' @details Suppose ordinal outcome \code{Y} takes three values \eqn{a,b,c,...,z}.
#' The pseudo hazard regression model specifies: \deqn{logit{P(Y=y|Y\ge y,X)}=\beta_y*X}
#' where \eqn{X} is a vector of covariates and \eqn{\beta_y} a vector of regression coefficients.
#' Specifically, the predicted probabilities from this model are
#' \deqn{P(Y=a|X)=\expit(\beta_a * X)}
#' \deqn{P(Y=b|X)=(1-P(Y=a|X))* \expit(\beta_b * X)}
#' \deqn{P(Y=c|X)=(1-P(Y=a|X))*(1-P(Y=b|X,Y\ge b))* \expit(\beta_b * X)}
#' and so on.
#' @title Pseudo hazard regression model for ordinal categoricals.
#' @param formula Object of class \code{"formula"} that describes the
#'                model to be fitted. The right hand side should be a
#'                categorical factor variable or be coercible to factor.
#' @param data Optional data.frame in which argument \code{formula} is evaluated.
#' @param ... Additional arguments that are passed to the call of \code{glm}
#'            which is used to estimate the model.
#' @examples
#' library(survival)
#' data(pbc)
#' f <- fit.categorical(edema~age+sex,data=pbc)
#' predict(f,newdata=pbc[1:3,])
#' @export
fit.categorical <- function(formula,data,...){
    pseudo.data <- get.pseudodata(formula,data)
    # introduce all interactions
    pseudo.formula <- update(formula,"pseudo.outcome~(.)*pseudo.factor")
    # fit pseudo model
    fit <- glm(pseudo.formula, data=pseudo.data, family=binomial())
    # add the call
    fit$catfit.call <- match.call()
    # add the outcome values
    fit$y <- attr(pseudo.data,"outcome.values")
    # add a new S3 class
    class(fit) <- c("catfit",class(fit))
    # return the fit
    return(fit)
}

get.pseudodata <- function(formula,data,atrisk=TRUE){
    d <- model.frame(formula,data)
    # should check that outcome is categorical 
    N <- NROW(d)
    d$id <- 1:N
    Y <- d[[1]]
    outcome <- names(d)[[1]]
    # HOV: need to deal with factors and factor levels
    # that cannot be coerced to numeric
    if (is.factor(Y))
        outcome.values <- levels(Y)
    else
        outcome.values <- sort(unique(Y))
    pseudo.data <- merge(d,data.frame(id=rep(1:N,each=length(outcome.values)),
                                      pseudo.aa=rep(outcome.values,N)),by="id")
    # should check that no variable in the original data has such name 
    pseudo.data$pseudo.outcome <- as.numeric(pseudo.data[[outcome]]==pseudo.data$pseudo.aa)
    pseudo.data$pseudo.factor <- factor(pseudo.data$pseudo.aa)
    if (atrisk[[1]]==TRUE){
        # subset to where the subjects are 'at-risk'
        pseudo.data <- pseudo.data[pseudo.data[[outcome]]>=pseudo.data$pseudo.aa,]
    }
    # add the outcome name and the unique values
    attr(pseudo.data,"outcome") <- outcome
    attr(pseudo.data,"outcome.values") <- outcome.values
    # return data
    return(pseudo.data)
}
#----------------------------------------------------------------------
### fit.categorical.R ends here
