#' @title Confidence upper bound of true coverage error or threshold selection based on cross-fit one-step corrected estimator (via grid search)
#' @name CVonestep
#' @export
#' @description Method to compute a confidence upper bound of true coverage or select a threshold (via grid search) for APAC prediction sets based on cross-fit one-step corrected estimators
#' @param A vector of population indicator. 1 for source population, 0 for target population
#' @param X data frame of covariates with each row being one observation
#' @param Y vector of dependent variable/outcome. For data from the target population (`A=0`), set the corresponding entries of `Y` to be `NA`
#' @param scores either a function assigning scores of `Y` given `X` trained using an independent dataset from source population or a vector of this function evaluated at observed `(X,Y)`, taking `NA` for observations from the target population. If it is a function, it must take input `(x,y)`, where `x` is one row of `X` (a data frame with one row) and `y` is a nonmissing value of `Y`, and output a scalar
#' @param candidate.tau a numeric vector of candidate thresholds, default to `c(scores,Inf)` (after `scores` is evaluated at observations if `scores` is a function).
#' @param error.bound desired bound on the prediction set coverage error between 0 and 1, default 0.05
#' @param conf.level desired level of confidence of low coverage error between 0.5 and 1, default to 0.95
#' @param nfolds number of folds for sample splitting, default to 5
#' @param g.control a named list containing options passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}} to estimate propensity score `g`. Must not specify `Y`, `X`, `newX` or `family`. Default to `list(SL.library=c("SL.glm","SL.gam","SL.randomForest"))`
#' @param Q.control a named list containing options passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}} to estimate conditional coverage error `Q`. Must not specify `Y`, `X`, `newX` or `family`. Default to `list(SL.library=c("SL.glm","SL.gam","SL.randomForest"))`
#' @param g.trunc truncation level of propensity score `g` from zero, default to 0.01
#' @param select.tau whether to select threshold tau (otherwise just reposrt estimates and confidence upper bounds of coverage error for all `candidate.tau`), default to `TRUE` if `length(candidate.tau)>1` and `FALSE` if `length(candidate.tau)==1`
#' @return If `select.tau==FALSE`, then a list with the following components:
#' \describe{
#' \item{`tau`}{Input tau}
#' \item{`error.CI.upper`}{The (approximate) confidence upper bound of coverage error corresponding to the input tau}
#' \item{`error.est`}{The point estimate of coverage error corresponding to the input tau}
#' }
#' Otherwise a list with the following components:
#' \describe{
#' \item{`tau`}{Selected threshold tau, the maximal tau with (approximate) confidence upper bound of coverage error lower than `error.bound`}
#' \item{`error.CI.upper`}{The (approximate) confidence upper bound of coverage error corresponding to the selected tau}
#' \item{`error.est`}{The point estimate of coverage error corresponding to the selected tau}
#' \item{`feasible.tau`}{The set of feasible thresholds tau defined by (approximate) confidence upper bounds of coverage errors being lower than `error.bound`}
#' \item{`feasible.tau.error.CI.upper`}{The (approximate) confidence upper bounds of coverage errors corresponding to `feasible.tau`}
#' \item{`feasible.tau.error.est`}{The point estimates of coverage errors corresponding to `feasible.tau`}
#' }
#' 
#' @section Warnings/Errors due to extreme candidate thresholds:
#' When extremely small/large thresholds are included in `candidata.tau`, it is common to receive warnings/errors from the machine learning algorithms used by \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}, because in such cases, almost all `Y` are included in (for small thresholds) or excluded from (for large thresholds) the corresponding prediction sets, leading to complaints from machine learning algorithms. This is usually not an issue because the resulting predictions are still quite accurate. We also strongly encourage the user to specify a lerner that can deal with such cases (e.g., `SL.glm`) in `Q.control`.
#' 
#' @examples
#' n<-100
#' expit<-function(x) 1/(1+exp(-x))
#' A<-rbinom(n,1,.5)
#' X<-data.frame(X=rnorm(n,sd=ifelse(A==1,1,.5)))
#' Y<-rbinom(n,1,expit(1+X$X))
#' scores<-dbinom(Y,1,expit(.08+1.1*X$X))
#' candidate.tau<-seq(0,.5,length.out=10)
#' CVonestep(A,X,Y,scores,candidate.tau,nfolds=2,
#'           g.control=list(SL.library="SL.glm"),
#'           Q.control=list(SL.library="SL.glm"))
CVonestep<-function(A,X,Y,scores,candidate.tau,error.bound=0.05,conf.level=.95,nfolds=5,g.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),Q.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),g.trunc=1e-2,select.tau=ifelse(length(candidate.tau)==1,FALSE,TRUE)){
    #check inputs
    if(any(!(A %in% c(0,1)))){
        stop("A must take values 0 or 1")
    }
    if(all(A==0) || all(A==1)){
        stop("A must contain both values 0 and 1")
    }
    n<-length(A)
    
    assert_that(is.data.frame(X))
    assert_that(nrow(X)==length(A))
    assert_that(length(Y)==length(A))
    assert_that(is.number(error.bound) && error.bound>=0)
    assert_that(is.number(conf.level) && conf.level>=0.5 && conf.level<1)
    
    if(is.function(scores)){
        scores<-sapply(1:n,function(i){
            if(A[i]==1){
                scores(X[i,],Y[i])
            }else{
                NA
            }
        })
    }else{
        assert_that(is.numeric(scores))
    }
    
    if(missing(candidate.tau)){
        candidate.tau<-c(scores,Inf)
    }
    
    assert_that(is.numeric(candidate.tau))
    assert_that(!anyNA(A))
    assert_that(!anyNA(X))
    if(any(is.na(Y[A==1]))){
        stop("Y in the source population (A=1) must not be NA")
    }
    if(any(is.na(scores[A==1]))){
        stop("scores in the source population (A=1) must not be NA")
    }
    
    candidate.tau<-sort(unique(candidate.tau))
    
    assert_that(is.count(nfolds),nfolds>1)
    assert_that(is.list(g.control))
    assert_that(is.list(Q.control))
    if(any(c("Y","X","newX","family") %in% names(g.control))){
        stop("Y, X, newX and family must not be specified in g.control")
    }
    if(any(c("Y","X","newX","family") %in% names(Q.control))){
        stop("Y, X newX and family must not be specified in Q.control")
    }
    assert_that(is.number(g.trunc),g.trunc>0,g.trunc<1)
    assert_that(is.flag(select.tau),!is.na(select.tau))
    
    folds<-create.folds(A,nfolds)
    fold.sizes<-sapply(folds,length)
    
    #nuisance function estimators
    Z<-do.call(cbind,lapply(candidate.tau,function(tau){
        as.numeric(scores<tau)
    }))
    g<-CV.est.g(A,X,folds,g.control,g.trunc)
    Q<-CV.est.Q(A,X,Z,folds,Q.control)
    
    #calculate coverage error estimator and asymptotic variance for each fold and tau
    sigma2.v<-psi.v<-matrix(nrow=nfolds,ncol=length(candidate.tau))
    for(tau.index in 1:length(candidate.tau)){
        for(v in 1:nfolds){
            A.train<-A[-folds[[v]]]
            A.test<-A[folds[[v]]]
            Z.test<-Z[folds[[v]],tau.index]
            g.test<-g[[v]]
            Q.test<-Q[[v]][,tau.index]
            gamma.train<-mean(A.train)
            gamma.test<-mean(A.test)
            w.test<-(1-g.test)/g.test*gamma.train/(1-gamma.train)
            
            naive.psi<-sum((1-A.test)*Q.test)/sum(1-A.test)
            IF<-ifelse(A.test,(Z.test-Q.test)*w.test/gamma.test,(Q.test-naive.psi)/(1-gamma.test))
            psi.v[v,tau.index]<-naive.psi+mean(IF)
            sigma2.v[v,tau.index]<-mean(IF^2)
        }
    }
    
    #combine folds
    psi<-t(fold.sizes)%*%psi.v/n
    psi<-pmin(pmax(psi,0),1) #project psi to [0,1]
    sigma2<-t(fold.sizes)%*%sigma2.v/n
    SE<-sqrt(sigma2/n)
    
    #calculate confidence upper bound and select tau if needed
    CI.upper<-psi+SE*qnorm(conf.level)
    if(!select.tau){
        list(tau=candidate.tau,error.CI.upper=as.numeric(CI.upper),error.est=as.numeric(psi))
    }else if(CI.upper[1]<=error.bound){
        tau.index<-find.max.true(CI.upper<error.bound)
        if(tau.index==length(CI.upper)){
            warning("Larest threshold is selected. Try a different set of candidate.tau with a larger max value.")
        }
        tau<-candidate.tau[tau.index]
        list(tau=tau,error.CI.upper=CI.upper[tau.index],error.est=psi[tau.index],
             feasible.tau=candidate.tau[1:tau.index],feasible.tau.error.CI.upper=CI.upper[1:tau.index],feasible.tau.error.est=psi[1:tau.index])
    }else{
        message("No candidate threshold tau satisfies the criterion! Try smaller candidate.tau!")
        invisible()
    }
}















# Confidence upper bound of true coverage error or threshold selection based on cross-fit one-step corrected estimator with known likelihood ratio (via grid search)
# Method to compute a confidence upper bound of true coverage or select a threshold (via grid search) for APAC prediction sets based on cross-fit one-step corrected estimators with known likelihood ratio
# A: vector of population indicator. 1 for source population, 0 for target population
# X: data frame of covariates with each row being one observation
# Y: vector of dependent variable/outcome. For data from the target population (`A=0`), set the corresponding entries of `Y` to be `NA`
# scores: either a function assigning scores of `Y` given `X` trained using an independent dataset from source population or a vector of this function evaluated at observed `(X,Y)`, taking `NA` for observations from the target population. If it is a function, it must take input `(x,y)`, where `x` is one row of `X` (a data frame with one row) and `y` is a nonmissing value of `Y`, and output a scalar
# candidate.tau: a numeric vector of candidate thresholds, default to `c(scores,Inf)` (after `scores` is evaluated at observations if `scores` is a function). If `candidate.tau` has length 1, then just compute the point estimate and confidence upper bound of true coverage error of this threshold.
# LR: known likelihood ratio function (ratio of the density of covariate `X` from target population to the density of covariate `X` from source population). Must be a function that takes in a row of `X` and outputs a non-negative number. Default to constant function 1, i.e., no covariate shift
# error.bound desired bound on the prediction set coverage error between 0 and 1, default 0.05
# conf.level: desired level of confidence of low coverage error between 0.5 and 1, default to 0.95
# nfolds: number of folds for sample splitting, default to 5
# Q.control: a named list containing options passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}} to estimate conditional coverage error `Q`. Must not specify `Y`, `X`, `newX` or `family`. Default to `list(SL.library=c("SL.glm","SL.gam","SL.randomForest"))`
# select.tau: whether to select threshold tau (otherwise just reposrt estimates and confidence upper bounds of coverage error for all `candidate.tau`), default to `TRUE` if `length(candidate.tau)>1` and `FALSE` if `length(candidate.tau)==1`
# Output:
# If `select.tau==FALSE`, then a list with the following components:
# \describe{
# \item{`tau`}{Input tau}
# \item{`error.CI.upper`}{The (approximate) confidence upper bound of coverage error corresponding to the input tau}
# \item{`error.est`}{The point estimate of coverage error corresponding to the input tau}
# }
# Otherwise a list with the following components:
# \describe{
# \item{`tau`}{Selected threshold tau, the maximal tau with (approximate) confidence upper bound of coverage error lower than `error.bound`}
# \item{`error.CI.upper`}{The (approximate) confidence upper bound of coverage error corresponding to the selected tau}
# \item{`error.est`}{The point estimate of coverage error corresponding to the selected tau}
# \item{`feasible.tau`}{The set of feasible thresholds tau defined by (approximate) confidence upper bounds of coverage errors being lower than `error.bound`}
# \item{`feasible.tau.error.CI.upper`}{The (approximate) confidence upper bounds of coverage errors corresponding to `feasible.tau`}
# \item{`feasible.tau.error.est`}{The point estimates of coverage errors corresponding to `feasible.tau`}
# }
# 
# Warnings/Errors due to extreme candidate thresholds:
# When extremely small/large thresholds are included in `candidata.tau`, it is common to receive warnings/errors from the machine learning algorithms used by \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}, because in such cases, almost all `Y` are included in (for small thresholds) or excluded from (for large thresholds) the corresponding prediction sets, leading to complaints from machine learning algorithms. This is usually not an issue because the resulting predictions are still quite accurate. We also strongly encourage the user to specify a lerner that can deal with such cases (e.g., `SL.glm`) in `Q.control`.
# 
# example:
# n<-100
# expit<-function(x) 1/(1+exp(-x))
# A<-rbinom(n,1,.5)
# X<-data.frame(X=rnorm(n,sd=ifelse(A==1,1,.5)))
# Y<-rbinom(n,1,expit(1+X$X))
# scores<-dbinom(Y,1,expit(.08+1.1*X$X))
# LR<-function(x) dnorm(x$X,sd=.5)/dnorm(x$X)
# candidate.tau<-seq(0,.5,length.out=10)
# APACpred:::CVonestep_knownLR(A,X,Y,scores,candidate.tau,LR,nfolds=2,Q.control=list(SL.library="SL.glm"))
CVonestep_knownLR<-function(A,X,Y,scores,candidate.tau,LR=function(x) 1,error.bound=0.05,conf.level=.95,nfolds=5,Q.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),select.tau=ifelse(length(candidate.tau)==1,FALSE,TRUE)){
    #check inputs
    if(any(!(A %in% c(0,1)))){
        stop("A must take values 0 or 1")
    }
    if(all(A==0) || all(A==1)){
        stop("A must contain both values 0 and 1")
    }
    n<-length(A)
    
    assert_that(is.data.frame(X))
    assert_that(nrow(X)==length(A))
    assert_that(length(Y)==length(A))
    assert_that(is.number(error.bound) && error.bound>=0)
    assert_that(is.number(conf.level) && conf.level>=0.5 && conf.level<1)
    
    if(is.function(scores)){
        scores<-sapply(1:n,function(i){
            if(A[i]==1){
                scores(X[i,],Y[i])
            }else{
                NA
            }
        })
    }else{
        assert_that(is.numeric(scores))
    }
    
    if(missing(candidate.tau)){
        candidate.tau<-c(scores,Inf)
    }
    
    assert_that(is.numeric(candidate.tau))
    assert_that(!anyNA(A))
    assert_that(!anyNA(X))
    if(any(is.na(Y[A==1]))){
        stop("Y in the source population (A=1) must not be NA")
    }
    if(any(is.na(scores[A==1]))){
        stop("scores in the source population (A=1) must not be NA")
    }
    
    candidate.tau<-sort(unique(candidate.tau))
    
    assert_that(is.count(nfolds),nfolds>1)
    assert_that(is.list(Q.control))
    if(any(c("Y","X","newX","family") %in% names(Q.control))){
        stop("Y, X newX and family must not be specified in Q.control")
    }
    
    folds<-create.folds(A,nfolds)
    fold.sizes<-sapply(folds,length)
    assert_that(is.flag(select.tau),!is.na(select.tau))
    
    #nuisance function estimators
    Z<-do.call(cbind,lapply(candidate.tau,function(tau){
        as.numeric(scores<tau)
    }))
    Q<-CV.est.Q(A,X,Z,folds,Q.control)
    w<-sapply(1:nrow(X),function(i) LR(X[i,,drop=FALSE]))
    assert_that(is.numeric(w),all(w>=0 & w<Inf))
    
    #calculate coverage error estimator and asymptotic variance for each fold and tau
    sigma2.v<-psi.v<-matrix(nrow=nfolds,ncol=length(candidate.tau))
    for(tau.index in 1:length(candidate.tau)){
        for(v in 1:nfolds){
            A.train<-A[-folds[[v]]]
            A.test<-A[folds[[v]]]
            Z.test<-Z[folds[[v]],tau.index]
            Q.test<-Q[[v]][,tau.index]
            X.test<-X[folds[[v]],,drop=FALSE]
            gamma.test<-mean(A.test)
            w.test<-w[folds[[v]]]
            
            naive.psi<-sum((1-A.test)*Q.test)/sum(1-A.test)
            IF<-ifelse(A.test,(Z.test-Q.test)*w.test/gamma.test,(Q.test-naive.psi)/(1-gamma.test))
            psi.v[v,tau.index]<-naive.psi+mean(IF)
            sigma2.v[v,tau.index]<-mean(IF^2)
        }
    }
    
    #combine folds
    psi<-t(fold.sizes)%*%psi.v/n
    psi<-pmin(pmax(psi,0),1) #project psi to [0,1]
    sigma2<-t(fold.sizes)%*%sigma2.v/n
    SE<-sqrt(sigma2/n)
    
    #calculate confidence upper bound and select tau if needed
    CI.upper<-psi+SE*qnorm(conf.level)
    if(!select.tau){
        list(tau=candidate.tau,error.CI.upper=as.numeric(CI.upper),error.est=as.numeric(psi))
    }else if(CI.upper[1]<=error.bound){
        tau.index<-find.max.true(CI.upper<error.bound)
        if(tau.index==length(CI.upper)){
            warning("Larest threshold is selected. Try a different set of candidate.tau with a larger max value.")
        }
        tau<-candidate.tau[tau.index]
        list(tau=tau,error.CI.upper=CI.upper[tau.index],error.est=psi[tau.index],
             feasible.tau=candidate.tau[1:tau.index],feasible.tau.error.CI.upper=CI.upper[1:tau.index],feasible.tau.error.est=psi[1:tau.index])
    }else{
        message("No candidate threshold tau satisfies the criterion! Try smaller candidate.tau!")
        invisible()
    }
}










#' @title Threshold selection based on cross-fit one-step corrected estimator via bisection search
#' @name CVonestep_bisec
#' @export
#' @description Method to select a threshold for APAC prediction sets based on cross-fit one-step corrected estimators via bisection serach
#' @param A vector of population indicator. 1 for source population, 0 for target population
#' @param X data frame of covariates with each row being one observation
#' @param Y vector of dependent variable/outcome. For data from the target population (`A=0`), set the corresponding entries of `Y` to be `NA`
#' @param scores either a function assigning scores of `Y` given `X` trained using an independent dataset from source population or a vector of this function evaluated at observed `(X,Y)`, taking `NA` for observations from the target population. If it is a function, it must take input `(x,y)`, where `x` is one row of `X` (a data frame with one row) and `y` is a nonmissing value of `Y`, and output a scalar
#' @param tau.interval a numeric vector specifying the initial interval of candidate thresholds for bisection search, default to `c(0,1)`, which is suitable when all scores lie in the interval \eqn{[0,1]}{[0,1]}
#' @param error.bound desired bound on the prediction set coverage error between 0 and 1, default 0.05
#' @param conf.level desired level of confidence of low coverage error between 0.5 and 1, default to 0.95
#' @param nfolds number of folds for sample splitting, default to 5
#' @param g.control a named list containing options passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}} to estimate propensity score `g`. Must not specify `Y`, `X`, `newX` or `family`. Default to `list(SL.library=c("SL.glm","SL.gam","SL.randomForest"))`
#' @param Q.control a named list containing options passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}} to estimate conditional coverage error `Q`. Must not specify `Y`, `X`, `newX` or `family`. Default to `list(SL.library=c("SL.glm","SL.gam","SL.randomForest"))`
#' @param g.trunc truncation level of propensity score `g` from zero, default to 0.01
#' @param tol stopping criterion for bisection search, default to 0.01. The search will stop if the length of the interval for thresholds is below `tol` and the function will output the end point with confidence upper bound no greater than `error.bound`
#' @return a list with the following components:
#' \describe{
#' \item{`tau`}{Selected tau}
#' \item{`error.CI.upper`}{The (approximate) confidence upper bound of coverage error corresponding to the select tau}
#' \item{`error.est`}{The point estimate of coverage error corresponding to the selected tau}
#' }
#' @section Warnings/Errors due to extreme candidate thresholds:
#' When extremely small/large thresholds are included in `candidata.tau`, it is common to receive warnings/errors from the machine learning algorithms used by \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}, because in such cases, almost all `Y` are included in (for small thresholds) or excluded from (for large thresholds) the corresponding prediction sets, leading to complaints from machine learning algorithms. This is usually not an issue because the resulting predictions are still quite accurate. We also strongly encourage the user to specify a lerner that can deal with such cases (e.g., `SL.glm`) in `Q.control`.
#' 
#' @examples
#' n<-100
#' expit<-function(x) 1/(1+exp(-x))
#' A<-rbinom(n,1,.5)
#' X<-data.frame(X=rnorm(n,sd=ifelse(A==1,1,.5)))
#' Y<-rbinom(n,1,expit(1+X$X))
#' scores<-dbinom(Y,1,expit(.08+1.1*X$X))
#' CVonestep_bisec(A,X,Y,scores,nfolds=2,
#'                 g.control=list(SL.library="SL.glm"),
#'                 Q.control=list(SL.library="SL.glm"))
CVonestep_bisec<-function(A,X,Y,scores,tau.interval=c(0,1),error.bound=0.05,conf.level=.95,nfolds=5,g.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),Q.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),g.trunc=1e-2,tol=1e-2){
    #check inputs
    if(any(!(A %in% c(0,1)))){
        stop("A must take values 0 or 1")
    }
    if(all(A==0) || all(A==1)){
        stop("A must contain both values 0 and 1")
    }
    n<-length(A)
    
    assert_that(is.data.frame(X))
    assert_that(nrow(X)==length(A))
    assert_that(length(Y)==length(A))
    assert_that(is.number(error.bound) && error.bound>=0)
    assert_that(is.number(conf.level) && conf.level>=0.5 && conf.level<1)
    
    if(is.function(scores)){
        scores<-sapply(1:n,function(i){
            if(A[i]==1){
                scores(X[i,],Y[i])
            }else{
                NA
            }
        })
    }else{
        assert_that(is.numeric(scores))
    }
    
    assert_that(is.numeric(tau.interval),length(tau.interval)==2,!anyNA(tau.interval),all(is.finite(tau.interval)))
    tau.interval<-sort(tau.interval)
    
    assert_that(!anyNA(A))
    assert_that(!anyNA(X))
    if(any(is.na(Y[A==1]))){
        stop("Y in the source population (A=1) must not be NA")
    }
    if(any(is.na(scores[A==1]))){
        stop("scores in the source population (A=1) must not be NA")
    }
    
    assert_that(is.count(nfolds),nfolds>1)
    assert_that(is.list(g.control))
    assert_that(is.list(Q.control))
    if(any(c("Y","X","newX","family") %in% names(g.control))){
        stop("Y, X, newX and family must not be specified in g.control")
    }
    if(any(c("Y","X","newX","family") %in% names(Q.control))){
        stop("Y, X newX and family must not be specified in Q.control")
    }
    assert_that(is.number(g.trunc),g.trunc>0,g.trunc<1)
    assert_that(is.number(tol),tol>0)
    
    folds<-create.folds(A,nfolds)
    fold.sizes<-sapply(folds,length)
    
    #nuisance function estimators
    Z<-do.call(cbind,lapply(tau.interval,function(tau){
        as.numeric(scores<tau)
    }))
    g<-CV.est.g(A,X,folds,g.control,g.trunc)
    Q<-CV.est.Q(A,X,Z,folds,Q.control)
    
    #calculate coverage error estimator and asymptotic variance for each fold and the end points of tau.interval
    sigma2.v<-psi.v<-matrix(nrow=nfolds,ncol=2)
    for(tau.index in 1:2){
        for(v in 1:nfolds){
            A.train<-A[-folds[[v]]]
            A.test<-A[folds[[v]]]
            Z.test<-Z[folds[[v]],tau.index]
            g.test<-g[[v]]
            Q.test<-Q[[v]][,tau.index]
            gamma.train<-mean(A.train)
            gamma.test<-mean(A.test)
            w.test<-(1-g.test)/g.test*gamma.train/(1-gamma.train)
            
            naive.psi<-sum((1-A.test)*Q.test)/sum(1-A.test)
            IF<-ifelse(A.test,(Z.test-Q.test)*w.test/gamma.test,(Q.test-naive.psi)/(1-gamma.test))
            psi.v[v,tau.index]<-naive.psi+mean(IF)
            sigma2.v[v,tau.index]<-mean(IF^2)
        }
    }
    psi<-t(fold.sizes)%*%psi.v/n
    psi<-pmin(pmax(psi,0),1) #project psi to [0,1]
    sigma2<-t(fold.sizes)%*%sigma2.v/n
    SE<-sqrt(sigma2/n)
    CI.upper<-psi+SE*qnorm(conf.level)
    
    tau.l<-tau.interval[1]
    tau.u<-tau.interval[2]
    psi.l<-as.numeric(psi[1])
    psi.u<-as.numeric(psi[2])
    CI.upper.l<-as.numeric(CI.upper[1])
    CI.upper.u<-as.numeric(CI.upper[2])
    sign.l<-sign(CI.upper.l-error.bound)
    sign.u<-sign(CI.upper.u-error.bound)
    
    if(sign.l==sign.u){
        stop("Bisection search failed! Both confidence upper bounds at tau.interval are above/below error.bound! Please rerun with a different tau.interval.")
    }
    
    #bisection search
    while(tau.u-tau.l>=tol){
        new.tau<-(tau.l+tau.u)/2
        
        Z<-matrix(as.numeric(scores<new.tau),ncol=1)
        Q<-CV.est.Q(A,X,Z,folds,Q.control)
        sigma2.v<-psi.v<-matrix(nrow=nfolds,ncol=1)
        for(v in 1:nfolds){
            A.train<-A[-folds[[v]]]
            A.test<-A[folds[[v]]]
            Z.test<-Z[folds[[v]],1]
            g.test<-g[[v]]
            Q.test<-Q[[v]][,1]
            gamma.train<-mean(A.train)
            gamma.test<-mean(A.test)
            w.test<-(1-g.test)/g.test*gamma.train/(1-gamma.train)
            
            naive.psi<-sum((1-A.test)*Q.test)/sum(1-A.test)
            IF<-ifelse(A.test,(Z.test-Q.test)*w.test/gamma.test,(Q.test-naive.psi)/(1-gamma.test))
            psi.v[v,1]<-naive.psi+mean(IF)
            sigma2.v[v,1]<-mean(IF^2)
        }
        
        new.psi<-as.numeric(t(fold.sizes)%*%psi.v/n)
        new.psi<-pmin(pmax(new.psi,0),1) #project psi to [0,1]
        new.sigma2<-as.numeric(t(fold.sizes)%*%sigma2.v/n)
        new.SE<-sqrt(new.sigma2/n)
        new.CI.upper<-new.psi+new.SE*qnorm(conf.level)
        new.sign<-sign(new.CI.upper-error.bound)
        
        if(new.sign==sign.l){
            tau.l<-new.tau
            psi.l<-new.psi
            CI.upper.l<-new.CI.upper
            sign.l<-new.sign
        }else{
            tau.u<-new.tau
            psi.u<-new.psi
            CI.upper.u<-new.CI.upper
            sign.u<-new.sign
        }
    }
    
    if(sign.u<0){
        list(tau=tau.u,error.CI.upper=CI.upper.u,error.est=psi.u)
    }else{
        list(tau=tau.l,error.CI.upper=CI.upper.l,error.est=psi.l)
    }
}








# Threshold selection based on cross-fit one-step corrected estimator via bisection search with known likelihood ratio
# Method to select a threshold for APAC prediction sets based on cross-fit one-step corrected estimators via bisection search with known likelihood ratio
# A: vector of population indicator. 1 for source population, 0 for target population
# X: data frame of covariates with each row being one observation
# Y: vector of dependent variable/outcome. For data from the target population (`A=0`), set the corresponding entries of `Y` to be `NA`
# scores: either a function assigning scores of `Y` given `X` trained using an independent dataset from source population or a vector of this function evaluated at observed `(X,Y)`, taking `NA` for observations from the target population. If it is a function, it must take input `(x,y)`, where `x` is one row of `X` (a data frame with one row) and `y` is a nonmissing value of `Y`, and output a scalar
# tau.interval: a numeric vector specifying the initial interval of candidate thresholds for bisection search, default to `c(0,1)`, which is suitable when all scores lie in the interval \eqn{[0,1]}{[0,1]}
# LR: known likelihood ratio function (ratio of the density of covariate `X` from target population to the density of covariate `X` from source population). Must be a function that takes in a row of `X` and outputs a non-negative number. Default to constant function 1, i.e., no covariate shift
# error.bound: desired bound on the prediction set coverage error between 0 and 1, default 0.05
# conf.level: desired level of confidence of low coverage error between 0.5 and 1, default to 0.95
# nfolds: number of folds for sample splitting, default to 5
# Q.control: a named list containing options passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}} to estimate conditional coverage error `Q`. Must not specify `Y`, `X`, `newX` or `family`. Default to `list(SL.library=c("SL.glm","SL.gam","SL.randomForest"))`
# @param tol stopping criterion for bisection search, default to 0.01. The search will stop if the length of the interval for thresholds is below `tol` and the function will output the end point with confidence upper bound no greater than `error.bound`
# 
# Output:
# a list with the following components:
# \describe{
# \item{`tau`}{Selected tau}
# \item{`error.CI.upper`}{The (approximate) confidence upper bound of coverage error corresponding to the select tau}
# \item{`error.est`}{The point estimate of coverage error corresponding to the selected tau}
# }
# 
# Warnings/Errors due to extreme candidate thresholds:
# When extremely small/large thresholds are included in `candidata.tau`, it is common to receive warnings/errors from the machine learning algorithms used by \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}, because in such cases, almost all `Y` are included in (for small thresholds) or excluded from (for large thresholds) the corresponding prediction sets, leading to complaints from machine learning algorithms. This is usually not an issue because the resulting predictions are still quite accurate. We also strongly encourage the user to specify a lerner that can deal with such cases (e.g., `SL.glm`) in `Q.control`.
# 
# example:
# n<-100
# expit<-function(x) 1/(1+exp(-x))
# A<-rbinom(n,1,.5)
# X<-data.frame(X=rnorm(n,sd=ifelse(A==1,1,.5)))
# Y<-rbinom(n,1,expit(1+X$X))
# scores<-dbinom(Y,1,expit(.08+1.1*X$X))
# LR<-function(x) dnorm(x$X,sd=.5)/dnorm(x$X)
# APACpred:::CVonestep_bisec(A,X,Y,scores,Q.control=list(SL.library="SL.glm"))
CVonestep_knownLR_bisec<-function(A,X,Y,scores,tau.interval=c(0,1),LR=function(x) 1,error.bound=0.05,conf.level=.95,nfolds=5,Q.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),tol=1e-2){
    #check inputs
    if(any(!(A %in% c(0,1)))){
        stop("A must take values 0 or 1")
    }
    if(all(A==0) || all(A==1)){
        stop("A must contain both values 0 and 1")
    }
    n<-length(A)
    
    assert_that(is.data.frame(X))
    assert_that(nrow(X)==length(A))
    assert_that(length(Y)==length(A))
    assert_that(is.number(error.bound) && error.bound>=0)
    assert_that(is.number(conf.level) && conf.level>=0.5 && conf.level<1)
    
    if(is.function(scores)){
        scores<-sapply(1:n,function(i){
            if(A[i]==1){
                scores(X[i,],Y[i])
            }else{
                NA
            }
        })
    }else{
        assert_that(is.numeric(scores))
    }
    
    assert_that(is.numeric(tau.interval),length(tau.interval)==2,!anyNA(tau.interval),all(is.finite(tau.interval)))
    tau.interval<-sort(tau.interval)
    
    assert_that(!anyNA(A))
    assert_that(!anyNA(X))
    if(any(is.na(Y[A==1]))){
        stop("Y in the source population (A=1) must not be NA")
    }
    if(any(is.na(scores[A==1]))){
        stop("scores in the source population (A=1) must not be NA")
    }
    
    assert_that(is.count(nfolds),nfolds>1)
    assert_that(is.list(Q.control))
    if(any(c("Y","X","newX","family") %in% names(Q.control))){
        stop("Y, X newX and family must not be specified in Q.control")
    }
    assert_that(is.number(tol),tol>0)
    
    folds<-create.folds(A,nfolds)
    fold.sizes<-sapply(folds,length)
    
    #nuisance function estimators
    Z<-do.call(cbind,lapply(tau.interval,function(tau){
        as.numeric(scores<tau)
    }))
    Q<-CV.est.Q(A,X,Z,folds,Q.control)
    w<-sapply(1:nrow(X),function(i) LR(X[i,,drop=FALSE]))
    assert_that(is.numeric(w),all(w>=0 & w<Inf))
    
    #calculate coverage error estimator and asymptotic variance for each fold and the end points of tau
    sigma2.v<-psi.v<-matrix(nrow=nfolds,ncol=2)
    for(tau.index in 1:2){
        for(v in 1:nfolds){
            A.train<-A[-folds[[v]]]
            A.test<-A[folds[[v]]]
            Z.test<-Z[folds[[v]],tau.index]
            Q.test<-Q[[v]][,tau.index]
            gamma.train<-mean(A.train)
            gamma.test<-mean(A.test)
            w.test<-w[folds[[v]]]
            
            naive.psi<-sum((1-A.test)*Q.test)/sum(1-A.test)
            IF<-ifelse(A.test,(Z.test-Q.test)*w.test/gamma.test,(Q.test-naive.psi)/(1-gamma.test))
            psi.v[v,tau.index]<-naive.psi+mean(IF)
            sigma2.v[v,tau.index]<-mean(IF^2)
        }
    }
    psi<-t(fold.sizes)%*%psi.v/n
    psi<-pmin(pmax(psi,0),1) #project psi to [0,1]
    sigma2<-t(fold.sizes)%*%sigma2.v/n
    SE<-sqrt(sigma2/n)
    CI.upper<-psi+SE*qnorm(conf.level)
    
    tau.l<-tau.interval[1]
    tau.u<-tau.interval[2]
    psi.l<-as.numeric(psi[1])
    psi.u<-as.numeric(psi[2])
    CI.upper.l<-as.numeric(CI.upper[1])
    CI.upper.u<-as.numeric(CI.upper[2])
    sign.l<-sign(CI.upper.l-error.bound)
    sign.u<-sign(CI.upper.u-error.bound)
    
    if(sign.l==sign.u){
        stop("Bisection search failed! Both confidence upper bounds at tau.interval are above/below error.bound! Please rerun with a different tau.interval.")
    }
    
    #bisection search
    while(tau.u-tau.l>=tol){
        new.tau<-(tau.l+tau.u)/2
        
        Z<-matrix(as.numeric(scores<new.tau),ncol=1)
        Q<-CV.est.Q(A,X,Z,folds,Q.control)
        sigma2.v<-psi.v<-matrix(nrow=nfolds,ncol=1)
        for(v in 1:nfolds){
            A.train<-A[-folds[[v]]]
            A.test<-A[folds[[v]]]
            Z.test<-Z[folds[[v]],1]
            Q.test<-Q[[v]][,1]
            gamma.train<-mean(A.train)
            gamma.test<-mean(A.test)
            w.test<-w[folds[[v]]]
            
            naive.psi<-sum((1-A.test)*Q.test)/sum(1-A.test)
            IF<-ifelse(A.test,(Z.test-Q.test)*w.test/gamma.test,(Q.test-naive.psi)/(1-gamma.test))
            psi.v[v,1]<-naive.psi+mean(IF)
            sigma2.v[v,1]<-mean(IF^2)
        }
        
        new.psi<-as.numeric(t(fold.sizes)%*%psi.v/n)
        new.psi<-pmin(pmax(new.psi,0),1) #project psi to [0,1]
        new.sigma2<-as.numeric(t(fold.sizes)%*%sigma2.v/n)
        new.SE<-sqrt(new.sigma2/n)
        new.CI.upper<-new.psi+new.SE*qnorm(conf.level)
        new.sign<-sign(new.CI.upper-error.bound)
        
        if(new.sign==sign.l){
            tau.l<-new.tau
            psi.l<-new.psi
            CI.upper.l<-new.CI.upper
            sign.l<-new.sign
        }else{
            tau.u<-new.tau
            psi.u<-new.psi
            CI.upper.u<-new.CI.upper
            sign.u<-new.sign
        }
    }
    
    if(sign.u<0){
        list(tau=tau.u,error.CI.upper=CI.upper.u,error.est=psi.u)
    }else{
        list(tau=tau.l,error.CI.upper=CI.upper.l,error.est=psi.l)
    }
}























# Similar to CVonestep_bisec but does bisection search over observed scores
# 
# example:
# n<-100
# expit<-function(x) 1/(1+exp(-x))
# A<-rbinom(n,1,.5)
# X<-data.frame(X=rnorm(n,sd=ifelse(A==1,1,.5)))
# Y<-rbinom(n,1,expit(1+X$X))
# scores<-dbinom(Y,1,expit(.08+1.1*X$X))
# APACpred:::CVonestep_bisec2(A,X,Y,scores,nfolds=2,g.control=list(SL.library="SL.glm"),Q.control=list(SL.library="SL.glm"))
CVonestep_bisec2<-function(A,X,Y,scores,error.bound=0.05,conf.level=.95,nfolds=5,g.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),Q.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),g.trunc=1e-2,tol=1e-2){
    #check inputs
    if(any(!(A %in% c(0,1)))){
        stop("A must take values 0 or 1")
    }
    if(all(A==0) || all(A==1)){
        stop("A must contain both values 0 and 1")
    }
    n<-length(A)
    
    assert_that(is.data.frame(X))
    assert_that(nrow(X)==length(A))
    assert_that(length(Y)==length(A))
    assert_that(is.number(error.bound) && error.bound>=0)
    assert_that(is.number(conf.level) && conf.level>=0.5 && conf.level<1)
    
    if(is.function(scores)){
        scores<-sapply(1:n,function(i){
            if(A[i]==1){
                scores(X[i,],Y[i])
            }else{
                NA
            }
        })
    }else{
        assert_that(is.numeric(scores))
    }
    
    assert_that(!anyNA(A))
    assert_that(!anyNA(X))
    if(any(is.na(Y[A==1]))){
        stop("Y in the source population (A=1) must not be NA")
    }
    if(any(is.na(scores[A==1]))){
        stop("scores in the source population (A=1) must not be NA")
    }
    
    assert_that(is.count(nfolds),nfolds>1)
    assert_that(is.list(g.control))
    assert_that(is.list(Q.control))
    if(any(c("Y","X","newX","family") %in% names(g.control))){
        stop("Y, X, newX and family must not be specified in g.control")
    }
    if(any(c("Y","X","newX","family") %in% names(Q.control))){
        stop("Y, X newX and family must not be specified in Q.control")
    }
    assert_that(is.number(g.trunc),g.trunc>0,g.trunc<1)
    assert_that(is.number(tol),tol>0)
    
    folds<-create.folds(A,nfolds)
    fold.sizes<-sapply(folds,length)
    
    sorted.scores<-c(sort(scores),Inf)
    tau.index.interval<-c(1,length(sorted.scores))
    
    #nuisance function estimators
    Z<-do.call(cbind,lapply(tau.index.interval,function(tau.index){
        as.numeric(scores<sorted.scores[tau.index])
    }))
    g<-CV.est.g(A,X,folds,g.control,g.trunc)
    Q<-CV.est.Q(A,X,Z,folds,Q.control)
    
    #calculate coverage error estimator and asymptotic variance for each fold and the end points of tau.interval
    sigma2.v<-psi.v<-matrix(nrow=nfolds,ncol=2)
    for(tau.index in 1:2){
        for(v in 1:nfolds){
            A.train<-A[-folds[[v]]]
            A.test<-A[folds[[v]]]
            Z.test<-Z[folds[[v]],tau.index]
            g.test<-g[[v]]
            Q.test<-Q[[v]][,tau.index]
            gamma.train<-mean(A.train)
            gamma.test<-mean(A.test)
            w.test<-(1-g.test)/g.test*gamma.train/(1-gamma.train)
            
            naive.psi<-sum((1-A.test)*Q.test)/sum(1-A.test)
            IF<-ifelse(A.test,(Z.test-Q.test)*w.test/gamma.test,(Q.test-naive.psi)/(1-gamma.test))
            psi.v[v,tau.index]<-naive.psi+mean(IF)
            sigma2.v[v,tau.index]<-mean(IF^2)
        }
    }
    psi<-t(fold.sizes)%*%psi.v/n
    psi<-pmin(pmax(psi,0),1) #project psi to [0,1]
    sigma2<-t(fold.sizes)%*%sigma2.v/n
    SE<-sqrt(sigma2/n)
    CI.upper<-psi+SE*qnorm(conf.level)
    
    tau.index.l<-tau.index.interval[1]
    tau.index.u<-tau.index.interval[2]
    psi.l<-as.numeric(psi[1])
    psi.u<-as.numeric(psi[2])
    CI.upper.l<-as.numeric(CI.upper[1])
    CI.upper.u<-as.numeric(CI.upper[2])
    sign.l<-sign(CI.upper.l-error.bound)
    sign.u<-sign(CI.upper.u-error.bound)
    
    if(sign.l==sign.u){
        stop("Bisection search failed! Both confidence upper bounds at tau.interval are above/below error.bound! Please rerun with a different tau.interval.")
    }
    
    #bisection search
    while(sorted.scores[tau.index.u]-sorted.scores[tau.index.l]>=tol && tau.index.u-tau.index.l>1){
        new.tau.index<-round((tau.index.l+tau.index.u)/2)
        
        Z<-matrix(as.numeric(scores<sorted.scores[new.tau.index]),ncol=1)
        Q<-CV.est.Q(A,X,Z,folds,Q.control)
        sigma2.v<-psi.v<-matrix(nrow=nfolds,ncol=1)
        for(v in 1:nfolds){
            A.train<-A[-folds[[v]]]
            A.test<-A[folds[[v]]]
            Z.test<-Z[folds[[v]],1]
            g.test<-g[[v]]
            Q.test<-Q[[v]][,1]
            gamma.train<-mean(A.train)
            gamma.test<-mean(A.test)
            w.test<-(1-g.test)/g.test*gamma.train/(1-gamma.train)
            
            naive.psi<-sum((1-A.test)*Q.test)/sum(1-A.test)
            IF<-ifelse(A.test,(Z.test-Q.test)*w.test/gamma.test,(Q.test-naive.psi)/(1-gamma.test))
            psi.v[v,1]<-naive.psi+mean(IF)
            sigma2.v[v,1]<-mean(IF^2)
        }
        
        new.psi<-as.numeric(t(fold.sizes)%*%psi.v/n)
        new.psi<-pmin(pmax(new.psi,0),1) #project psi to [0,1]
        new.sigma2<-as.numeric(t(fold.sizes)%*%sigma2.v/n)
        new.SE<-sqrt(new.sigma2/n)
        new.CI.upper<-new.psi+new.SE*qnorm(conf.level)
        new.sign<-sign(new.CI.upper-error.bound)
        
        if(new.sign==sign.l){
            tau.index.l<-new.tau.index
            psi.l<-new.psi
            CI.upper.l<-new.CI.upper
            sign.l<-new.sign
        }else{
            tau.index.u<-new.tau.index
            psi.u<-new.psi
            CI.upper.u<-new.CI.upper
            sign.u<-new.sign
        }
    }
    
    if(sign.u<0){
        list(tau=sorted.scores[tau.index.u],error.CI.upper=CI.upper.u,error.est=psi.u)
    }else{
        list(tau=sorted.scores[tau.index.l],error.CI.upper=CI.upper.l,error.est=psi.l)
    }
}









# Similar to CVonestep_knownLR_bisec, but does bisection search over observed scores
# 
# example:
# n<-100
# expit<-function(x) 1/(1+exp(-x))
# A<-rbinom(n,1,.5)
# X<-data.frame(X=rnorm(n,sd=ifelse(A==1,1,.5)))
# Y<-rbinom(n,1,expit(1+X$X))
# scores<-dbinom(Y,1,expit(.08+1.1*X$X))
# LR<-function(x) dnorm(x$X,sd=.5)/dnorm(x$X)
# APACpred:::CVonestep_knownLR_bisec2(A,X,Y,scores,Q.control=list(SL.library="SL.glm"))
CVonestep_knownLR_bisec2<-function(A,X,Y,scores,LR=function(x) 1,error.bound=0.05,conf.level=.95,nfolds=5,Q.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),tol=1e-2){
    #check inputs
    if(any(!(A %in% c(0,1)))){
        stop("A must take values 0 or 1")
    }
    if(all(A==0) || all(A==1)){
        stop("A must contain both values 0 and 1")
    }
    n<-length(A)
    
    assert_that(is.data.frame(X))
    assert_that(nrow(X)==length(A))
    assert_that(length(Y)==length(A))
    assert_that(is.number(error.bound) && error.bound>=0)
    assert_that(is.number(conf.level) && conf.level>=0.5 && conf.level<1)
    
    if(is.function(scores)){
        scores<-sapply(1:n,function(i){
            if(A[i]==1){
                scores(X[i,],Y[i])
            }else{
                NA
            }
        })
    }else{
        assert_that(is.numeric(scores))
    }
    
    assert_that(!anyNA(A))
    assert_that(!anyNA(X))
    if(any(is.na(Y[A==1]))){
        stop("Y in the source population (A=1) must not be NA")
    }
    if(any(is.na(scores[A==1]))){
        stop("scores in the source population (A=1) must not be NA")
    }
    
    assert_that(is.count(nfolds),nfolds>1)
    assert_that(is.list(Q.control))
    if(any(c("Y","X","newX","family") %in% names(Q.control))){
        stop("Y, X newX and family must not be specified in Q.control")
    }
    assert_that(is.number(tol),tol>0)
    
    folds<-create.folds(A,nfolds)
    fold.sizes<-sapply(folds,length)
    
    sorted.scores<-c(sort(scores),Inf)
    tau.index.interval<-c(1,length(sorted.scores))
    
    #nuisance function estimators
    Z<-do.call(cbind,lapply(tau.index.interval,function(tau.index){
        as.numeric(scores<sorted.scores[tau.index])
    }))
    Q<-CV.est.Q(A,X,Z,folds,Q.control)
    w<-sapply(1:nrow(X),function(i) LR(X[i,,drop=FALSE]))
    assert_that(is.numeric(w),all(w>=0 & w<Inf))
    
    #calculate coverage error estimator and asymptotic variance for each fold and the end points of tau.interval
    sigma2.v<-psi.v<-matrix(nrow=nfolds,ncol=2)
    for(tau.index in 1:2){
        for(v in 1:nfolds){
            A.train<-A[-folds[[v]]]
            A.test<-A[folds[[v]]]
            Z.test<-Z[folds[[v]],tau.index]
            Q.test<-Q[[v]][,tau.index]
            gamma.train<-mean(A.train)
            gamma.test<-mean(A.test)
            w.test<-w[folds[[v]]]
            
            naive.psi<-sum((1-A.test)*Q.test)/sum(1-A.test)
            IF<-ifelse(A.test,(Z.test-Q.test)*w.test/gamma.test,(Q.test-naive.psi)/(1-gamma.test))
            psi.v[v,tau.index]<-naive.psi+mean(IF)
            sigma2.v[v,tau.index]<-mean(IF^2)
        }
    }
    psi<-t(fold.sizes)%*%psi.v/n
    psi<-pmin(pmax(psi,0),1) #project psi to [0,1]
    sigma2<-t(fold.sizes)%*%sigma2.v/n
    SE<-sqrt(sigma2/n)
    CI.upper<-psi+SE*qnorm(conf.level)
    
    tau.index.l<-tau.index.interval[1]
    tau.index.u<-tau.index.interval[2]
    psi.l<-as.numeric(psi[1])
    psi.u<-as.numeric(psi[2])
    CI.upper.l<-as.numeric(CI.upper[1])
    CI.upper.u<-as.numeric(CI.upper[2])
    sign.l<-sign(CI.upper.l-error.bound)
    sign.u<-sign(CI.upper.u-error.bound)
    
    if(sign.l==sign.u){
        stop("Bisection search failed! Both confidence upper bounds at tau.interval are above/below error.bound! Please rerun with a different tau.interval.")
    }
    
    #bisection search
    while(sorted.scores[tau.index.u]-sorted.scores[tau.index.l]>=tol && tau.index.u-tau.index.l>1){
        new.tau.index<-round((tau.index.l+tau.index.u)/2)
        
        Z<-matrix(as.numeric(scores<sorted.scores[new.tau.index]),ncol=1)
        Q<-CV.est.Q(A,X,Z,folds,Q.control)
        sigma2.v<-psi.v<-matrix(nrow=nfolds,ncol=1)
        for(v in 1:nfolds){
            A.train<-A[-folds[[v]]]
            A.test<-A[folds[[v]]]
            Z.test<-Z[folds[[v]],1]
            Q.test<-Q[[v]][,1]
            gamma.train<-mean(A.train)
            gamma.test<-mean(A.test)
            w.test<-w[folds[[v]]]
            
            naive.psi<-sum((1-A.test)*Q.test)/sum(1-A.test)
            IF<-ifelse(A.test,(Z.test-Q.test)*w.test/gamma.test,(Q.test-naive.psi)/(1-gamma.test))
            psi.v[v,1]<-naive.psi+mean(IF)
            sigma2.v[v,1]<-mean(IF^2)
        }
        
        new.psi<-as.numeric(t(fold.sizes)%*%psi.v/n)
        new.psi<-pmin(pmax(new.psi,0),1) #project psi to [0,1]
        new.sigma2<-as.numeric(t(fold.sizes)%*%sigma2.v/n)
        new.SE<-sqrt(new.sigma2/n)
        new.CI.upper<-new.psi+new.SE*qnorm(conf.level)
        new.sign<-sign(new.CI.upper-error.bound)
        
        if(new.sign==sign.l){
            tau.index.l<-new.tau.index
            psi.l<-new.psi
            CI.upper.l<-new.CI.upper
            sign.l<-new.sign
        }else{
            tau.index.u<-new.tau.index
            psi.u<-new.psi
            CI.upper.u<-new.CI.upper
            sign.u<-new.sign
        }
    }
    
    if(sign.u<0){
        list(tau=sorted.scores[tau.index.u],error.CI.upper=CI.upper.u,error.est=psi.u)
    }else{
        list(tau=sorted.scores[tau.index.l],error.CI.upper=CI.upper.l,error.est=psi.l)
    }
}
