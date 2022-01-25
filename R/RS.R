#' @title Confidence upper bound of true coverage error or threshold selection based on Rejection Sampling & Binomial Proportion confidence upper bound
#' @name RS
#' @export
#' @description Method to compute a confidence upper bound of true coverage or select a threshold for APAC prediction sets based on Rejection Sampling & Binomial Proportion confidence upper bound
#' @param A vector of population indicator. 1 for source population, 0 for target population
#' @param X data frame of covariates with each row being one observation
#' @param Y vector of dependent variable/outcome. For data from the target population (`A=0`), set the corresponding entries of `Y` to be `NA`
#' @param scores either a function assigning scores of `Y` given `X` trained using an independent dataset from source population or a vector of this function evaluated at observed `(X,Y)`, taking `NA` for observations from the target population. If it is a function, it must take input `(x,y)`, where `x` is one row of `X` (a data frame with one row) and `y` is a nonmissing value of `Y`, and output a scalar
#' @param candidate.tau a numeric vector of candidate thresholds, default to `c(scores,Inf)` (after `scores` is evaluated at observations if `scores` is a function). If `candidate.tau` has length 1, then just compute the point estimate and confidence upper bound of true coverage error of this threshold.
#' @param LR.bound known upper bound on likelihood ratio between target population and source population. As long as `LR.bound` is a valid upper bound, smaller values lead to better performance. If is `NULL`, will use an ad hoc choice, the maximum value of estimated likelihood ratio at observations in the testing data. Default to `NULL`
#' @param error.bound desired bound on the prediction set coverage error between 0 and 1, default 0.05
#' @param conf.level desired level of confidence of low coverage error between 0.5 and 1, default to 0.95
#' @param train.prop proportion of training data used to estimate nuisance functions, default to 0.5
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
#' When extremely small/large thresholds are included in `candidata.tau`, it is common to receive warnings/errors from the machine learning algorithms used by \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}, because in such cases, almost all `Y` are included in (for small thresholds) or excluded from (for large thresholds) the corresponding prediction sets, leading to complaints from machine learning algorithms. This is usually not an issue because the resulting predictions are still quite accurate.
#' 
#' @examples
#' n<-100
#' expit<-function(x) 1/(1+exp(-x))
#' A<-rbinom(n,1,.5)
#' X<-data.frame(X=rnorm(n,sd=ifelse(A==1,1,.5)))
#' Y<-rbinom(n,1,expit(1+X$X))
#' scores<-dbinom(Y,1,expit(.08+1.1*X$X))
#' candidate.tau<-seq(0,.5,length.out=10)
#' LR.bound<-4
#' RS(A,X,Y,scores,candidate.tau,LR.bound,
#'    g.control=list(SL.library="SL.glm"),
#'    Q.control=list(SL.library="SL.glm"))
RS<-function(A,X,Y,scores,candidate.tau,LR.bound=NULL,error.bound=0.05,conf.level=.95,train.prop=0.5,g.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),Q.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),g.trunc=1e-2,select.tau=ifelse(length(candidate.tau)==1,FALSE,TRUE)){
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
    
    assert_that(is.number(train.prop),train.prop>0,train.prop<1)
    assert_that(is.list(g.control))
    assert_that(is.list(Q.control))
    if(any(c("Y","X","newX") %in% names(g.control))){
        stop("Y, X and newX must not be specified in g.control")
    }
    if(any(c("Y","X","newX") %in% names(Q.control))){
        stop("Y, X and newX must not be specified in Q.control")
    }
    if(!is.null(LR.bound)){
        assert_that(is.number(LR.bound),LR.bound>0,LR.bound<Inf)
    }else{
        message("LR.bound is NULL. An ad hoc LR.bound will be selected.")
    }
    assert_that(is.number(g.trunc),g.trunc>0,g.trunc<1)
    assert_that(is.flag(select.tau),!is.na(select.tau))
    
    folds<-data.split(A,train.prop)
    fold.sizes<-sapply(folds,length)
    
    #nuisance function estimators
    Z<-do.call(cbind,lapply(candidate.tau,function(tau){
        as.numeric(scores<tau)
    }))
    A.train<-A[folds$train]
    A.test<-A[folds$test]
    Z.train<-Z[folds$train,,drop=FALSE]
    Z.test<-Z[folds$test,,drop=FALSE]
    X.train<-X[folds$train,,drop=FALSE]
    X.test<-X[folds$test,,drop=FALSE]
    gamma<-mean(A.train)
    
    if(all(A.train==0)){
        g<-rep(0,fold.sizes["test"])
    }else if(all(A.train==1)){
        g<-rep(1,fold.sizes["test"])
    }else{
        args<-c(list(Y=A.train,X=X.train,newX=X.test,family=binomial()),g.control)
        SL.model<-do.call(SuperLearner,args)
        g<-pmax(as.numeric(predict(SL.model)$pred),g.trunc)
    }
    
    Q<-matrix(nrow=fold.sizes["test"],ncol=length(candidate.tau))
    for(tau.index in 1:length(candidate.tau)){
        if(all(Z.train[A.train==1,tau.index]==0)){
            Q[,tau.index]<-rep(0,fold.sizes["test"])
        }else if(all(Z.train[A.train==1,tau.index]==1)){
            Q[,tau.index]<-rep(1,fold.sizes["test"])
        }else{
            args<-c(list(Y=Z.train[A.train==1,tau.index],X=X.train[A.train==1,,drop=FALSE],newX=X.test,family=binomial()),Q.control)
            SL.model<-do.call(SuperLearner,args)
            Q[,tau.index]<-predict(SL.model)$pred
        }
    }
    
    
    #rejection sampling
    w<-(1-g)/g*gamma/(1-gamma)
    w<-ifelse(is.na(w) | is.infinite(w),Inf,w)
    if(is.null(LR.bound)){
        LR.bound<-max(w)
        assert_that(!is.na(LR.bound),is.finite(LR.bound))
    }else if(any(w>LR.bound)){
        if(any(is.infinite(w))){
            stop("Estimated likelihood ratio being Inf! Check whether target population is dominated by source population!")
        }else{
            warning(paste0("Estimated likelihood ratio exceeding specified LR.bound!\nThe specified LR.bound is ",LR.bound," but the max estimated LR is ",max(w),".\nTruncating estimated likelihood ratio.\nBetter to rerun with a larger bound."))
            w<-pmin(w,LR.bound)
        }
    }
    accepted<-A.test==1 & runif(fold.sizes["test"])<=w/LR.bound
    
    #calculate one part of IF
    IF1<-ifelse(A.test==1,-w/gamma/mean(w[A.test==1]),1/(1-gamma))
    
    #calculate point estimate and confidence upper bound
    n.trial<-sum(accepted)
    if(n.trial==0){
        warning("No sample accepted in rejection smapling. Setting confidence upper bound to 1. Try rerunning with a smaller LR.bound.")
        psi<-rep(NA,length(candidate.tau))
        CI.upper<-rep(1,length(candidate.tau))
    }else{
        psi<-CI.upper<-numeric(length(candidate.tau))
        for(tau.index in 1:length(candidate.tau)){
            n.suc<-sum(Z.test[accepted,tau.index])
            
            naive.psi<-n.suc/n.trial
            IF.to.correct<-Q[,tau.index]*IF1
            psi[tau.index]<-naive.psi+mean(IF.to.correct)
            
            IF.train<-(A.train-gamma)/gamma/(1-gamma)*psi[tau.index]
            IF.test<-(ifelse(accepted,LR.bound*(Z.test[,tau.index]-psi[tau.index]),0)+A.test*(w-1)*psi[tau.index])/gamma+IF.to.correct
            SE<-sqrt(mean(IF.train^2)/fold.sizes["train"]+mean(IF.test^2)/fold.sizes["test"])
            SE<-as.numeric(SE)
            CI.upper[tau.index]<-psi[tau.index]+SE*qnorm(conf.level)
        }
        psi<-pmin(pmax(psi,0),1) #project psi to [0,1]
        CI.upper<-pmin(pmax(CI.upper,0),1) #project CI.upper to [0,1]
    }
    
    #select tau if needed
    if(!select.tau){
        list(tau=candidate.tau,error.CI.upper=as.numeric(CI.upper),error.est=as.numeric(psi))
    }else if(any(CI.upper<=error.bound)){
        feasible.indicator<-CI.upper<=error.bound
        if(all(feasible.indicator)){
            warning("Larest threshold is selected. Try a different set of candidate.tau with a larger max value.")
        }
        tau<-max(candidate.tau[feasible.indicator])
        tau.index<-which(candidate.tau==tau)
        list(tau=tau,error.CI.upper=CI.upper[tau.index],error.est=psi[tau.index],
             feasible.tau=candidate.tau[feasible.indicator],feasible.tau.error.CI.upper=CI.upper[feasible.indicator],feasible.tau.error.est=psi[feasible.indicator])
    }else{
        message("No candidate threshold tau satisfies the criterion! Try smaller candidate.tau!")
        invisible()
    }
}

















# Confidence upper bound of true coverage error or threshold selection based on Rejection Sampling & Binomial Proportion confidence upper bound with known likelihood ratio
# Method to compute a confidence upper bound of true coverage or select a threshold for APAC prediction sets based on Rejection Sampling & Binomial Proportion confidence upper bound with known likelihood ratio
# A: vector of population indicator. 1 for source population, 0 for target population
# X: data frame of covariates with each row being one observation
# Y: vector of dependent variable/outcome. For data from the target population (`A=0`), set the corresponding entries of `Y` to be `NA`
# scores: either a function assigning scores of `Y` given `X` trained using an independent dataset from source population or a vector of this function evaluated at observed `(X,Y)`, taking `NA` for observations from the target population. If it is a function, it must take input `(x,y)`, where `x` is one row of `X` (a data frame with one row) and `y` is a nonmissing value of `Y`, and output a scalar
# candidate.tau: a numeric vector of candidate thresholds, default to `c(scores,Inf)` (after `scores` is evaluated at observations if `scores` is a function). If `candidate.tau` has length 1, then just compute the point estimate and confidence upper bound of true coverage error of this threshold.
# LR: known likelihood ratio function (ratio of the density of covariate `X` from target population to the density of covariate `X` from source population). Must be a function that takes in a row of `X` and outputs a non-negative number. Default to constant function 1, i.e., no covariate shift
# LR.bound: known upper bound on likelihood ratio between target population and source population. As long as `LR.bound` is a valid upper bound, smaller values lead to better performance. If is `NULL`, will use an ad hoc choice, the maximum value of likelihood ratio at all observations. Default to `NULL`
# error.bound: desired bound on the prediction set coverage error between 0 and 1, default 0.05
# conf.level: desired level of confidence of low coverage error between 0.5 and 1, default to 0.95
# CI.method: a string indicating the method to compute binomial confidence upper bound that is passed to \code{\link[DescTools:BinomCI]{DescTools::BinomCI}}, default to `"wilson"` (Wilson score). Alternatively, `CI.method` may be a function that takes in three arguments in order (number of successes, number of trials, confidence level) and outputs a vector of with two numbers (point estimate, confidence upper bound)
# select.tau: whether to select threshold tau (otherwise just reposrt estimates and confidence upper bounds of coverage error for all `candidate.tau`), default to `TRUE` if `length(candidate.tau)>1` and `FALSE` if `length(candidate.tau)==1`
# 
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
# example:
# n<-100
# expit<-function(x) 1/(1+exp(-x))
# A<-rbinom(n,1,.5)
# X<-data.frame(X=rnorm(n,sd=ifelse(A==1,1,.5)))
# Y<-rbinom(n,1,expit(1+X$X))
# scores<-dbinom(Y,1,expit(.08+1.1*X$X))
# LR<-function(x) dnorm(x$X,sd=.5)/dnorm(x$X)
# candidate.tau<-seq(0,.5,length.out=10)
# LR.bound<-2
# APACpred:::RS_knownLR(A,X,Y,scores,candidate.tau,LR,LR.bound,g.control=list(SL.library="SL.glm"),Q.control=list(SL.library="SL.glm"))
RS_knownLR<-function(A,X,Y,scores,candidate.tau,LR=function(x) 1,LR.bound=NULL,error.bound=0.05,conf.level=.95,CI.method="wilson",select.tau=ifelse(length(candidate.tau)==1,FALSE,TRUE)){
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
    
    if(!is.null(LR.bound)){
        assert_that(is.number(LR.bound),LR.bound>0,LR.bound<Inf)
    }else{
        message("LR.bound is NULL. An ad hoc LR.bound will be selected.")
    }
    assert_that(is.string(CI.method) || is.function(CI.method))
    assert_that(is.flag(select.tau),!is.na(select.tau))
    
    #nuisance function estimators
    Z<-do.call(cbind,lapply(candidate.tau,function(tau){
        as.numeric(scores<tau)
    }))
    w<-sapply(1:nrow(X),function(i) LR(X[i,,drop=FALSE]))
    if(is.null(LR.bound)){
        LR.bound<-max(w)
        assert_that(!is.na(LR.bound),is.finite(LR.bound))
    }else{
        assert_that(is.numeric(w),all(w>=0 & w<=LR.bound))
    }
    gamma<-mean(A)
    
    #rejection sampling
    accepted<-A==1 & runif(n)<=w/LR.bound
    
    #calculate point estimate and confidence upper bound
    n.trial<-sum(accepted)
    if(n.trial==0){
        warning("No sample accepted in rejection smapling. Setting confidence upper bound to 1. Try rerunning with a smaller LR.bound.")
        psi<-rep(NA,length(candidate.tau))
        CI.upper<-rep(1,length(candidate.tau))
    }else{
        psi<-CI.upper<-numeric(length(candidate.tau))
        for(tau.index in 1:length(candidate.tau)){
            n.suc<-sum(Z[accepted,tau.index])
            
            if(is.function(CI.method)){
                dummy<-CI.method(n.suc,n.trial,conf.level)
                psi[tau.index]<-dummy[1]
                CI.upper[tau.index]<-dummy[2]
            }else{
                dummy<-BinomCI(n.suc,n.trial,conf.level,sides="right",method=CI.method)
                psi[tau.index]<-dummy[1]
                CI.upper[tau.index]<-dummy[3]
            }
        }
        psi<-pmin(pmax(psi,0),1) #project psi to [0,1]
        CI.upper<-pmin(pmax(CI.upper,0),1) #project CI.upper to [0,1]
    }
    
    #select tau if needed
    if(!select.tau){
        list(tau=candidate.tau,error.CI.upper=as.numeric(CI.upper),error.est=as.numeric(psi))
    }else if(any(CI.upper<=error.bound)){
        feasible.indicator<-CI.upper<=error.bound
        if(all(feasible.indicator)){
            warning("Larest threshold is selected. Try a different set of candidate.tau with a larger max value.")
        }
        tau<-max(candidate.tau[feasible.indicator])
        tau.index<-which(candidate.tau==tau)
        list(tau=tau,error.CI.upper=CI.upper[tau.index],error.est=psi[tau.index],
             feasible.tau=candidate.tau[feasible.indicator],feasible.tau.error.CI.upper=CI.upper[feasible.indicator],feasible.tau.error.est=psi[feasible.indicator])
    }else{
        message("No candidate threshold tau satisfies the criterion! Try smaller candidate.tau!")
        invisible()
    }
}










#' @title Confidence upper bound of true coverage error or threshold selection based on Rejection Sampling & Binomial Proportion confidence upper bound via bisection search
#' @name RS_bisec
#' @export
#' @description Method to compute a confidence upper bound of true coverage or select a threshold for APAC prediction sets based on Rejection Sampling & Binomial Proportion confidence upper bound via bisection search
#' @param A vector of population indicator. 1 for source population, 0 for target population
#' @param X data frame of covariates with each row being one observation
#' @param Y vector of dependent variable/outcome. For data from the target population (`A=0`), set the corresponding entries of `Y` to be `NA`
#' @param scores either a function assigning scores of `Y` given `X` trained using an independent dataset from source population or a vector of this function evaluated at observed `(X,Y)`, taking `NA` for observations from the target population. If it is a function, it must take input `(x,y)`, where `x` is one row of `X` (a data frame with one row) and `y` is a nonmissing value of `Y`, and output a scalar
#' @param tau.interval a numeric vector specifying the initial interval of candidate thresholds for bisection search, default to `c(0,1)`, which is suitable when all scores lie in the interval \eqn{[0,1]}{[0,1]}
#' @param LR.bound known upper bound on likelihood ratio between target population and source population. As long as `LR.bound` is a valid upper bound, smaller values lead to better performance. If is `NULL`, will use an ad hoc choice, the maximum value of estimated likelihood ratio at observations in the testing data. Default to `NULL`
#' @param error.bound desired bound on the prediction set coverage error between 0 and 1, default 0.05
#' @param conf.level desired level of confidence of low coverage error between 0.5 and 1, default to 0.95
#' @param train.prop proportion of training data used to estimate nuisance functions, default to 0.5
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
#' When extremely small/large thresholds are included in `candidata.tau`, it is common to receive warnings/errors from the machine learning algorithms used by \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}, because in such cases, almost all `Y` are included in (for small thresholds) or excluded from (for large thresholds) the corresponding prediction sets, leading to complaints from machine learning algorithms. This is usually not an issue because the resulting predictions are still quite accurate.
#' 
#' @examples
#' n<-100
#' expit<-function(x) 1/(1+exp(-x))
#' A<-rbinom(n,1,.5)
#' X<-data.frame(X=rnorm(n,sd=ifelse(A==1,1,.5)))
#' Y<-rbinom(n,1,expit(1+X$X))
#' scores<-dbinom(Y,1,expit(.08+1.1*X$X))
#' LR.bound<-4
#' RS_bisec(A,X,Y,scores,LR.bound=LR.bound,
#'          g.control=list(SL.library="SL.glm"),
#'          Q.control=list(SL.library="SL.glm"))
RS_bisec<-function(A,X,Y,scores,tau.interval=c(0,1),LR.bound=NULL,error.bound=0.05,conf.level=.95,train.prop=0.5,g.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),Q.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),g.trunc=1e-2,tol=1e-2){
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
    
    assert_that(is.number(train.prop),train.prop>0,train.prop<1)
    assert_that(is.list(g.control))
    assert_that(is.list(Q.control))
    if(any(c("Y","X","newX") %in% names(g.control))){
        stop("Y, X and newX must not be specified in g.control")
    }
    if(any(c("Y","X","newX") %in% names(Q.control))){
        stop("Y, X and newX must not be specified in Q.control")
    }
    if(!is.null(LR.bound)){
        assert_that(is.number(LR.bound),LR.bound>0,LR.bound<Inf)
    }else{
        message("LR.bound is NULL. An ad hoc LR.bound will be selected.")
    }
    assert_that(is.number(g.trunc),g.trunc>0,g.trunc<1)
    assert_that(is.number(tol),tol>0)
    
    folds<-data.split(A,train.prop)
    fold.sizes<-sapply(folds,length)
    
    #nuisance function estimators
    Z<-do.call(cbind,lapply(tau.interval,function(tau){
        as.numeric(scores<tau)
    }))
    A.train<-A[folds$train]
    A.test<-A[folds$test]
    Z.train<-Z[folds$train,,drop=FALSE]
    Z.test<-Z[folds$test,,drop=FALSE]
    X.train<-X[folds$train,,drop=FALSE]
    X.test<-X[folds$test,,drop=FALSE]
    gamma<-mean(A.train)
    
    if(all(A.train==0)){
        g<-rep(0,fold.sizes["test"])
    }else if(all(A.train==1)){
        g<-rep(1,fold.sizes["test"])
    }else{
        args<-c(list(Y=A.train,X=X.train,newX=X.test,family=binomial()),g.control)
        SL.model<-do.call(SuperLearner,args)
        g<-pmax(as.numeric(predict(SL.model)$pred),g.trunc)
    }
    
    Q<-matrix(nrow=fold.sizes["test"],ncol=2)
    for(tau.index in 1:2){
        if(all(Z.train[A.train==1,tau.index]==0)){
            Q[,tau.index]<-rep(0,fold.sizes["test"])
        }else if(all(Z.train[A.train==1,tau.index]==1)){
            Q[,tau.index]<-rep(1,fold.sizes["test"])
        }else{
            args<-c(list(Y=Z.train[A.train==1,tau.index],X=X.train[A.train==1,,drop=FALSE],newX=X.test,family=binomial()),Q.control)
            SL.model<-do.call(SuperLearner,args)
            Q[,tau.index]<-predict(SL.model)$pred
        }
    }
    
    
    #rejection sampling
    w<-(1-g)/g*gamma/(1-gamma)
    w<-ifelse(is.na(w) | is.infinite(w),Inf,w)
    if(is.null(LR.bound)){
        LR.bound<-max(w)
        assert_that(!is.na(LR.bound),is.finite(LR.bound))
    }else if(any(w>LR.bound)){
        if(any(is.infinite(w))){
            stop("Estimated likelihood ratio being Inf! Check whether target population is dominated by source population!")
        }else{
            warning(paste0("Estimated likelihood ratio exceeding specified LR.bound!\nThe specified LR.bound is ",LR.bound," but the max estimated LR is ",max(w),".\nTruncating estimated likelihood ratio.\nBetter to rerun with a larger bound."))
            w<-pmin(w,LR.bound)
        }
    }
    accepted<-A.test==1 & runif(fold.sizes["test"])<=w/LR.bound
    
    #calculate one part of IF
    IF1<-ifelse(A.test==1,-w/gamma/mean(w[A.test==1]),1/(1-gamma))
    
    #calculate point estimate and confidence upper bound at end points of tau.interval
    n.trial<-sum(accepted)
    if(n.trial==0){
        warning("No sample accepted in rejection smapling. Setting confidence upper bound to 1. Try rerunning with a smaller LR.bound.")
        psi<-rep(NA,2)
        CI.upper<-rep(1,2)
    }else{
        psi<-CI.upper<-numeric(2)
        for(tau.index in 1:2){
            n.suc<-sum(Z.test[accepted,tau.index])
            
            naive.psi<-n.suc/n.trial
            IF.to.correct<-Q[,tau.index]*IF1
            psi[tau.index]<-naive.psi+mean(IF.to.correct)
            
            IF.train<-(A.train-gamma)/gamma/(1-gamma)*psi[tau.index]
            IF.test<-(ifelse(accepted,LR.bound*(Z.test[,tau.index]-psi[tau.index]),0)+A.test*(w-1)*psi[tau.index])/gamma+IF.to.correct
            SE<-sqrt(mean(IF.train^2)/fold.sizes["train"]+mean(IF.test^2)/fold.sizes["test"])
            SE<-as.numeric(SE)
            CI.upper[tau.index]<-psi[tau.index]+SE*qnorm(conf.level)
        }
        psi<-pmin(pmax(psi,0),1) #project psi to [0,1]
        CI.upper<-pmin(pmax(CI.upper,0),1) #project CI.upper to [0,1]
    }
    
    if(any(CI.upper==error.bound)){
        tau<-max(tau.interval[CI.upper==error.bound])
        tau.index<-which(tau.interval==tau)
        return(list(tau=tau,error.CI.upper=CI.upper[tau.index],error.est=psi[tau.index]))
    }
    
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
        Z.train<-Z[folds$train,,drop=FALSE]
        Z.test<-Z[folds$test,,drop=FALSE]
        
        Q<-matrix(nrow=fold.sizes["test"],ncol=1)
        if(all(Z.train[A.train==1,1]==0)){
            Q[,1]<-rep(0,fold.sizes["test"])
        }else if(all(Z.train[A.train==1,1]==1)){
            Q[,1]<-rep(1,fold.sizes["test"])
        }else{
            args<-c(list(Y=Z.train[A.train==1,1],X=X.train[A.train==1,,drop=FALSE],newX=X.test,family=binomial()),Q.control)
            SL.model<-do.call(SuperLearner,args)
            Q[,1]<-predict(SL.model)$pred
        }
        
        n.suc<-sum(Z.test[accepted,1])
        
        naive.psi<-n.suc/n.trial
        IF.to.correct<-Q[,1]*IF1
        new.psi<-naive.psi+mean(IF.to.correct)
        
        IF.train<-(A.train-gamma)/gamma/(1-gamma)*new.psi
        IF.test<-(ifelse(accepted,LR.bound*(Z.test-new.psi),0)+A.test*(w-1)*new.psi)/gamma+IF.to.correct
        SE<-sqrt(mean(IF.train^2)/fold.sizes["train"]+mean(IF.test^2)/fold.sizes["test"])
        SE<-as.numeric(SE)
        new.CI.upper<-new.psi+SE*qnorm(conf.level)
        
        new.sign<-sign(new.CI.upper-error.bound)
        
        if(new.sign==0){
            return(list(tau=new.tau,error.CI.upper=new.CI.upper,error.est=new.psi))
        }else if(new.sign==sign.l){
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










# Confidence upper bound of true coverage error or threshold selection based on Rejection Sampling & Binomial Proportion confidence upper bound via bisection search with known likelihood ratio
# Method to compute a confidence upper bound of true coverage or select a threshold for APAC prediction sets based on Rejection Sampling & Binomial Proportion confidence upper bound via bisection search with known likelihood ratio
# A: vector of population indicator. 1 for source population, 0 for target population
# X: data frame of covariates with each row being one observation
# Y: vector of dependent variable/outcome. For data from the target population (`A=0`), set the corresponding entries of `Y` to be `NA`
# scores: either a function assigning scores of `Y` given `X` trained using an independent dataset from source population or a vector of this function evaluated at observed `(X,Y)`, taking `NA` for observations from the target population. If it is a function, it must take input `(x,y)`, where `x` is one row of `X` (a data frame with one row) and `y` is a nonmissing value of `Y`, and output a scalar
# tau.interval: a numeric vector specifying the initial interval of candidate thresholds for bisection search, default to `c(0,1)`, which is suitable when all scores lie in the interval \eqn{[0,1]}{[0,1]}
# LR: known likelihood ratio function (ratio of the density of covariate `X` from target population to the density of covariate `X` from source population). Must be a function that takes in a row of `X` and outputs a non-negative number. Default to constant function 1, i.e., no covariate shift
# LR.bound: known upper bound on likelihood ratio between target population and source population. As long as `LR.bound` is a valid upper bound, smaller values lead to better performance. If is `NULL`, will use an ad hoc choice, the maximum value of estimated likelihood ratio at al observations. Default to `NULL`
# error.bound: desired bound on the prediction set coverage error between 0 and 1, default 0.05
# conf.level: desired level of confidence of low coverage error between 0.5 and 1, default to 0.95
# CI.method: a string indicating the method to compute binomial confidence upper bound that is passed to \code{\link[DescTools:BinomCI]{DescTools::BinomCI}}, default to `"wilson"` (Wilson score). Alternatively, `CI.method` may be a function that takes in three arguments in order (number of successes, number of trials, confidence level) and outputs a vector of with two numbers (point estimate, confidence upper bound)
# tol: stopping criterion for bisection search, default to 0.01. The search will stop if the length of the interval for thresholds is below `tol` and the function will output the end point with confidence upper bound no greater than `error.bound`
# 
# Output:
# a list with the following components:
# \describe{
# \item{`tau`}{Selected tau}
# \item{`error.CI.upper`}{The (approximate) confidence upper bound of coverage error corresponding to the select tau}
# \item{`error.est`}{The point estimate of coverage error corresponding to the selected tau}
# }
# 
# example:
# n<-100
# expit<-function(x) 1/(1+exp(-x))
# A<-rbinom(n,1,.5)
# X<-data.frame(X=rnorm(n,sd=ifelse(A==1,1,.5)))
# Y<-rbinom(n,1,expit(1+X$X))
# scores<-dbinom(Y,1,expit(.08+1.1*X$X))
# LR<-function(x) dnorm(x$X,sd=.5)/dnorm(x$X)
# LR.bound<-2
# APACpred:::RS_knownLR_bisec(A,X,Y,scores,LR=LR,LR.bound=LR.bound,g.control=list(SL.library="SL.glm"),Q.control=list(SL.library="SL.glm"))
RS_knownLR_bisec<-function(A,X,Y,scores,tau.interval=c(0,1),LR=function(x) 1,LR.bound=NULL,error.bound=0.05,conf.level=.95,CI.method="wilson",tol=1e-2){
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
    
    if(!is.null(LR.bound)){
        assert_that(is.number(LR.bound),LR.bound>0,LR.bound<Inf)
    }else{
        message("LR.bound is NULL. An ad hoc LR.bound will be selected.")
    }
    assert_that(is.string(CI.method) || is.function(CI.method))
    assert_that(is.number(tol),tol>0)
    
    #nuisance function estimators
    Z<-do.call(cbind,lapply(tau.interval,function(tau){
        as.numeric(scores<tau)
    }))
    w<-sapply(1:nrow(X),function(i) LR(X[i,,drop=FALSE]))
    if(is.null(LR.bound)){
        LR.bound<-max(w)
        assert_that(is.number(LR.bound),LR.bound>0,LR.bound<Inf)
    }else{
        assert_that(is.numeric(w),all(w>=0 & w<=LR.bound))
    }
    
    #rejection sampling
    accepted<-A==1 & runif(n)<=w/LR.bound
    
    #calculate point estimate and confidence upper bound at the end points of tau.interval
    n.trial<-sum(accepted)
    if(n.trial==0){
        warning("No sample accepted in rejection smapling. Setting confidence upper bound to 1. Try rerunning with a smaller LR.bound.")
        psi<-rep(NA,2)
        CI.upper<-rep(1,2)
    }else{
        psi<-CI.upper<-numeric(2)
        for(tau.index in 1:2){
            n.suc<-sum(Z[accepted,tau.index])
            
            if(is.function(CI.method)){
                dummy<-CI.method(n.suc,n.trial,conf.level)
                psi[tau.index]<-dummy[1]
                CI.upper[tau.index]<-dummy[2]
            }else{
                dummy<-BinomCI(n.suc,n.trial,conf.level,sides="right",method=CI.method)
                psi[tau.index]<-dummy[1]
                CI.upper[tau.index]<-dummy[3]
            }
        }
        psi<-pmin(pmax(psi,0),1) #project psi to [0,1]
        CI.upper<-pmin(pmax(CI.upper,0),1) #project CI.upper to [0,1]
    }
    
    if(any(CI.upper==error.bound)){
        tau<-max(tau.interval[CI.upper==error.bound])
        tau.index<-which(tau.interval==tau)
        return(list(tau=tau,error.CI.upper=CI.upper[tau.index],error.est=psi[tau.index]))
    }
    
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
        
        n.suc<-sum(Z[accepted,1])
        
        if(is.function(CI.method)){
            dummy<-CI.method(n.suc,n.trial,conf.level)
            new.psi<-dummy[1]
            new.CI.upper<-dummy[2]
        }else{
            dummy<-BinomCI(n.suc,n.trial,conf.level,sides="right",method=CI.method)
            new.psi<-dummy[1]
            new.CI.upper<-dummy[3]
        }
        new.psi<-pmin(pmax(new.psi,0),1) #project psi to [0,1]
        new.CI.upper<-pmin(pmax(new.CI.upper,0),1) #project CI.upper to [0,1]
        new.sign<-sign(new.CI.upper-error.bound)
        
        if(new.sign==0){
            return(list(tau=new.tau,error.CI.upper=new.CI.upper,error.est=new.psi))
        }else if(new.sign==sign.l){
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





















# Similar to RS_bisec, but does bisection search over observed scores
# 
# example:
# n<-100
# expit<-function(x) 1/(1+exp(-x))
# A<-rbinom(n,1,.5)
# X<-data.frame(X=rnorm(n,sd=ifelse(A==1,1,.5)))
# Y<-rbinom(n,1,expit(1+X$X))
# scores<-dbinom(Y,1,expit(.08+1.1*X$X))
# LR.bound<-4
# APACpred:::RS_bisec2(A,X,Y,scores,LR.bound=LR.bound,g.control=list(SL.library="SL.glm"),Q.control=list(SL.library="SL.glm"))
RS_bisec2<-function(A,X,Y,scores,LR.bound=NULL,error.bound=0.05,conf.level=.95,train.prop=0.5,g.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),Q.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),g.trunc=1e-2,tol=1e-2){
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
    
    assert_that(is.number(train.prop),train.prop>0,train.prop<1)
    assert_that(is.list(g.control))
    assert_that(is.list(Q.control))
    if(any(c("Y","X","newX") %in% names(g.control))){
        stop("Y, X and newX must not be specified in g.control")
    }
    if(any(c("Y","X","newX") %in% names(Q.control))){
        stop("Y, X and newX must not be specified in Q.control")
    }
    if(!is.null(LR.bound)){
        assert_that(is.number(LR.bound),LR.bound>0,LR.bound<Inf)
    }else{
        message("LR.bound is NULL. An ad hoc LR.bound will be selected.")
    }
    assert_that(is.number(g.trunc),g.trunc>0,g.trunc<1)=
    assert_that(is.number(tol),tol>0)
    
    folds<-data.split(A,train.prop)
    fold.sizes<-sapply(folds,length)
    
    sorted.scores<-c(sort(scores),Inf)
    tau.index.interval<-c(1,length(sorted.scores))
    
    #nuisance function estimators
    Z<-do.call(cbind,lapply(tau.index.interval,function(tau.index){
        as.numeric(scores<sorted.scores[tau.index])
    }))
    A.train<-A[folds$train]
    A.test<-A[folds$test]
    Z.train<-Z[folds$train,,drop=FALSE]
    Z.test<-Z[folds$test,,drop=FALSE]
    X.train<-X[folds$train,,drop=FALSE]
    X.test<-X[folds$test,,drop=FALSE]
    gamma<-mean(A.train)
    
    if(all(A.train==0)){
        g<-rep(0,fold.sizes["test"])
    }else if(all(A.train==1)){
        g<-rep(1,fold.sizes["test"])
    }else{
        args<-c(list(Y=A.train,X=X.train,newX=X.test,family=binomial()),g.control)
        SL.model<-do.call(SuperLearner,args)
        g<-pmax(as.numeric(predict(SL.model)$pred),g.trunc)
    }
    
    Q<-matrix(nrow=fold.sizes["test"],ncol=2)
    for(tau.index in 1:2){
        if(all(Z.train[A.train==1,tau.index]==0)){
            Q[,tau.index]<-rep(0,fold.sizes["test"])
        }else if(all(Z.train[A.train==1,tau.index]==1)){
            Q[,tau.index]<-rep(1,fold.sizes["test"])
        }else{
            args<-c(list(Y=Z.train[A.train==1,tau.index],X=X.train[A.train==1,,drop=FALSE],newX=X.test,family=binomial()),Q.control)
            SL.model<-do.call(SuperLearner,args)
            Q[,tau.index]<-predict(SL.model)$pred
        }
    }
    
    
    #rejection sampling
    w<-(1-g)/g*gamma/(1-gamma)
    w<-ifelse(is.na(w) | is.infinite(w),Inf,w)
    if(is.null(LR.bound)){
        LR.bound<-max(w)
        assert_that(!is.na(LR.bound),is.finite(LR.bound))
    }else if(any(w>LR.bound)){
        if(any(is.infinite(w))){
            stop("Estimated likelihood ratio being Inf! Check whether target population is dominated by source population!")
        }else{
            warning(paste0("Estimated likelihood ratio exceeding specified LR.bound!\nThe specified LR.bound is ",LR.bound," but the max estimated LR is ",max(w),".\nTruncating estimated likelihood ratio.\nBetter to rerun with a larger bound."))
            w<-pmin(w,LR.bound)
        }
    }
    accepted<-A.test==1 & runif(fold.sizes["test"])<=w/LR.bound
    
    #calculate one part of IF
    IF1<-ifelse(A.test==1,-w/gamma/mean(w[A.test==1]),1/(1-gamma))
    
    #calculate point estimate and confidence upper bound at the end points of tau.interval
    n.trial<-sum(accepted)
    if(n.trial==0){
        warning("No sample accepted in rejection smapling. Setting confidence upper bound to 1. Try rerunning with a smaller LR.bound.")
        psi<-rep(NA,2)
        CI.upper<-rep(1,2)
    }else{
        psi<-CI.upper<-numeric(2)
        for(tau.index in 1:2){
            n.suc<-sum(Z.test[accepted,tau.index])
            
            naive.psi<-n.suc/n.trial
            IF.to.correct<-Q[,tau.index]*IF1
            psi[tau.index]<-naive.psi+mean(IF.to.correct)
            
            IF.train<-(A.train-gamma)/gamma/(1-gamma)*psi[tau.index]
            IF.test<-(ifelse(accepted,LR.bound*(Z.test[,tau.index]-psi[tau.index]),0)+A.test*(w-1)*psi[tau.index])/gamma+IF.to.correct
            SE<-sqrt(mean(IF.train^2)/fold.sizes["train"]+mean(IF.test^2)/fold.sizes["test"])
            SE<-as.numeric(SE)
            CI.upper[tau.index]<-psi[tau.index]+SE*qnorm(conf.level)
        }
        psi<-pmin(pmax(psi,0),1) #project psi to [0,1]
        CI.upper<-pmin(pmax(CI.upper,0),1) #project CI.upper to [0,1]
    }
    
    if(any(CI.upper==error.bound)){
        tau<-max(sorted.scores[CI.upper==error.bound])
        tau.index<-which(sorted.scores[tau.index.interval]==tau)
        return(list(tau=tau,error.CI.upper=CI.upper[tau.index],error.est=psi[tau.index]))
    }
    
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
        Z.train<-Z[folds$train,,drop=FALSE]
        Z.test<-Z[folds$test,,drop=FALSE]
        
        Q<-matrix(nrow=fold.sizes["test"],ncol=1)
        if(all(Z.train[A.train==1,1]==0)){
            Q[,1]<-rep(0,fold.sizes["test"])
        }else if(all(Z.train[A.train==1,1]==1)){
            Q[,1]<-rep(1,fold.sizes["test"])
        }else{
            args<-c(list(Y=Z.train[A.train==1,1],X=X.train[A.train==1,,drop=FALSE],newX=X.test,family=binomial()),Q.control)
            SL.model<-do.call(SuperLearner,args)
            Q[,1]<-predict(SL.model)$pred
        }
        
        n.suc<-sum(Z.test[accepted,1])
        
        naive.psi<-n.suc/n.trial
        IF.to.correct<-Q*IF1
        new.psi<-naive.psi+mean(IF.to.correct)
        
        IF.train<-(A.train-gamma)/gamma/(1-gamma)*new.psi
        IF.test<-(ifelse(accepted,LR.bound*(Z.test-new.psi),0)+A.test*(w-1)*new.psi)/gamma+IF.to.correct
        SE<-sqrt(mean(IF.train^2)/fold.sizes["train"]+mean(IF.test^2)/fold.sizes["test"])
        SE<-as.numeric(SE)
        new.CI.upper<-new.psi+SE*qnorm(conf.level)
        
        new.sign<-sign(new.CI.upper-error.bound)
        
        if(new.sign==0){
            return(list(tau=sorted.scores[new.tau.index],error.CI.upper=new.CI.upper,error.est=new.psi))
        }else if(new.sign==sign.l){
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
