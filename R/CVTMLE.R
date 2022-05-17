#' @title Confidence upper bound of true coverage error or threshold selection based on cross-validated targeted minimum loss-based estimator (via grid search)
#' @name CVTMLE
#' @export
#' @description Method to compute a confidence upper bound of true coverage or select a threshold (via grid search) for APAC prediction sets based on cross-validated targeted minimum loss-based estimators
#' @param A vector of population indicator. 1 for source population, 0 for target population
#' @param X data frame of covariates with each row being one observation
#' @param Y vector of dependent variable/outcome. For data from the target population (`A=0`), set the corresponding entries of `Y` to be `NA`
#' @param scores either a function assigning scores of `Y` given `X` trained using an independent dataset from source population or a vector of this function evaluated at observed `(X,Y)`, taking `NA` for observations from the target population. If it is a function, it must take input `(x,y)`, where `x` is one row of `X` (a data frame with one row) and `y` is a nonmissing value of `Y`, and output a scalar
#' @param candidate.tau a numeric vector of candidate thresholds, default to `c(scores,Inf)` (after `scores` is evaluated at observations if `scores` is a function). If `candidate.tau` has length 1, then just compute the point estimate and confidence upper bound of true coverage error of this threshold.
#' @param error.bound desired bound on the prediction set coverage error between 0 and 1, default 0.05
#' @param conf.level desired level of confidence of low coverage error between 0.5 and 1, default to 0.95
#' @param nfolds number of folds for sample splitting, default to 5
#' @param g.control a named list containing options passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}} to estimate propensity score `g`. Must not specify `Y`, `X`, `newX` or `family`. Default to `list(SL.library=c("SL.glm","SL.gam","SL.randomForest"))`
#' @param Q.control a named list containing options passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}} to estimate conditional coverage error `Q`. Must not specify `Y`, `X`, `newX` or `family`. Default to `list(SL.library=c("SL.glm","SL.gam","SL.randomForest"))`
#' @param g.trunc truncation level of initially estimated propensity score `g` from zero, default to 0.01
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
#' CVTMLE(A,X,Y,scores,candidate.tau,nfolds=2,
#'        g.control=list(SL.library="SL.glm"),
#'        Q.control=list(SL.library="SL.glm"))
CVTMLE<-function(A,X,Y,scores,candidate.tau,error.bound=0.05,conf.level=.95,nfolds=5,g.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),Q.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),g.trunc=1e-2,select.tau=ifelse(length(candidate.tau)==1,FALSE,TRUE)){
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
    
    g<-CV.est.g(A,X,folds,g.control,g.trunc)
    
    psi<-SE<-CI.upper<-numeric(length(candidate.tau))
    for(tau.index in 1:length(candidate.tau)){
        tau<-candidate.tau[tau.index]
        
        #nuisance function estimators
        Z<-cbind(as.numeric(scores<tau))
        Q<-CV.est.Q(A,X,Z,folds,Q.control)
        
        #calculate coverage error estimator and asymptotic variance for each fold and the end points of tau.interval
        sigma2.v<-psi.v<-matrix(nrow=nfolds,ncol=1)
        for(v in 1:nfolds){
            A.test<-A[folds[[v]]]
            Z.test<-Z[folds[[v]],1]
            g.test<-g[[v]]
            Q.test<-Q[[v]][,1]
            gamma.test<-mean(A.test)
            w.test<-(1-g.test)/g.test*gamma.test/(1-gamma.test)
            
            if(all(Z.test[A.test==1]==1)){
                Qhat<-rep(1,fold.sizes[v])
            }else if(all(Z.test[A.test==1]==0)){
                Qhat<-rep(0,fold.sizes[v])
            }else{
                clever.covariate<-w.test/gamma.test
                Qhat<-tryCatch({
                    if(any(Q.test==0 | Q.test==1)){
                        warning("")
                    }
                    offset<-logit(Q.test)
                    logistic.model<-glm(Z.test~clever.covariate-1,family=binomial(),subset=A.test==1,offset=offset)
                    beta<-coef(logistic.model)
                    expit(offset+beta*clever.covariate)
                },warning=function(...){
                    offset<-Q.test
                    linear.model<-lm(Z.test~clever.covariate-1,subset=A.test==1,offset=offset)
                    beta<-coef(linear.model)
                    offset+beta*clever.covariate
                })
            }
            
            psi.v[v,1]<-sum((1-A.test)*Qhat)/sum(1-A.test)
            IF<-ifelse(A.test,(Z.test-Qhat)*w.test/gamma.test,(Qhat-psi.v[v,1])/(1-gamma.test))
            sigma2.v[v,1]<-mean(IF^2)
        }
        
        psi[tau.index]<-t(fold.sizes)%*%psi.v/n
        psi[tau.index]<-pmin(pmax(psi[tau.index],0),1) #project psi to [0,1]
        sigma2<-t(fold.sizes)%*%sigma2.v/n
        SE[tau.index]<-sqrt(sigma2/n)
        CI.upper[tau.index]<-psi[tau.index]+SE[tau.index]*qnorm(conf.level)
        
        if(select.tau && CI.upper[tau.index]>=error.bound){
            tau.index<-tau.index-1
            break
        }
    }
    
    #select tau if needed
    if(!select.tau){
        list(tau=candidate.tau,error.CI.upper=as.numeric(CI.upper),error.est=as.numeric(psi))
    }else if(tau.index>0){
        if(tau.index==length(candidate.tau)){
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
