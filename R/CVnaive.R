#for simulation only
#naive cross-fit plug-in method (i.e., CVonestep without one-step correction)
#arguments see CVonestep
CVnaive<-function(A,X,Y,scores,candidate.tau,error.bound=0.05,conf.level=.95,nfolds=5,g.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),Q.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),g.trunc=1e-2,select.tau=ifelse(length(candidate.tau)==1,FALSE,TRUE)){
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
            
            naive.psi<-sum((1-A.test)*Q.test)/sum(1-A.test)
            IF<-ifelse(A.test,(Z.test-Q.test)*w.test/gamma.test,(Q.test-naive.psi)/(1-gamma.test))
            psi.v[v,1]<-naive.psi
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





#for simulation only
#naive cross-fit plug-in method based on weighting (i.e., plug in estimated LR without one-step correction)
#arguments see CVonestep
CVnaive.weight<-function(A,X,Y,scores,candidate.tau,error.bound=0.05,conf.level=.95,nfolds=5,g.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),Q.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),g.trunc=1e-2,select.tau=ifelse(length(candidate.tau)==1,FALSE,TRUE)){
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
            
            naive.psi<-sum(ifelse(A.test==1,Z.test,0))/sum(A.test)
            IF<-ifelse(A.test,(Z.test-Q.test)*w.test/gamma.test,(Q.test-naive.psi)/(1-gamma.test))
            psi.v[v,1]<-naive.psi
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


