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
            psi.v[v,tau.index]<-naive.psi
            sigma2.v[v,tau.index]<-mean(IF^2)
        }
    }
    
    #combine folds
    psi<-t(fold.sizes)%*%psi.v/n
    psi<-pmin(pmax(psi,0),1) #project psi to [0,1]
    sigma2<-t(fold.sizes)%*%sigma2.v/n
    SE<-sqrt(sigma2/n)
    
    #calculate confidence upper limit and select tau if needed
    CI.upper<-psi+SE*qnorm(conf.level)
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

















#naive cross-fit plug-in method (i.e., CVonestep_bisec without one-step correction)
#arguments see CVonestep_bisec
CVnaive_bisec<-function(A,X,Y,scores,tau.interval=c(0,1),error.bound=0.05,conf.level=.95,nfolds=5,g.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),Q.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),g.trunc=1e-2,tol=1e-2){
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
    
    #nuisance estimators
    Z<-do.call(cbind,lapply(tau.interval,function(tau){
        as.numeric(scores<tau)
    }))
    g<-CV.est.g(A,X,folds,g.control,g.trunc)
    Q<-CV.est.Q(A,X,Z,folds,Q.control)
    
    #calculate coverage error estimator and asymptotic variance for each fold and the end point of tau.interval
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
            psi.v[v,tau.index]<-naive.psi
            sigma2.v[v,tau.index]<-mean(IF^2)
        }
    }
    psi<-t(fold.sizes)%*%psi.v/n
    psi<-pmin(pmax(psi,0),1) #project psi to [0,1]
    sigma2<-t(fold.sizes)%*%sigma2.v/n
    SE<-sqrt(sigma2/n)
    CI.upper<-psi+SE*qnorm(conf.level)
    
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
            psi.v[v,1]<-naive.psi
            sigma2.v[v,1]<-mean(IF^2)
        }
        
        new.psi<-as.numeric(t(fold.sizes)%*%psi.v/n)
        new.psi<-pmin(pmax(new.psi,0),1) #project psi to [0,1]
        new.sigma2<-as.numeric(t(fold.sizes)%*%sigma2.v/n)
        new.SE<-sqrt(new.sigma2/n)
        new.CI.upper<-new.psi+new.SE*qnorm(conf.level)
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
    
    if(sign.l<0){
        list(tau=tau.l,error.CI.upper=CI.upper.l,error.est=psi.l)
    }else{
        list(tau=tau.u,error.CI.upper=CI.upper.u,error.est=psi.u)
    }
}





#naive cross-fit plug-in method (i.e., CVonestep_bisec2 without one-step correction)
#arguments see CVonestep_bisec2
CVnaive_bisec2<-function(A,X,Y,scores,error.bound=0.05,conf.level=.95,nfolds=5,g.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),Q.control=list(SL.library=c("SL.glm","SL.gam","SL.randomForest")),g.trunc=1e-2,tol=1e-2){
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
    
    #nuisance estimators
    Z<-do.call(cbind,lapply(tau.index.interval,function(tau.index){
        as.numeric(scores<sorted.scores[tau.index])
    }))
    g<-CV.est.g(A,X,folds,g.control,g.trunc)
    Q<-CV.est.Q(A,X,Z,folds,Q.control)
    
    #calculate coverage error estimator and asymptotic variance for each fold and the end points of tau
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
            psi.v[v,tau.index]<-naive.psi
            sigma2.v[v,tau.index]<-mean(IF^2)
        }
    }
    psi<-t(fold.sizes)%*%psi.v/n
    psi<-pmin(pmax(psi,0),1) #project psi to [0,1]
    sigma2<-t(fold.sizes)%*%sigma2.v/n
    SE<-sqrt(sigma2/n)
    CI.upper<-psi+SE*qnorm(conf.level)
    
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
            psi.v[v,1]<-naive.psi
            sigma2.v[v,1]<-mean(IF^2)
        }
        
        new.psi<-as.numeric(t(fold.sizes)%*%psi.v/n)
        new.psi<-pmin(pmax(new.psi,0),1) #project psi to [0,1]
        new.sigma2<-as.numeric(t(fold.sizes)%*%sigma2.v/n)
        new.SE<-sqrt(new.sigma2/n)
        new.CI.upper<-new.psi+new.SE*qnorm(conf.level)
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
