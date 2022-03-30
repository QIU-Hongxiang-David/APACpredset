logit<-binomial()$linkfun
expit<-binomial()$linkinv

#create k folds of a vector id with A distributed roughly evenly in each fold
create.folds<-function(A,k){
    if(sum(A)<k || sum(1-A)<k){
        warning("Too few observations from source or target population for sample splitting!")
    }
    id<-1:length(A)
    id1<-id[A==1]
    id0<-id[A==0]
    
    order1<-sample.int(length(id1))
    d1<-suppressWarnings(data.frame(cbind(id1[order1],1:k)))
    names(d1)<-c("id","fold.id")
    
    order0<-sample.int(length(id0))
    d0<-suppressWarnings(data.frame(cbind(id0[order0],k:1)))
    names(d0)<-c("id","fold.id")
    
    d<-rbind(d0,d1)
    lapply(tapply(d$id,d$fold.id,identity,simplify=FALSE),sort)
}

#split data into training and testing sets with desired proportion
data.split<-function(A,train.prop){
    n<-length(A)
    if(sum(A)<2 || sum(1-A)<2){
        warning("Too few observations from source or target population for training data for sample splitting!")
    }
    id<-1:length(A)
    id1<-id[A==1]
    id0<-id[A==0]
    train<-sort(c(sample(id1,max(1,round(length(id1)*train.prop))),sample(id0,max(1,round(length(id0)*train.prop)))))
    test<-sort(setdiff(id,train))
    list(train=train,test=test)
}

#initial estimate of nuisance functions g and Q evaluated at data in each fold
CV.est.g<-function(A,X,folds,g.control,g.trunc){
    lapply(1:length(folds),function(v){
        A.train<-A[-folds[[v]]]
        if(all(A.train==0)){
            rep(0,length(folds[[v]]))
        }else if(all(A.train==1)){
            rep(1,length(folds[[v]]))
        }else{
            X.train<-X[-folds[[v]],,drop=FALSE]
            X.test<-X[folds[[v]],,drop=FALSE]
            args<-c(list(Y=A.train,X=X.train,newX=X.test,family=binomial()),g.control)
            SL.model<-do.call(SuperLearner,args)
            pmax(as.numeric(predict(SL.model)$pred),g.trunc)
        }
    })
}
CV.est.Q<-function(A,X,Z,folds,Q.control){
    n.tau<-ncol(Z)
    lapply(1:length(folds),function(v){
        A.train<-A[-folds[[v]]]
        X.train<-X[-folds[[v]],,drop=FALSE]
        X.test<-X[folds[[v]],,drop=FALSE]
        
        out<-matrix(nrow=length(folds[[v]]),ncol=n.tau)
        
        for(tau.index in 1:n.tau){
            Z.train<-Z[-folds[[v]],tau.index]
            if(all(Z.train[A.train==1]==0)){
                out[,tau.index]<-rep(0,length(folds[[v]]))
            }else if(all(Z.train[A.train==1]==1)){
                out[,tau.index]<-rep(1,length(folds[[v]]))
            }else{
                args<-c(list(Y=Z.train[A.train==1],X=X.train[A.train==1,,drop=FALSE],newX=X.test,family=binomial()),Q.control)
                SL.model<-do.call(SuperLearner,args)
                out[,tau.index]<-predict(SL.model)$pred
            }
        }
        out
    })
}

#find the index i of x such that x[j] (j<=i) are all TRUE
find.max.true<-function(x){
    for(i in 1:length(x)){
        if(!x[i]){
            return(i-1)
        }
    }
    length(x)
}
