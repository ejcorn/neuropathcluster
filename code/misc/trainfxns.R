classbalance.downsample <- function(y,p){
  # y: binary outcome variable
  # p: desired proportion of 1's 
  # will downsample 
  y <- as.numeric(y)
  if(mean(y) < p){ # undersample negative class
    print('undersampling 0\'s')
    num.controls.to.remove <- sum(!y) - (sum(y)/p - sum(y))
    remove.idx <- sample(which(y == 0),size=num.controls.to.remove,replace = F)
  } else if(mean(y) > p) {
    print('undersampling 1\'s')
    num.positives.to.remove <- sum(y) - (sum(!y)/p - sum(!y))
    remove.idx <- sample(which(y == 1),size=num.positives.to.remove,replace = F)
  }
  return(remove.idx)
}


train.test.RF.class <- function(x,y,nreps=10,trainfrac=0.7,p.ds=0.5,ntree=500,mtry=ncol(x)-1){
  # y has k columns of different labels, i.e. different clusters or different diseases
  # x contains all predictor variables
  # p.ds is the class balance to achieve by downsampling
  # testing data
  #x=df;y=dx;nreps=10;trainfrac=0.7;p.ds=0.5;ntree=500;mtry=ncol(x)-1
  k <- ncol(y)
  tte <- list()  
  for(k.i in 1:k){
    n.i <- colnames(y)[k.i]
    for(i in 1:nreps){
      y.i <- y[,k.i]
      remove.idx <- classbalance.downsample(y.i,p.ds)
      df <- x[-remove.idx,]
      y.i <- y.i[-remove.idx]
      list[df.train,y.i.train,df.test,y.i.test,train.idx] <- train.test.split(x=df,y=y.i,trainfrac=trainfrac)
      m.t <- randomForest(y=as.factor(y.i.train),
                          x=df.train,
                          ntree=ntree,
                          mtry=mtry,
                          importance = TRUE, norm.votes = TRUE, proximity = TRUE,keep.forest=TRUE)
      test.ground_truth <- factor(y.i.test,levels=c(0,1))
      train.ground_truth <- factor(y.i.train,levels=c(0,1))
      test.pred <- predict(m.t,as.matrix(df.test),type="prob")[,'1']
      train.pred <- predict(m.t,as.matrix(df.train),type="prob")[,'1']
      # tte stores the results of every model run
      #tte[[n.i]][[i]] <- train.test.eval(test.pred,test.ground_truth,train.pred,train.ground_truth)
      tte[[n.i]][[i]] <- test.eval(test.pred,test.ground_truth)
      tte[[n.i]][[i]]$FeatureImportance <- m.t$importance[,'MeanDecreaseAccuracy']
    }
  }
  return(tte)
}

undo.dummy.df <- function(df){
  # go from dummy df (k non-overlapping binary columns)
  # to a single column with k classes
  multiclass.df <- sapply(1:nrow(df), function(i) colnames(df)[which(as.logical(df[i,]))])
  return(multiclass.df)
}

train.test.RF.multiclass <- function(x,y,nreps=10,trainfrac=0.7,p.ds=0.5,ntree=500,mtry=ncol(x)-1){
  # y has k columns of different labels, i.e. different clusters or different diseases
  # x contains all predictor variables
  # p.ds is the class balance to achieve by downsampling
  # testing data
  #x=df;y=dx;nreps=10;trainfrac=0.7;p.ds=0.5;ntree=500;mtry=ncol(x)-1
  y <- data.frame(y=undo.dummy.df(y))
  k <- ncol(y)
  tte <- list()  
  for(i in 1:nreps){
    list[df.train,y.i.train,df.test,y.i.test,train.idx] <- train.test.split.multiclass(x=df,y=y,trainfrac=trainfrac)
    m.t <- randomForest(y=as.factor(y.i.train),
                        x=df.train,
                        ntree=ntree,
                        mtry=mtry,
                        importance = TRUE, norm.votes = TRUE, proximity = TRUE,keep.forest=TRUE)
    test.ground_truth <- factor(y.i.test,levels=c(0,1))
    train.ground_truth <- factor(y.i.train,levels=c(0,1))
    test.pred <- predict(m.t,as.matrix(df.test),type="prob")[,'1']
    train.pred <- predict(m.t,as.matrix(df.train),type="prob")[,'1']
    # tte stores the results of every model run
    #tte[[n.i]][[i]] <- train.test.eval(test.pred,test.ground_truth,train.pred,train.ground_truth)
    tte[[n.i]][[i]] <- test.eval(test.pred,test.ground_truth)
    tte[[n.i]][[i]]$FeatureImportance <- m.t$importance[,'MeanDecreaseAccuracy']
  }
  return(tte)
}

train.test.split <- function(x,y,trainfrac){
  # assumes first column of x is INDDID
  # split into training and testing samples
  train.idx <- createDataPartition(x$INDDID,times = 1,p = trainfrac)$Resample1
  return(list(df.train=x[train.idx,-1],y.i.train=y[train.idx],df.test=x[-train.idx,-1],y.i.test=y[-train.idx],train.idx=train.idx))
}

train.test.split.multiclass <- function(x,y,trainfrac){
  # assumes first column of x is INDDID
  # split into training and testing samples
  train.idx <- createDataPartition(x$INDDID,times = 1,p = trainfrac)$Resample1
  return(list(df.train=x[train.idx,-1],y.i.train=y[train.idx,],df.test=x[-train.idx,-1],y.i.test=y[-train.idx,],train.idx=train.idx))
}

test.eval <- function(test.pred,test.ground_truth){
  # perform a series of evaluations of model performance on test set
  # test.ground_truth: factor vector of binary true class labels, for test set
  # test.pred: numeric vector of probabilities, for test set

  tte <- list()  
  # confusionMatrix for threshold of 0.5  
  tte$test.ground_truth <- test.ground_truth  
  tte$test.pred <- test.pred

  # roc.curve for test set
  roc.obj <- roc(response=test.ground_truth,predictor = as.numeric(test.pred))
  tte$roc.obj <- roc.obj
  tte$AUC <- auc(roc.obj)   

  tte$test.pred.factor <- factor(round(test.pred),levels=c(0,1))  
  tte$OverallAccuracy <- mean(tte$test.pred.factor == test.ground_truth)
  tte$cm <- confusionMatrix(positive='1',data=tte$test.pred.factor,reference=test.ground_truth)
  tte$Sensitivity <- tte$cm$byClass['Sensitivity']
  tte$Specificity <- tte$cm$byClass['Specificity']
  
  # find best accuracy on test set
  thrsh.rng <- seq(0,1,length.out=100)
  acc.by.thrsh <- sapply(thrsh.rng, function(thrsh.i) 
    mean(factor(as.numeric(test.pred > thrsh.i),levels=c(0,1)) == test.ground_truth))
  tte$test.acc.best <- max(acc.by.thrsh)
  tte$test.thrsh.best <- thrsh.rng[which.max(acc.by.thrsh)]
  return(tte)

}


train.test.eval <- function(test.pred,test.ground_truth,train.pred,train.ground_truth){
  # perform a series of evaluations of model performance on test set
  # based on training set thresholds
  # test.ground_truth: factor vector of binary true class labels, for test set
  # test.pred: numeric vector of probabilities, for test set
  # train.ground_truth: factor vector of binary true class labels, for training set
  # train.pred: numeric vector of probabilities, for training set
  tte <- list()
  tte$test.ground_truth <- test.ground_truth
  tte$train.ground_truth <- train.ground_truth
  tte$test.pred <- test.pred
  tte$train.pred <- train.pred

  # confusionMatrix for threshold of 0.5
  tte <- confusionMatrix(positive='1',as.factor(round(test.pred)),test.ground_truth)
  
  # roc.curve for test set
  roc.obj <- roc(response=test.ground_truth,predictor = as.numeric(test.pred))
  tte$roc.obj <- roc.obj
  tte$AUC <- auc(roc.obj)   

  # find threshold giving highest mean(sensitivity and specificity) on TRAINING set
  roc.train <- roc(response=train.ground_truth,predictor = as.numeric(train.pred))
  snspmean.idx <- which.max(0.5*(roc.train$sensitivities + roc.train$specificities))
  train.thresh.snspmean <- roc.train$thresholds[snspmean.idx]
  print(paste('Training Threshold:',train.thresh.snspmean))
  # find sens., spec., accuracy values on TEST set using best threshold from TRAINING set  
  tte$test.pred.best <- as.factor(as.numeric(test.pred > train.thresh.snspmean))
  tte$OverallAccuracy.snspmean <- mean(tte$test.pred.best == test.ground_truth)
  cm <- confusionMatrix(positive='1',data=tte$test.pred.best,reference=test.ground_truth)
  tte$Sensitivity.snspmean <- cm$byClass['Sensitivity']
  tte$Specificity.snspmean <- cm$byClass['Specificity']

  # find best possible accuracy on test set
  thrsh.rng <- seq(0,1,length.out=100)
  acc.by.thrsh <- sapply(thrsh.rng, function(thrsh.i) 
    as.factor(as.numeric(test.pred > thrsh.i)) == test.ground_truth)
  tte$test.acc.best <- max(acc.by.thrsh)
  tte$test.thrsh.best <- thrsh.rng[which.max(acc.by.thrsh)]
  return(tte)

}

tune.RF.grid <- function(x,y){
  control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
  # random seed based on current time in msec
  s = (print(as.numeric(Sys.time())*1000,digits=15) %% 1)* 1000
  set.seed(s)
  metric <- 'Accuracy'
  tunegrid <- expand.grid(.mtry=c(1:(ncol(x)-1)))
  y <- as.factor(y)
  dataset <- cbind(data.frame(y=y),x)
  rf_gridsearch <- train(y~., data=dataset, method="rf", metric=metric, tuneGrid=tunegrid, trControl=control)
  return(rf_gridsearch)
}

train.tune.test.RF.class <- function(x,y,nreps=10,trainfrac=0.7,p.ds=0.5,ntree=500,mtry=ncol(x)-1){
  # y has k columns of different labels, i.e. different clusters or different diseases
  # x contains all predictor variables
  # p.ds is the class balance to achieve by downsampling
  # testing data
  #x=df;y=dx;nreps=10;trainfrac=0.7;p.ds=0.5;ntree=500;mtry=ncol(x)-1
  k <- ncol(y)
  tte <- list()
  for(k.i in 1:k){
    n.i <- colnames(y)[k.i]
    for(i in 1:nreps){
      y.i <- y[,k.i]
      remove.idx <- classbalance.downsample(y.i,p.ds)
      df <- x[-remove.idx,]
      y.i <- y.i[-remove.idx]
      list[df.train,y.i.train,df.test,y.i.test,train.idx] <- train.test.split(x=df,y=y.i,trainfrac=trainfrac)
      rf.tune <- tune.RF.grid(df.train,y.i.train)
      mtry <- as.numeric(rf.tune$bestTune['mtry'])
      print(paste('best mtry:',mtry))
      m.t <- randomForest(y=as.factor(y.i.train),
                          x=df.train,
                          ntree=ntree,
                          mtry=mtry,
                          importance = TRUE, norm.votes = TRUE, proximity = TRUE,keep.forest=TRUE)
      test.ground_truth <- factor(y.i.test,levels=c(0,1))
      train.ground_truth <- factor(y.i.train,levels=c(0,1))
      test.pred <- predict(m.t,as.matrix(df.test),type="prob")[,'1']
      train.pred <- predict(m.t,as.matrix(df.train),type="prob")[,'1']
      # tte stores the results of every model run
      #tte[[n.i]][[i]] <- train.test.eval(test.pred,test.ground_truth,train.pred,train.ground_truth)
      tte[[n.i]][[i]] <- test.eval(test.pred,test.ground_truth)
      tte[[n.i]][[i]]$FeatureImportance <- m.t$importance[,'MeanDecreaseAccuracy']      
    }
  }
  return(tte)
}

cv.lasso <- function(x,y){
  # choose lambda for lasso with binary outcome
  m <- cv.glmnet(x=as.matrix(x),y=as.matrix(y),family = 'binomial',standardize=TRUE)
  return(m)
}

lasso <- function(x,y,lambda){
  m <- glmnet(x=as.matrix(x),y=as.matrix(y),lambda=lambda,family = 'binomial',standardize=TRUE)
  return(m)
}


train.test.Lasso <- function(x,y,nreps=10,trainfrac=0.7,p.ds=0.5,ntree=500,mtry=ncol(x)-1){
  # y has k columns of different labels, i.e. different clusters or different diseases
  # x contains all predictor variables
  # p.ds is the class balance to achieve by downsampling
  # testing data
  #x=df;y=dx;nreps=10;trainfrac=0.7;p.ds=0.5
  #x=df;y=clusters;nreps=100;trainfrac=0.7;p.ds=0.5
  k <- ncol(y)
  tte <- list()  
  for(k.i in 1:k){
    n.i <- colnames(y)[k.i]
    for(i in 1:nreps){
      y.i <- y[,k.i]
      remove.idx <- classbalance.downsample(y.i,p.ds)
      df <- x[-remove.idx,]
      y.i <- y.i[-remove.idx]
      list[df.train,y.i.train,df.test,y.i.test,train.idx] <- train.test.split(x=df,y=y.i,trainfrac=trainfrac)
      cv <- cv.lasso(y=y.i.train,x=df.train)
      lambda <- cv$lambda.min      
      m.t <- lasso(y=as.factor(y.i.train),
               x=df.train,
               lambda=lambda)      
      test.ground_truth <- factor(y.i.test,levels=c(0,1))
      train.ground_truth <- factor(y.i.train,levels=c(0,1))
      # type = 'response' converts predicted log odds ratio to a probability
      test.pred <- predict(m.t,as.matrix(df.test),type='response')
      train.pred <- predict(m.t,as.matrix(df.train),type='response')
      # tte stores the results of every model run
      #tte[[n.i]][[i]] <- train.test.eval(test.pred,test.ground_truth,train.pred,train.ground_truth)
      tte[[n.i]][[i]] <- test.eval(test.pred,test.ground_truth)
      tte[[n.i]][[i]]$FeatureImportance <- coefficients(m.t)
    }
  }
  return(tte)
}

get.continuous.mask <- function(x){
  # x is a data frame which may contain both continuous and interval/categorical predictors
  # For my purposes, a continuous predictor is one with more than 3 levels 
  # (allele tables have 3 levels, 0,1,2)
  return(sapply(1:ncol(x), function(i) length(unique(x[,i]))) > 3)
}
standardize.continuous <- function(x){
  # x is a data frame which may contain both continuous and interval/categorical predictors
  # this function will standardize only the continuous predictors
  cont.mask <- get.continuous.mask(x)
  x[,cont.mask] <- scale(x[,cont.mask],center=T)
  return(x)
}

train.test.old.GLM <- function(x,y,nreps=10,trainfrac=0.7,p.ds=0.5){
  # y has k columns of different labels, i.e. different clusters or different diseases
  # x contains all predictor variables
  # p.ds is the class balance to achieve by downsampling

  #x=df;y=dx;nreps=10;trainfrac=0.7;p.ds=0.5
  #x=df;y=clusters;nreps=10;trainfrac=0.7;p.ds=0.5
  k <- ncol(y)
  tte <- list()  
  for(k.i in 1:k){
    n.i <- colnames(y)[k.i]
    for(i in 1:nreps){
      y.i <- y[,k.i]
      remove.idx <- classbalance.downsample(y.i,p.ds)
      df <- x[-remove.idx,]
      y.i <- y.i[-remove.idx]
      list[df.train,y.i.train,df.test,y.i.test,train.idx] <- train.test.split(x=df,y=y.i,trainfrac=trainfrac)
      df.glm <- cbind(data.frame(y=y.i.train),df.train)
      m.t <- glm(y~.,data=df.glm,family='binomial')
      test.ground_truth <- factor(y.i.test,levels=c(0,1))
      train.ground_truth <- factor(y.i.train,levels=c(0,1))
      test.pred <- predict(m.t,df.test,type='response')
      train.pred <- predict(m.t,df.train,type='response')
      tte[[n.i]][[i]] <- test.eval(test.pred,test.ground_truth)
      #tte[[n.i]][[i]]$vif <- vif(m.t) # compute variance inflation factor for each model
      # extract feature weights by standardizing only continuous predictors
      tte[[n.i]][[i]]$FeatureImportance <-
        coefficients(glm(y~.,data=standardize.continuous(df.glm),family='binomial'))

    }
  }
  return(tte)
}

train.test.GLM <- function(x,y,nreps=10,trainfrac=0.7,p.ds=0.5){
  # y has k columns of different labels, i.e. different clusters or different diseases
  # x contains all predictor variables
  # p.ds is the class balance to achieve by downsampling
  # this function evaluates performance on test set with original class balance
  # 
  #x=df;y=dx;nreps=10;trainfrac=0.7;p.ds=0.5
  #x=df;y=clusters;nreps=10;trainfrac=0.7;p.ds=0.5
  k <- ncol(y)
  tte <- list()  
  for(k.i in 1:k){
    n.i <- colnames(y)[k.i] # k are clusters or diseases
    for(i in 1:nreps){
      y.i <- y[,k.i] # select class labels for particular one vs. all cluster/disease prediction task
      # split into training and testing data sets
      list[df.train,y.i.train,df.test,y.i.test,train.idx] <- train.test.split(x=x,y=y.i,trainfrac=trainfrac)
      # balance training set by downsampling majority class
      remove.idx <- classbalance.downsample(y.i.train,p.ds)
      df.train <- df.train[-remove.idx,]
      y.i.train <- y.i.train[-remove.idx]
      # train glm to predict balanced class labels
      df.glm <- cbind(data.frame(y=y.i.train),df.train)
      m.t <- glm(y~.,data=df.glm,family='binomial')
      # evaluate model performance on test set with original class balance
      test.ground_truth <- factor(y.i.test,levels=c(0,1))
      test.pred <- predict(m.t,df.test,type='response')
      tte[[n.i]][[i]] <- test.eval(test.pred,test.ground_truth)      
      # extract feature weights by standardizing only continuous predictors
      tte[[n.i]][[i]]$FeatureImportance <-
        coefficients(glm(y~.,data=standardize.continuous(df.glm),family='binomial'))

    }
  }
  return(tte)
}
