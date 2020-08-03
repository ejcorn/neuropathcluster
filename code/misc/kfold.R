# this code creates a modified version of the glm and rf methods to use with caret train function
# this version downsamples to predetermined class balance during training
# but uses original class balance during testing

library(caret) # load caret package
# start with original glm
ds.glm <- getModelInfo(model='glm',regex=FALSE)[[1]]
# copy and paste original glm$fit() function, but add a line to class balance
ds.glm$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
  dat <- if(is.data.frame(x)) x else as.data.frame(x)

  # downsample and name y-variable '.outcome'
  # if there is a column in dat called 'WCB_LABELS', then do a multi class balance instead
  if('WCB_LABELS' %in% colnames(dat)){
    
    dat <- cbind(dat,OUTCOME_TMP=as.factor(y)) # add real y as outcome
    # now remove WCB_LABELS from dat, feed it in as the class balancing factor, label it to remove in next line
    dat <- downSample(dat[,setdiff(colnames(dat),'WCB_LABELS')],as.factor(dat[,'WCB_LABELS']),yname='.REMOVE') # 
    colnames(dat)[which(colnames(dat)=='OUTCOME_TMP')] <- '.outcome' # rename y .outcome for caret
    dat <- dat[,setdiff(colnames(dat),'.REMOVE')]
  } else {dat <- downSample(dat,as.factor(y),yname='.outcome')}

  if(length(levels(y)) > 2) stop("glm models can only use 2-class outcomes")

  theDots <- list(...)
  if(!any(names(theDots) == "family"))
  {
    theDots$family <- if(is.factor(y)) binomial() else gaussian()
  }

  ## pass in any model weights
  if(!is.null(wts)) theDots$weights <- wts
  modelArgs <- c(list(formula = as.formula(".outcome ~ ."), data = dat), theDots)

  out <- do.call("glm", modelArgs)
  ## When we use do.call(), the call infformation can contain a ton of
  ## information. Inlcuding the contenst of the data. We eliminate it.
  out$call <- NULL
  out
}

# glm$predict() is untouched

# now modify twoClassSummary function so that positive level can be specified

customTwoClassSummary <- function (data, lev = NULL, model = NULL, positive = 'X1',negative='X0') 
{
    lvls <- levels(data$obs)
    if (length(lvls) > 2) 
        stop(paste("Your outcome has", length(lvls), "levels. The twoClassSummary() function isn't appropriate."))
    caret:::requireNamespaceQuietStop("ModelMetrics")
    if (!all(levels(data[, "pred"]) == lvls)) 
        stop("levels of observed and predicted data do not match")
    rocAUC <- ModelMetrics::auc(ifelse(data$obs == lev[2], 0, 
        1), data[, lvls[1]])
    out <- c(rocAUC, sensitivity(data[, "pred"], data[, "obs"], 
        positive = positive), specificity(data[, "pred"], data[, "obs"], negative = negative),
        precision(data[,"pred"], data[,"obs"],relevant=positive),
        recall(data[,"pred"], data[,"obs"],relevant=positive))
    names(out) <- c("ROC", "Sens", "Spec",'Prec','Rec')
    out
}

kfold.GLM <- function(x,y,nreps=10,k.folds=3){
  # y has k columns of different labels, i.e. different clusters or different diseases
  # x contains all predictor variables, and an optional column called WCB_LABELS to enforce 
  # class balancing in the training process
  # x will be split into training and testing samples, then the training sample is downsampled to achieve class balance
  # k is the number of folds
  # nreps is the number of repeats of k-fold
  # 
  #x=df;y=dx;nreps=10;k.folds=3;p.ds=0.5
  #x=df;y=clusters;nreps=10;k.folds=10;p.ds=0.5
  data_ctrl <- trainControl(method = "repeatedcv", number = k.folds, repeats=nreps,
                          summaryFunction= function(...) customTwoClassSummary(...), classProbs= TRUE,
                          savePredictions='all')
  k <- ncol(y)
  tte <- list()  
  for(k.i in 1:k){
    n.i <- colnames(y)[k.i] # k are clusters or diseases
    print(paste('Training',n.i))
    y.i <- y[,k.i] # select class labels for particular one vs. all cluster/disease prediction task
    df.glm <- cbind(data.frame(y=make.names(as.factor(y.i))),x[,-1]) # need to remove first column of x (INDDIDs)
    model_caret <- train(y ~ .,   # model to fit
                         data = df.glm, metric = 'ROC',
                         trControl = data_ctrl, method = ds.glm) 
    
    # save performance of each fold and the feature weights of the best model
    tte[[n.i]] <- kfold.eval.glm(model_caret,y.i)

  }
  return(tte)
}

kfold.eval.glm <- function(model_caret,y){    
	# model_caret: caret training object using repeated k-fold cross validation
	# y: full sample of ground truth labels used for training
	tte <- list()
  	# get names of each rep-fold combo:
  	rep.fold.names <- unique(model_caret$pred$Resample)
  	k.folds <- model_caret$control$number
  	n.reps <- model_caret$control$repeats

	for(i in 1:(n.reps*k.folds)){
		
    	tte[[i]] <- list()
    	tte[[i]]$Sensitivity <- model_caret$resample$Sens[i]
    	tte[[i]]$Specificity <- model_caret$resample$Spec[i]
    	tte[[i]]$AUC <- model_caret$resample$ROC[i]
      tte[[i]]$Precision <- model_caret$resample$Prec[i]
      tte[[i]]$Recall <- model_caret$resample$Rec[i]
    	# get ROC curve for entire Rep, rather than each fold
    	# get names of all folds for each reps
		  rep.fold.names.k <- rep.fold.names[roundDownToMultiple(i,k.folds):roundUpToMultiple(i,k.folds)]
    	# extract sample indices for each rep
    	resampleIndex <- model_caret$pred$rowIndex[model_caret$pred$Resample %in% rep.fold.names.k]
    	# extract predicted values for positive class, which is X1
    	test.pred <- model_caret$pred$X1[model_caret$pred$Resample %in% rep.fold.names.k]
    	# extract ground truth values for each rep in correct order
    	test.ground_truth <- y[resampleIndex]
    	# compute ROC curve for each fold
    	roc.obj <- roc(response=test.ground_truth,predictor = as.numeric(test.pred))
      tte[[i]]$roc.obj <- roc.obj
      # get the index of best Youdens J (sens+spec) for choosing the right threshold
      # ideally would do this in training set
      YJ.index <- which.max(roc.obj$sensitivities + roc.obj$specificities - 1) 
      tte[[i]]$Sensitivity.YJ <- roc.obj$sensitivities[YJ.index]
      tte[[i]]$Specificity.YJ <- roc.obj$specificities[YJ.index]
      tte[[i]]$Threshold.YJ <- roc.obj$thresholds[YJ.index]
    	# compute PR curve for each fold      
      tte[[i]]$pr.obj <- pr.curve(scores.class0 = test.pred[test.ground_truth==1],scores.class1 = test.pred[test.ground_truth==0],curve=T)
      tte[[i]]$AUPRC <- tte[[i]]$pr.obj$auc.integral
    	# get weights for each fold for rep then average later
    	# extract feature weights by standardizing only continuous predictors
      n <- max(model_caret$pred$rowIndex) # get number of samples
      # get indices of held out data for fold j of rep i
    	held.out.idx <- model_caret$pred$rowIndex[model_caret$pred$Resample %in% rep.fold.names[i]]
      # get indices of training data for fold j of rep i (i.e. folds 1:k excluding j)
      resampleIndex <- c(1:n)[!1:n %in% held.out.idx]
      # get training data set for each fold/rep and downsample
      df.train <- model_caret$trainingData[resampleIndex,]
      df.train <- df.train[,setdiff(colnames(df.train),'WCB_LABELS')]
      df.train <- standardize.continuous(downSample(df.train,df.train$.outcome,list=FALSE,yname='.outcome'))
    	tte[[i]]$FeatureImportance <- 
    		coefficients(glm(.outcome~.,data=df.train,family='binomial'))
    }
    return(tte)
}

#####################
### Random Forest ###
#####################

# start with original random forest
ds.rf <- getModelInfo('rf',regex=FALSE)[[1]]
ds.rf$fit <- function(x, y, wts, param, lev, last, classProbs, ...){
  # downsample x and y before training

  if('WCB_LABELS' %in% colnames(x)){    
    dat <- cbind(as.data.frame(x),OUTCOME_TMP=as.factor(y)) # add real y as outcome
    # now remove WCB_LABELS from dat, feed it in as the class balancing factor, label it to remove in next line
    dat <- downSample(dat[,setdiff(colnames(dat),'WCB_LABELS')],as.factor(dat[,'WCB_LABELS']),yname='.REMOVE') # 
    colnames(dat)[which(colnames(dat)=='OUTCOME_TMP')] <- 'y' # rename y .outcome for caret
    df <- dat[,setdiff(colnames(dat),'.REMOVE')]
    #print(df)
  } else {df <- downSample(x,as.factor(y),yname='y')}
  
  if(length(unique(df[,'y'])) ==1){df <- downSample(x,as.factor(y),yname='y')} # if weighted downsampling got rid of positive class, use binary downsampling

  x <- subset(df,select= -c(y))
  y <- df$y
  randomForest::randomForest(x, y, importance = TRUE, norm.votes = TRUE, proximity = TRUE,keep.forest=TRUE,mtry = param$mtry, ...)
}

# ds.rf$predict() remains untouched so no downsampling is applied to test set

kfold.RF <- function(x,y,nreps=10,k.folds=3,ntree=500,mtry=ncol(x)-1){
  # y has k columns of different labels, i.e. different clusters or different diseases
  # x contains all predictor variables
  # x will be split into training and testing samples, then the training sample is downsampled to achieve class balance
  # k is the number of folds
  # nreps is the number of repeats of k-fold
  # 
  #x=df;y=dx;nreps=10;k.folds=3;p.ds=0.5;ntree=500;mtry=ncol(x)-1
  #x=df;y=clusters;nreps=10;k.folds=10;p.ds=0.5
  if('WCB_LABELS' %in% colnames(x)){mtry <- mtry-1} # not a real variable so exclude from max mtry
  data_ctrl <- trainControl(method = "repeatedcv", number = k.folds, repeats=nreps,
                          summaryFunction= function(...) customTwoClassSummary(...), classProbs= TRUE,
                          savePredictions='all',returnResamp = 'all')
  k <- ncol(y)
  tte <- list()  
  for(k.i in 1:k){
    n.i <- colnames(y)[k.i] # k are clusters or diseases
    print(paste('Training',n.i))
    y.i <- y[,k.i] # select class labels for particular one vs. all cluster/disease prediction task
    df.rf <- cbind(data.frame(y=make.names(as.factor(y.i))),x[,-1]) # need to remove first column of x (INDDIDs)
    model_caret <- train(y ~ .,   # model to fit
                         data = df.rf, metric = 'ROC', tuneGrid = data.frame(.mtry=mtry),
                         trControl = data_ctrl, method = ds.rf,ntree=500,returnResamp=TRUE) 
    
    # save performance of each fold and the feature weights of the best model
    tte[[n.i]] <- kfold.eval.RF(model_caret,y.i)

  }
  return(tte)
}

kfold.eval.RF <- function(model_caret,y){    
  # model_caret: caret training object using repeated k-fold cross validation
  # y: full sample of ground truth labels used for training
  tte <- list()
  # get names of each rep-fold combo:
  rep.fold.names <- unique(model_caret$pred$Resample)
  k.folds <- model_caret$control$number
  n.reps <- model_caret$control$repeats
  # get feature importance for best model... ideally would get this for every rep and average
  # but seems impossible to do with caret
  finalImportance <- model_caret$finalModel$importance[,'MeanDecreaseAccuracy']

  for(i in 1:(n.reps*k.folds)){
    
      tte[[i]] <- list()
      tte[[i]]$Sensitivity <- model_caret$resample$Sens[i]
      tte[[i]]$Specificity <- model_caret$resample$Spec[i]
      tte[[i]]$AUC <- model_caret$resample$ROC[i]
      tte[[i]]$Precision <- model_caret$resample$Prec[i]
      tte[[i]]$Recall <- model_caret$resample$Rec[i]      
      # get ROC curve for entire Rep, rather than each fold
      # get names of all folds for each reps
    rep.fold.names.k <- rep.fold.names[roundDownToMultiple(i,k.folds):roundUpToMultiple(i,k.folds)]
      # extract sample indices for each rep
      resampleIndex <- model_caret$pred$rowIndex[model_caret$pred$Resample %in% rep.fold.names.k]
      # extract predicted values for positive class
      test.pred <- model_caret$pred$X1[model_caret$pred$Resample %in% rep.fold.names.k]
      # extract ground truth values for each rep in correct order
      test.ground_truth <- y[resampleIndex]

      # compute ROC curve for each fold
      roc.obj <- roc(response=test.ground_truth,predictor = as.numeric(test.pred))
      tte[[i]]$roc.obj <- roc.obj
      YJ.index <- which.max(roc.obj$sensitivities + roc.obj$specificities - 1) 
      tte[[i]]$Sensitivity.YJ <- roc.obj$sensitivities[YJ.index]
      tte[[i]]$Specificity.YJ <- roc.obj$specificities[YJ.index]
      tte[[i]]$Threshold.YJ <- roc.obj$thresholds[YJ.index]
      # compute PR curve for each fold      
      #tte[[i]]$pr.obj <- pr.curve(scores.class0 = test.pred[test.ground_truth==1],scores.class1 = test.pred[test.ground_truth==0],curve=T)
      #tte[[i]]$AUPRC <- tte[[i]]$pr.obj$auc.integral
      # get feature importance for best model... ideally would get this for every rep and average
      # but seems impossible to do with caret
      tte[[i]]$FeatureImportance <- finalImportance
    }
    return(tte)
}


caret_oos <- function(model_caret,y){    
  # model_caret: caret training object using repeated k-fold cross validation
  # y: full sample of ground truth labels used for training
  tte <- list()
  # get names of each rep-fold combo:
  rep.fold.names <- unique(model_caret$pred$Resample)
  k.folds <- model_caret$control$number
  n.reps <- model_caret$control$repeats  

  for(i in 1:(n.reps*k.folds)){
      # get ROC curve for entire Rep, rather than each fold
      # get names of all folds for each reps
      tte[[i]] <- list()
      rep.fold.names.k <- rep.fold.names[roundDownToMultiple(i,k.folds):roundUpToMultiple(i,k.folds)]
      # extract sample indices for each rep
      resampleIndex <- model_caret$pred$rowIndex[model_caret$pred$Resample %in% rep.fold.names.k]
      # extract predicted values for positive class
      tte[[i]]$test.pred <- model_caret$pred$pred[model_caret$pred$Resample %in% rep.fold.names.k]
      # extract ground truth values for each rep in correct order
      tte[[i]]$test.ground_truth <- y[resampleIndex]
    }
  return(tte)
}
