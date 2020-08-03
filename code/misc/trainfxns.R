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

extract.coefs <- function(caret_obj,p.meth){
  # for given caret prediction method string p.meth, extract variable importance from model object
  # INPUTS:
  # caret_obj: caret model object output from train
  # p.meth: character specifying caret prediction method
  #
  # OUTPUTS:
  # 1xN matrix of coefficients. rf: random forest importance, lm: standardized betas where p <0.05 regression
  finalModel <- caret_obj$finalModel
  if(p.meth=='rf'){coefs <- t(finalModel$importance)/min(finalModel$importance) - 1 # normalize to minimum
  } else if(p.meth=='lm'){
    fm <- summary(lm.beta(finalModel))
    coefs <- t(fm$coefficients[,'Standardized',drop=FALSE] * (fm$coefficients[,'Pr(>|t|)',drop=FALSE] < 0.05)) 
  } else if(p.meth == 'glmnet'){
    coefs <- t(as.matrix(coef(finalModel,finalModel$lambdaOpt)))
  }
  return(coefs)
}

extract.bestTune.metric <- function(caret_obj,metric='Rsquared'){
  # INPUTS:
  # caret_obj: output of caret train() function
  # 
  # OUTPUTS:
  # get distribution of performances specified in metric over multiple resamples

  hyperparams <- caret_obj$bestTune # what tuning parameters were tried, select best
  res <- caret_obj$resample
  for(param in names(hyperparams)){
    res <- res[res[,param] == hyperparams[[param]],]  # subset results by those using parameter in best tune
  }
  return(res[,metric])  
}