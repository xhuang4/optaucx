#' Youden score for a specific cutpoint
#'
#' @param data input data frame.
#' @param yvar response variable name.
#' @param xvar composite score from optAUC.
#' @param dir direction of cut.
#' @param cutoff a specific cutpoint for which the score needs to be computed.
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff.
#'
#' @return Youden score for the cutpoint specified.
#'
seq.find.score.Youden<-function(data,yvar,xvar,dir,cutoff,min.sigp.prcnt)
{
  if (dir==">") {
    id <- data[, xvar]>cutoff
  } else {
    id <- data[, xvar]<cutoff
  }
  nsubj<-nrow(data)
  n.grp <- table(id)
  sigp.prcnt<-sum(id)/nsubj
  y<-data[,yvar]

  if ((sigp.prcnt>min.sigp.prcnt) & (length(n.grp)==2) & (length(table(y))==2))
  {
    table.roc<-as.matrix(table(id,y))
    se<-table.roc[2,2]/(table.roc[2,2]+table.roc[1,2])
    sp<-table.roc[1,1]/(table.roc[1,1]+table.roc[2,1])
    score<-se+sp
  } else {
    score<-NA
  }
  return(score)
}

#' Find cutoff for predictive case using Youden score.
#'
#' @param data input data frame.
#' @param yvar response variable name.
#' @param censorvar censoring variable name only for TTE data , censor=0,event=1)  - default = NULL.
#' @param xvar composite score from optAUC.
#' @param trtvar treatment variable name (for predictive signature). Set trt.vec to NULL for prognostic signature.
#' @param trtref coding (in the column of trtvar) for treatment arm.
#' @param type type of response variable: "s" survival; "b" binary (default).
#' @param dir direction of cut.
#' @param nsubj number of subjects.
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff.
#'
#' @return The cutpoint corresponding to maximum Youden score.
#'
seq.find.cutoff.Youden<-function(data, yvar, censorvar, xvar, trtvar, trtref, type, dir, nsubj, min.sigp.prcnt)
{
  cut.vec <- sort(unique(data[, xvar]))
  cut.vec <- quantile(cut.vec,  prob=seq(0.05, 0.95, 0.05),type=3)
  cut.vec=sort(unique(cut.vec))

  if (type=="b")
  {
    if (is.null(trtvar))
    {
      cut.score <- unlist(lapply(cut.vec, function(x) seq.find.score.Youden(data, yvar, xvar, dir, x, min.sigp.prcnt)))
      tmp<-max((cut.score),na.rm=T)
      cut.val <- cut.vec[which(cut.score==tmp)]
    } else {
      trt.id<-1*(data[,trtvar]==trtref)
      data.trt<-data[trt.id==1,]
      data.ctl<-data[trt.id==0,]
      cut.score.trt <- unlist(lapply(cut.vec, function(x) seq.find.score.Youden(data.trt, yvar, xvar, dir, x, min.sigp.prcnt)))
      cut.score.ctl <- unlist(lapply(cut.vec, function(x) seq.find.score.Youden(data.ctl, yvar, xvar, dir, x, min.sigp.prcnt)))
      cut.score<-cut.score.trt-cut.score.ctl
      tmp<-max((cut.score),na.rm=T)
      cut.val <- cut.vec[which(cut.score==tmp)]
    }
  }
  return(cut.val)
}

#' Function for obtaining cutpoints based on Youden score.
#'
#' @param data the input data frame.
#' @param yvar response variable name.
#' @param xvar composite score from optAUC.
#' @param trtvar treatment variable name (for predictive signature). Set trt.vec to NULL for prognostic signature.
#' @param trtref coding (in the column of trtvar) for treatment arm.
#' @param censorvar censoring variable name only for TTE data , censor=0,event=1)  - default = NULL.
#' @param type type of response variable: "s" survival; "b" binary (default).
#' @param n.boot number of bootstrap for batting or Youden procedure.
#' @param des.res the desired response. "larger": prefer larger response (default) "smaller": prefer smaller response.
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff.
#'
#' @return this function returns a list containing
#' \item{res}{prediction rule based on optAUC score.}
#' \item{train.stat}{subgroup p-value and group statistics for training data.}
#' \item{train.plot}{Treatment-subgroup interaction plot for training data (only for predictive signature).}
#'
Youden.func<-function(data,yvar,xvar,trtvar,trtref,censorvar,type,n.boot,des.res,min.sigp.prcnt)
{
  niter <- n.boot
  nsubj<-nrow(data)

  if (!is.null(trtvar))
  {
    if (type=="s") {
      model.text <- paste("Surv(", yvar, ", ", censorvar, ")~",trtvar,"+",trtvar, "*", xvar, sep="")
      model.formula <- as.formula(model.text)
      model.temp <- coxph(model.formula, data=data)
      coef <- try(summary(model.temp)$coefficient[3,1],silent=TRUE)
    } else {
      model.text <- paste(yvar, "~", trtvar,"+",trtvar, "*", xvar, sep="")
      model.formula <- as.formula(model.text)
      model.temp <- glm(model.formula, family=binomial,data=data)
      coef <- try(summary(model.temp)$coefficient[4,1],silent=TRUE)
    }
  } else {
    if (type=="s") {
      model.text <- paste("Surv(", yvar, ", ", censorvar, ")~", xvar, sep="")
      model.formula <- as.formula(model.text)
      model.temp <- coxph(model.formula, data=data)
      coef <- try(summary(model.temp)$coefficient[1],silent=TRUE)    ##  Postive Hazard means greater risk
    } else {
      model.text <- paste(yvar, "~", xvar, sep="")
      model.formula <- as.formula(model.text)
      model.temp <- glm(model.formula, family=binomial,data=data)
      coef <- try(summary(model.temp)$coefficient[2,1],silent=TRUE)
    }
  }

  if (class(coef)=="try-error") {
    dir <- NA
    model.pval <- NA
    cutoff.med <- NA
    cut.result <- c(xvar, dir, cutoff.med, model.pval)
    return(cut.result)
  }

  if(des.res=="smaller"){
    if (type=="s") dir=ifelse(coef>=0,">","<")
    if (type=="b") dir=ifelse(coef>=0,"<",">")
  } else {
    if (type=="s") dir=ifelse(coef>=0,"<",">")
    if (type=="b") dir=ifelse(coef>=0,">","<")
  }

  cutoff.vec <- NULL
  if (!is.na(dir)) {
    for(i in 1:niter) {
      train.id <- sample(1:nsubj, nsubj, replace=TRUE)
      data.train <- data[train.id,]
      cutoff.temp <- seq.find.cutoff.Youden(data=data.train, yvar=yvar, censorvar=censorvar, xvar=xvar,
                                            trtvar=trtvar, trtref=trtref, type=type, dir=dir,
                                            nsubj=nsubj, min.sigp.prcnt=min.sigp.prcnt)
      cutoff.vec <- c(cutoff.vec, cutoff.temp)
    }
    cutoff.med <- median(cutoff.vec, na.rm=TRUE)
    cut.result <- cbind(xvar, dir, cutoff.med)
    colnames(cut.result) <- c("variable", "direction", "threshold")
  }

  if(is.null(censorvar)) {
    censor.vec=NULL
  }else{
    censor.vec=data[censorvar]
  }
  if(is.null(trtvar)) {
    trt.vec=NULL
  }else{
    trt.vec=data[trtvar]
  }

  pred.data=pred.seqlr(data[xvar], cut.result)
  train.stat=evaluate.results(data[yvar], pred.data, censor.vec=censor.vec,trt.vec=trt.vec, trtref=trtref, type=type)

  if (!is.null(trtvar)){
    train.plot=interaction.plot(train.stat, type=type, main="Interaction Plot (optAUC, Train)", trt.lab=c("Trt.", "Ctrl."))
  } else {
    train.plot=NULL
  }
  return(list(res=cut.result,train.stat=train.stat,train.plot=train.plot))
}

#' Function for obtaining cutpoints using univariate BATTing.
#'
#' @param data the input data frame.
#' @param yvar response variable name.
#' @param xvar composite score from optAUC.
#' @param censorvar censoring variable name only for TTE data , censor=0,event=1)  - default = NULL.
#' @param trtvar treatment variable name (for predictive signature). Set trt.vec to NULL for prognostic signature.
#' @param trtref coding (in the column of trtvar) for treatment arm.
#' @param type type of response variable: "s" survival; "b" binary (default).
#' @param des.res the desired response. "larger": prefer larger response (default) "smaller": prefer smaller response.
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff.
#' @param n.boot number of bootstrap for batting or Youden procedure.
#'
#' @return this function returns a list containing
#' \item{res}{prediction rule based on optAUC score.}
#' \item{train.stat}{subgroup p-value and group statistics for training data.}
#' \item{train.plot}{Treatment-subgroup interaction plot for training data (only for predictive signature).}
#'
  BATTing.func<-function(data, yvar, xvar, censorvar, trtvar, trtref, type, des.res, min.sigp.prcnt, n.boot)
{
  if(is.null(censorvar)) {
    censor.vec=NULL
  }else{
    censor.vec=data[censorvar]
  }
  if(is.null(trtvar)) {
    trt.vec=NULL
  }else{
    trt.vec=data[trtvar]
  }

  res=seqlr.batting(y=data[yvar], x=data[xvar], censor.vec=censor.vec,
                    trt.vec=trt.vec, trtref=trtref, type=type,
                    n.boot=n.boot, des.res=des.res,
                    min.sigp.prcnt=min.sigp.prcnt)
  pred.data=pred.seqlr(data[xvar], res)
  train.stat=evaluate.results(data[yvar], pred.data, censor.vec=censor.vec,
                              trt.vec=trt.vec, trtref=trtref, type=type)
  if (!is.null(trtvar)){
    train.plot=interaction.plot(train.stat, type=type, main="Interaction Plot (optAUC, Train)", trt.lab=c("Trt.", "Ctrl."))
  } else {
    train.plot=NULL
  }
  res<-rbind(res[,-4])
  colnames(res)<-c("variable","direction","threshold")
  return(list(res=res,train.stat=train.stat,train.plot=train.plot))
}

#' Fit function for optAUCX, fit linear composite and find threshold for subgroup
#'
#' @param data the input data frame
#' @param yvar response variable name
#' @param yvar.ref response variable reference
#' @param xvars vector of predictor variable names
#' @param censorvar censoring variable name only for TTE data , censor=0,event=1)  - default = NULL.
#' @param trtvar treatment variable name (for predictive signature). Set trt.vec to NULL for prognostic signature.
#' @param type type of response variable: "s" survival; "b" binary (default).
#' @param trtref coding (in the column of trtvar) for treatment arm.
#' @param des.res the desired response. "larger": prefer larger response (default) "smaller": prefer smaller response
#' @param method.subgrp method for subgroup identification: "BATTing" or "Youden"
#' @param scale logical variable indicating if scaling is required
#' @param xselect it is a logical flag. TRUE: use LASSO for variable selection; FALSE: use all predictors to calculate the composite score. Default = TRUE
#' @param cv.iter.xselect it is the number of iterations of cross-validation used for xselect
#' @param k.fold.xselect it is the number of folds of cross-validation used for xselect
#' @param method.xselect it is the method used for variable selection: CV (default) cross validated AUC/ABC; "aAUC" approximate AUC (only for prognostic case)
#' @param pre.filter it is either "opt" or an integer between 1 and total number of predictors. Default is "opt", no prefiltering conducted, optimized number of predictors selected; An integer: min(opt, integer) of predictors selected
#' @param filter.method method for filtering: NULL, no prefiltering, "univariate", univaraite filtering; "glmnet", glmnet filtering, "unicart": univariate rpart filtering for prognostic case
#' @param n.boot number of bootstrap for batting or Youden procedure
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff.
#'
#' @return This function returns a list containing following entries
#' \item{optAUC}{output entries as returned by optAUCX function.}
#' \item{subgrp}{a list containing prediction rule based on optAUC score, pvalue and group summary statistics of each subgroup and interaction plot (only for predictive case)}
#'
#' @export
#' @examples
#' library(optAUCX)
#' library(MASS)
#' library(stringr)
#'
optAUCX.fit <- function(data, yvar, yvar.ref=1, xvars, censorvar=NULL, trtvar=NULL, type="b", trtref=1, des.res="larger", method.subgrp="Youden", scale=TRUE, xselect=TRUE, cv.iter.xselect=20, k.fold.xselect=5, method.xselect="CV", pre.filter="opt", filter.method=NULL, n.boot=10, min.sigp.prcnt=0.2) {
  res.optAUC <-optAUCX(outcome=yvar, predictor=xvars, censorvar=censorvar, treatment=trtvar, type=type, data=data, resp.ref=yvar.ref, trt.ref=trtref, xselect=xselect, scale=scale, n.boot=n.boot, cv.iter.xselect=cv.iter.xselect, k.fold.xselect=k.fold.xselect, method.xselect=method.xselect, pre.filter=pre.filter, filter.method=filter.method)
  if (method.subgrp=="BATTing"){
    res.subgrp <- BATTing.func(data=res.optAUC$data, yvar=yvar, xvar="score",censorvar=censorvar, type=type, trtvar = trtvar, trtref = trtref, des.res=des.res, min.sigp.prcnt = min.sigp.prcnt, n.boot=n.boot)
  } else if (method.subgrp=="Youden") {
    if (type=="b")
    {
      res.subgrp <- Youden.func(data=res.optAUC$data,yvar=yvar,xvar="score",trtvar=trtvar,trtref=trtref,censorvar=censorvar,type=type,n.boot=n.boot,des.res=des.res, min.sigp.prcnt=min.sigp.prcnt)
    } else {
      stop("Youden method only works for binary outcome")
    }
  } else {
    stop("subgroup method can either be 'BATTing' or 'Youden'")
  }
  res.subgrp$res<-list(beta=res.optAUC$beta,rule=res.subgrp$res)
  return(list(optAUC=res.optAUC,subgrp=res.subgrp))
}

#' The wrapper function for optAUCX, to be passed to kfold.cv.
#'
#' @param data data frame equal to cbind(y, x), where y and x are inputs to optAUCX.
#' @param args list containing all other input arguments to optAUCX. except for x and y. Also contains xvars=names(x) and yvar=names(y).
#'
#' @return prediction rule returned by optAUCX.
#'
optAUCX.wrapper <- function(data, args){

  # Unpack args list.
  for (name in names(args)) {
    assign(name, args[[name]])
  }

  res <- optAUCX.fit(data=data, yvar=yvar, yvar.ref=yvar.ref, xvars=xvars, censorvar=censorvar,
                     trtvar=trtvar, type=type, trtref=trtref, des.res=des.res, method.subgrp=method.subgrp,
                     scale=scale, xselect=xselect, cv.iter.xselect=cv.iter.xselect,
                     k.fold.xselect=k.fold.xselect, method.xselect=method.xselect, pre.filter=pre.filter,
                     filter.method=filter.method, n.boot=n.boot)

  return(res)
}

#' Prediction function for optAUC CV.
#' @description Predict outcome from optAUCX fit for cross-validation data given trained model in optAUCX.fit. This is intended to be used with kfold.cv.
#'
#' @param data data frame equal to cbind(y, x), where y and x are inputs to optAUCX.
#' @param predict.rule prediction rule returned by optAUCX.fit.
#' @param args prediction rule arguments.
#'
#' @return  Data frame containing two columns: y and pred.class, a logical vector indicating the prediction for each row of data.
pred.optAUCX.CV <- function(data, predict.rule, args){

  # Input arguments:
  #     data:
  #     res - optAUCX fit result lists returned by optAUCX.fit function.
  #     args - other arguments - unused for optAUCX.GREP.CV, but kfold.cv includes it.
  #
  # Output:

    select.x <- names(predict.rule$subgrp$res$beta)
    cutoff <- predict.rule$subgrp$res$rule[3]
    dir <- predict.rule$subgrp$res$rule[2]
    #   new.x <- scale(x=data[,select.x],center=predict.rule$mean[select.x], scale=predict.rule$sd[select.x])
    new.x<-as.matrix(data[,select.x])
    new.score <- new.x%*% matrix(predict.rule$subgrp$res$beta,ncol=1)
    if (dir=="<") {
      id <- new.score<cutoff
    } else {
      id <- new.score>cutoff
    }
    pred.id <- id

    pred.data <- cbind(data, data.frame(pred.class=pred.id))
    # Note!!: If you don't cbind pred.class with the data in the current fold, the
    # row numbers of the current fold will be lost.
    pred.data <- subset(pred.data, select="pred.class")
    # Now that the row numbers are attached, you can get rid of marker.data.
    pred.data
}


#' CV function for optAUCX.fit
#'
#' @param data the input data frame.
#' @param yvar response variable name.
#' @param yvar.ref response variable reference.
#' @param xvars vector of predictor variable names.
#' @param censorvar censoring variable name only for TTE data , censor=0,event=1)  - default = NULL.
#' @param trtvar treatment variable name (for predictive signature). Set trt.vec to NULL for prognostic signature.
#' @param type type of response variable: "s" survival; "b" binary (default).
#' @param trtref coding (in the column of trtvar) for treatment arm.
#' @param des.res the desired response. "larger": prefer larger response (default) "smaller": prefer smaller response.
#' @param method.subgrp method for subgroup identification: "BATTing" or "Youden".
#' @param scale logical variable indicating if scaling is required.
#' @param xselect a logical flag. TRUE: use LASSO for variable selection; FALSE: use all predictors to calculate the composite score. Default = TRUE.
#' @param cv.iter.xselect the number of iterations of cross-validation used for xselect.
#' @param k.fold.xselect the number of folds of cross-validation used for xselect
#' @param method.xselect the method used for variable selection: CV (default) cross validated AUC/ABC; "aAUC" approximate AUC (only for prognostic case).
#' @param pre.filter it is either "opt" or an integer between 1 and total number of predictors. Default is "opt", no prefiltering conducted, optimized number of predictors selected; An integer: min(opt, integer) of predictors selected
#' @param filter.method method for filtering: NULL, no prefiltering, "univariate", univaraite filtering; "glmnet", glmnet filtering, "unicart": univariate rpart filtering for prognostic case.
#' @param n.boot number of bootstrap for batting or Youden procedure.
#' @param cv.iter.eval the number of iterations for outer cross-validation.
#' @param k.fold.eval the number of folds for outer cross-validation.
#' @param max.iter Maximum number of successful iterations required for CV.
#'
#' @return a list containing with following entries:
#' \item{stats.summary}{Summary of performance statistics.}
#' \item{pred.classes}{Data frame containing the predictive clases (TRUE/FALSE) for each iteration.}
#' \item{folds}{Data frame containing the fold indices (index of the fold for each row) for each iteration.}
#' \item{sig.list}{List of length cv.iter * k.fold containing the signature generated at each of the  k folds, for all iterations.}
#' \item{error.log}{List of any error messages that are returned at an iteration.}
#' \item{interplot}{Treatment*subgroup interaction plot for predictive case}
#'
#' @export
#'
cv.optAUCX.fit <- function(data, yvar, yvar.ref=1, xvars, censorvar=NULL, trtvar=NULL, type="b", trtref=1, des.res="larger", method.subgrp="BATTing", scale=TRUE, xselect=TRUE, cv.iter.xselect=20, k.fold.xselect=5, method.xselect="CV", pre.filter="opt", filter.method=NULL, n.boot=10, cv.iter.eval=20,k.fold.eval=5,max.iter=500)
{
  if (type=="b") strata=data[,yvar]
  if (type=="s") strata=data[, censorvar]
  if (type=="c") strata=NULL

  model.Rfunc <- "optAUCX.wrapper"
  # Create a list containing arguments to GREP.fit
  model.Rfunc.args <- list(yvar=yvar, yvar.ref=yvar.ref, xvars=xvars, censorvar=censorvar, trtvar=trtvar, type=type, trtref=trtref, des.res=des.res, method.subgrp=method.subgrp, scale=scale, xselect=xselect, cv.iter.xselect=cv.iter.xselect, k.fold.xselect=k.fold.xselect, method.xselect=method.xselect, pre.filter=pre.filter, filter.method=filter.method, n.boot=n.boot)
  if(is.null(censorvar)){
    censor.vec <- NULL
  }else{
    censor.vec=data[,censorvar,drop=FALSE]
  }
  if(is.null(trtvar)){
    trt.vec <- NULL
  }else{
    trt.vec=data[,trtvar,drop=FALSE]
  }

  # List of arguments for model.Rfunc. Must be packaged into
  # a list to pass to kfold.cv
  predict.Rfunc <- "pred.optAUCX.CV"
  # List of arguments for predict.Rfunc.
  predict.Rfunc.args <- list(yvar=yvar, xvars=xvars)
  res <- kfold.cv(data=data, model.Rfunc=model.Rfunc, model.Rfunc.args=model.Rfunc.args, predict.Rfunc=predict.Rfunc, predict.Rfunc.args=predict.Rfunc.args, k.fold=k.fold.eval, cv.iter=cv.iter.eval, strata=strata, max.iter=max.iter)
  if (length(res) > 0) {
    stats <- evaluate.cv.results(cv.data=res$cv.data, y=data[,yvar,drop=FALSE], censor.vec=censor.vec, trt.vec=trt.vec, type=type)
    summary <- summarize.cv.stats(stats$raw.stats, trtvar, type)
    interplot=interaction.plot(data.eval=summary, type=type, main="Interaction Plot", trt.lab=c("Trt.", "Ctrl."))

    results <- list(stats.summary=summary, pred.classes=stats$pred.classes, folds=stats$folds, sig.list=res$sig.list, raw.stats=stats$raw.stats, error.log=res$error.log, interplot=interplot)
  } else {
    results <- "No successful cross-validations."
  }

  return(results)

}


#' Prediction function for optAUCX.fit
#'
#' @param x input predictors matrix.
#' @param res.optAUC.fit Prediction rule returned by optAUCX.fit.
#'
#' @return a logical vector indicating the prediction for each row of data.
#' @export
pred.optAUCX<-function(x,res.optAUC.fit)
{
  score<-as.matrix(x[,names(res.optAUC.fit$beta), drop=FALSE])%*%as.matrix(res.optAUC.fit$beta)
  data.score<-cbind(x,score)
  pred.data=pred.seqlr(data.score["score"], res.optAUC.fit$rule)
  return(pred.data)
}


