############################################################
# optAUC for prognostic and predictive signature
############################################################

##########################  Function part ##################


#' Function for AUC when input is X and Y.
#'
#' @param X nondiseased sample
#' @param Y diseased sample
#'
#' @return AUC for given dieased and non-diseased sample
#'
AUC.emp <- function(X, Y){
  mean(outer(Y,X, FUN=function(y,x) 1*(y>x)))
}


#' This function is for calculating C index.

#' @param yvar is the outcome variable name
#' @param score is the composite score
#' @param censorvar is the censoring variable name
#' @param data is the data frame used
#'
#' @return C index for the survival outcome.
#' @author Xin Huang, Yan Sun, Lu Tian, Saptarshi Chatterjee, Viswanath Devanarayan
#'
C.index <- function(yvar, score, censorvar, data){
  Zi <- data[, score, drop=FALSE]
  Zj <- data[, score, drop=FALSE]
  Ti <- data[, yvar, drop=FALSE]
  Tj <- data[, yvar, drop=FALSE]
  Cs <- data[, censorvar, drop=FALSE]
  Z<-(Zi[rep(1:nrow(Zi), each=nrow(Zj)), ] - Zj[rep(1:nrow(Zj), nrow(Zi)),])
  Ts <- (Ti[rep(1:nrow(Ti), each=nrow(Tj)), ] - Tj[rep(1:nrow(Tj), nrow(Ti)),])
  dels <- Cs[rep(1:nrow(Cs), each=nrow(Cs)), ]
  sum(1*(Z<0)*(Ts>0)*dels/sum(1*(Ts>0)*dels))
}

#' Internal function for hessAUC.
#'
#' @param z (m x n) x p data matrix as defined in betahat/betahatX function.
#' @param beta coefficient estimates from betahat/betahatX function.
#'
#' @return Hessian matrix components.
#'
hessAUC.sub <- function(z, beta){
  z%*%t(z)*c((exp(t(beta)%*%(z))/(1+exp(t(beta)%*%(z)))^2))
}

#' function for Hessian matrix of AUC.
#'
#' @param beta coefficient estimates from betahat/betahatX function.
#' @param Z (m x n) x p data matrix as defined in betahat/betahatX function.
#' @param w inverse probability weighting for missing data (default is 1).
#'
#' @return Hessian matrix of AUC.
#'
hessAUC <- function(beta, Z, w=1){
  #  Z <- Y[rep(1:nrow(Y), each=nrow(X)), ] - X[rep(1:nrow(X), nrow(Y)),]
  tmp <- w*apply(Z, MARGIN=1, FUN=hessAUC.sub, beta)
  matrix(rowMeans(tmp), ncol(Z), ncol(Z))
  }

#' function of grad_square in the GCV
#'
#' @param z  (m x n) x p data matrix as defined in betahat/betahatX function.
#' @param beta coefficient estimates from betahat/betahatX function.
#'
#' @return grad_square in the GCV
#'
grad.sub <- function(z, beta){
  (z)*(-1/(1+exp(t(beta)%*%(z))))
}

#' Internal function for HIC calculation
#'
#' @param beta coefficient estimates from betahat/betahatX function.
#' @param X m x p data matrix for m non-diseased subjects with p markers.
#' @param Y n x p data matrix for n diseased subjects with p markers.
#' @param Z0 (m x n) x p data matrix as defined in betahat/betahatX function.
#' @param w inverse probability weighting for missing data (default is 1).
#' @param weights how much weight each predictor should get.
#'
#' @return gradient square for the GCV.
#'
gradsqr <- function(beta, X, Y, Z0, w=1, weights=1){
  m<-nrow(X)
  n<-nrow(Y)

  #  Z0 <- Y[rep(1:nrow(Y), each=nrow(X)), ] - X[rep(1:nrow(X), nrow(Y)),]
  tmp <- w*apply(Z0, MARGIN=1, FUN=grad.sub, beta)
  ################################################
  index.x <- rep(1:nrow(X), nrow(Y))[weights==1]
  an0 <- cbind(t(tmp), index.x)
  an1 <- do.call("rbind", by(an0[,1:ncol(X)], an0[,"index.x"], FUN=colSums))
  an2 <- apply(an1, MARGIN=1, FUN=function(x) x%*%t(x))
  a1 <- matrix(rowSums(an2), ncol(X), ncol(X))/(m*n)^2

  ################################################

  index.y <- rep(1:nrow(Y), each=nrow(X))[weights==1]
  am0 <- cbind(t(tmp), index.y)
  am1 <- do.call("rbind", by(am0[,1:ncol(X)], am0[,"index.y"], FUN=colSums))
  am2 <- apply(am1, MARGIN=1, FUN=function(x) x%*%t(x))
  a2 <- matrix(rowSums(am2), ncol(X), ncol(X))/(m*n)^2

  #################################################
  tmp.a <- apply(tmp, MARGIN=2, FUN=function(x) x%*%t(x))
  a <- matrix(rowSums(tmp.a), ncol(X), ncol(X))/(m*n)^2

  #################################################

  a0 <- a1+a2-a
  a0
}

#' Function for HIC calculation.
#'
#' @param beta coefficient estimates from betahat/betahatX function.
#' @param X m x p data matrix for m non-diseased subjects with p markers.
#' @param Y m x p data matrix for m diseased subjects with p markers.
#' @param Z (m x n) x p data matrix as defined in betahat/betahatX function.
#' @param weights weight of the variables
#'
#' @return A numeric value with corresponding HIC.
#'
HIC <- function(beta, X, Y, Z, weights){
  inv<-try(ginv(hessAUC(beta, Z)),silent=T)
  if(class(inv)=="try-error")
  {
    return(NA)
  } else {
    return(sum(diag(inv%*%gradsqr(beta, X, Y, Z, weights=weights))))
  }
}



#' Internal function for selecting candidate predictors for composite score.
#'
#' @param data the pseudo dataset as used in betahat or betahatX function.
#' @param type type of response variable: "b" binary; "s" survival.
#' @param yvar simulated binary variable name used in betahat or betahatX function.
#' @param xvars vector of variable names for predictors (covariates) assuming data is provided.
#' @param data.orig the original dataset used for training.
#' @param yvar.orig the original outcome variable (column) name.
#' @param censorvar name of censoring variable (1: event; 0: censor), default = NULL.
#' @param trtref the coding (in the column of treatment) for treatment arm.
#' @param n.boot the number of bootstrap for the MCglmnet procedure.
#' @param cv.iter.xselect the number of iterations of cross-validation used for MCglmnet.
#' @param k.fold.xselect the number of folds required for cross-validated AUC or ABC.
#' @param trtvar treatment variable name.
#' @param method.xselect the method used for variable selection: CV (default) cross validated AUC or ABC; "aAUC" approximate AUC (Only for prognostic case).
#'
#' @return the variables selected for composite score.
#'
#' @export
#'
xselect.glmnet <- function(data,type,yvar,xvars,data.orig,yvar.orig,censorvar=NULL,trtvar=NULL,trtref=1,n.boot=50,cv.iter.xselect=20,k.fold.xselect=5, method.xselect="CV"){

  pmax=length(xvars)
  library(glmnet)
  family="binomial"

  var.sel.all=NULL
  # fit cv.glmnet to get all choices of nonzero estimates for different values of tuning parameter in LASSO
  # get the unique nonzero estimates excluding the intercept
  beta.select <- betahat.fold.func(
    yvar=yvar.orig,
    xvars=xvars,
    censorvar=censorvar,
    trtvar=trtvar,
    trtref=trtref,
    data=data.orig,
    type=type,
    pmax=pmax,
    family="binomial",
    nfolds = 5
  )

  if (!is.null(trtvar)){  ## predictive

    fit.cv.nzero.mat <- beta.select$betahat.mat

    if (method.xselect=="aABC") {
      aabc.cv.list<-NULL
      for (k in 1:dim(fit.cv.nzero.mat)[2])
      {
        xvars.nzero<-names(which(fit.cv.nzero.mat[,k]!=0))
        # cross validated aAUC for selected predictors
        aabc.cv<-aabc.cv.func(data=data.orig,yvar=yvar.orig,xvars=xvars.nzero,trtvar=trtvar,censorvar=censorvar,type=type)
        aabc.cv.list<-c(aabc.cv.list,aabc.cv)
      }
      aabcmax.idx<-which(aabc.cv.list==max(aabc.cv.list,na.rm=T))
      xvars.aabc<-xvars[which(fit.cv.nzero.mat[,aabcmax.idx]!=0)]
      return(xvars.aabc)
    } else if (method.xselect=="CV") {

      xvars.abc<-rownames(coef(beta.select$fit.cv, s = "lambda.1se"))[nonzeroCoef(coef(beta.select$fit.cv, s = "lambda.1se"))]
      return(xvars.abc)

    } else {
      stop("Treatment variable and xselect method do not match")
    }
  } else {  ## prognostic

    fit.cv.nzero.mat <- beta.select$betahat.mat

    if (method.xselect=="aAUC")
    {
      aauc.cv.list<-NULL
      for (k in 1:dim(fit.cv.nzero.mat)[2])
      {
        xvars.nzero<-names(which(fit.cv.nzero.mat[,k]!=0))
        # cross validated aAUC for selected predictors
        aauc.cv<-aauc.cv.func(data=data.orig,yvar=yvar.orig,xvars=xvars.nzero,censorvar=censorvar,type=type)
        aauc.cv.list<-c(aauc.cv.list,aauc.cv)
      }
      # predictors with maximum cross validated aAUC
      aaucmax.idx<-which(aauc.cv.list==max(aauc.cv.list,na.rm=T))
      xvars.aauc<-xvars[which(fit.cv.nzero.mat[,aaucmax.idx]!=0)]
      return(xvars.aauc)
    } else if (method.xselect=="CV") {

      xvars.auc<-rownames(coef(beta.select$fit.cv, s = "lambda.1se"))[nonzeroCoef(coef(beta.select$fit.cv, s = "lambda.1se"))]

      return(xvars.auc)
    } else {
      stop("Treatment variable and xselect method do not match")
    }
  }
}

### following functions are for optAUCX

#' optAUC estimates for prognostic case.
#'
#' @description Fit optAUC to the data in prognostic case. Estimates along with the model bias correction are returned.
#'
#' @param yvar outcome variable (column) name for response variable assuming data is provided.
#' @param xvars vector of variable names for predictors (covariates) assuming data is provided.
#' @param censorvar name of censoring variable (1: event; 0: censor), default = NULL.
#' @param data the data frame for training dataset.
#' @param type the type of response variable: "b" binary; "s" survival.
#' @param xselect a logical flag. TRUE: use Cross Validated AUC or aAUC for variable selection; FALSE: use all predictors to calculate the composite score. Default = FALSE
#' @param cv.iter.xselect the number of iterations of cross-validation used for MCglmnet.
#' @param k.fold.xselect the number of folds required for cross-validated AUC.
#' @param method.xselect the method used for variable selection: CV (default) cross validated AUC; "aAUC" approximate AUC.
#'
#' @return This function returns a list of the following:
#' \item{beta}{the estimated linear coefficients.}
#' \item{bias}{the HIC for the selected model.}
#'
betahat <- function(yvar, xvars, censorvar=NULL, data, type="b", xselect=FALSE, cv.iter.xselect=20, k.fold.xselect=5, method.xselect="CV"){
  if(type=="b"){
    X <- data[data[,yvar]==0, xvars, drop=FALSE]
    Y <- data[data[,yvar]==1, xvars, drop=FALSE]
    Z<-(Y[rep(1:nrow(Y), each=nrow(X)), ] - X[rep(1:nrow(X), nrow(Y)),])
    weights <- 1
  }else if(type=="s"){
    Zi <- data[, xvars, drop=FALSE]
    Zj <- data[, xvars, drop=FALSE]
    Ti <- data[, yvar, drop=FALSE]
    Tj <- data[, yvar, drop=FALSE]
    Cs <- data[, censorvar, drop=FALSE]
    Z<-(Zj[rep(1:nrow(Zj), nrow(Zi)),] - Zi[rep(1:nrow(Zi), each=nrow(Zj)), ])
    Ts <- (Ti[rep(1:nrow(Ti), each=nrow(Tj)), ] - Tj[rep(1:nrow(Tj), nrow(Ti)),])
    dels <- Cs[rep(1:nrow(Cs), each=nrow(Cs)), ]
    weights <- 1*(Ts>0)*dels
    Z <- Z[weights==1,]
  }else{
    stop("optAUCX only handles binary (type=b) and survial (type=s) outcomes. \n")
  }
  yb.sim<-rbinom(nrow(Z),1,0.5)
  data.new<-cbind(yb.sim,Z)

  Y.new<-data.new[yb.sim==1,]
  X.new<-data.new[yb.sim==0,]

  X.new[,-1]<--X.new[,-1]
  data.fin<-rbind(Y.new,X.new)
  data.fin <- as.data.frame(data.fin)

  if(xselect==TRUE){
    var.sel <- xselect.glmnet(data=data.fin, type=type, yvar="yb.sim", xvars=names(data.fin[-1]),censorvar=censorvar,cv.iter.xselect=cv.iter.xselect, k.fold.xselect=k.fold.xselect, method.xselect=method.xselect, data.orig=data, yvar.orig=yvar)
  }else{
    var.sel <- xvars
  }
  fit.formula<-as.formula(paste("yb.sim~",paste(var.sel,sep="",collapse="+"),"-1"))
  beta.proposal<-glm(fit.formula,data=data.fin,family=binomial)$coefficients
  #  beta.hat <- theta2beta(beta2theta(beta.proposal))
  beta.hat <- beta.proposal
  names(beta.hat) <- var.sel


  if (length(var.sel)>1)
  {
    if(type=="b"){
      correct <- try(HIC(beta=beta.hat, X=X[,var.sel], Y=Y[,var.sel], Z=Z[,var.sel], weights=weights),silent=T)
      if (class(correct)=="try-error") correct<-NA
    }else if(type=="s"){
      correct <- try(HIC(beta=beta.hat, X=Zi[,var.sel], Y=Zj[,var.sel], Z=Z[,var.sel], weights=weights),silent = T)
      if (class(correct)=="try-error") correct<-NA
    }
    list(beta.hat=beta.hat, bias=correct)
  } else {
    list(beta.hat=beta.hat, bias=0)
  }
}


#' optAUC estimates for predictive case.
#'
#' @description Fit optAUC to the data in predictive case. Estimates along with the model bias correction are returned.
#'
#' @param yvar outcome variable (column) name for response variable assuming data is provided.
#' @param xvars vector of variable names for predictors (covariates) assuming data is provided.
#' @param censorvar name of censoring variable (1: event; 0: censor), default = NULL.
#' @param trtref coding (in the column of treatment) for treatment arm.
#' @param data data frame for training dataset.
#' @param type type of response variable: "b" binary; "s" survival.
#' @param xselect a logical flag. TRUE: use Cross Validated ABC for variable selection; FALSE: use all predictors to calculate the composite score. Default = FALSE.
#' @param cv.iter.xselect the number of iterations of cross-validation used for MCglmnet.
#' @param k.fold.xselect the number of folds required for cross-validated ABC.
#' @param method.xselect the method used for variable selection: CV (default) cross validated ABC.
#' @param trtvar traetment variable name.
#'
#' @return This function returns a list of the following:
#' \item{beta}{the estimated linear coefficients.}
#' \item{bias}{the HIC for the selected model.}
#'
betahatX <- function(yvar, xvars, censorvar=NULL, trtvar, trtref, data, type="b", xselect=FALSE, cv.iter.xselect=20, k.fold.xselect=5, method.xselect="CV"){
  #  X.trt: non-responder in treatment arm
  #  Y.trt: responder in treatment arm
  #  X.ctl: non-responder in control arm
  #  Y.ctl: responder in control arm

  if(type=="b"){
    X.trt <- data[data[,trtvar]==1 & data[,yvar]==0, xvars, drop=FALSE]
    Y.trt <- data[data[,trtvar]==1 & data[,yvar]==1, xvars, drop=FALSE]
    X.ctl <- data[data[,trtvar]==0 & data[,yvar]==0, xvars, drop=FALSE]
    Y.ctl <- data[data[,trtvar]==0 & data[,yvar]==1, xvars, drop=FALSE]

    Z.trt <- Y.trt[rep(1:nrow(Y.trt), each=nrow(X.trt)), ] - X.trt[rep(1:nrow(X.trt), nrow(Y.trt)),]
    Z.ctl <- X.ctl[rep(1:nrow(X.ctl), nrow(Y.ctl)),] - Y.ctl[rep(1:nrow(Y.ctl), each=nrow(X.ctl)), ]
    Z <- rbind(Z.trt, Z.ctl)
    weights<-1
  }else if(type=="s"){
    Zi.trt <- data[data[,trtvar]==1, xvars, drop=FALSE]
    Zj.trt <- data[data[,trtvar]==1, xvars, drop=FALSE]
    Ti.trt <- data[data[,trtvar]==1, yvar, drop=FALSE]
    Tj.trt <- data[data[,trtvar]==1, yvar, drop=FALSE]
    Cs.trt <- data[data[,trtvar]==1, censorvar, drop=FALSE]
    Z.trt<-(Zj.trt[rep(1:nrow(Zj.trt), nrow(Zi.trt)),] - Zi.trt[rep(1:nrow(Zi.trt), each=nrow(Zj.trt)), ])
    Ts.trt <- (Ti.trt[rep(1:nrow(Ti.trt), each=nrow(Tj.trt)), ] - Tj.trt[rep(1:nrow(Tj.trt), nrow(Ti.trt)),])
    dels.trt <- Cs.trt[rep(1:nrow(Cs.trt), each=nrow(Cs.trt)), ]
    weights.trt <- 1*(Ts.trt>0)*dels.trt
    Z.trt <- Z.trt[weights.trt==1,]

    Zi.ctl <- data[data[,trtvar]==0, xvars, drop=FALSE]
    Zj.ctl <- data[data[,trtvar]==0, xvars, drop=FALSE]
    Ti.ctl <- data[data[,trtvar]==0, yvar, drop=FALSE]
    Tj.ctl <- data[data[,trtvar]==0, yvar, drop=FALSE]
    Cs.ctl <- data[data[,trtvar]==0, censorvar, drop=FALSE]
    Z.ctl<-(Zi.ctl[rep(1:nrow(Zi.ctl), each=nrow(Zj.ctl)), ] - Zj.ctl[rep(1:nrow(Zj.ctl), nrow(Zi.ctl)),])
    Ts.ctl <- (Ti.ctl[rep(1:nrow(Ti.ctl), each=nrow(Tj.ctl)), ] - Tj.ctl[rep(1:nrow(Tj.ctl), nrow(Ti.ctl)),])
    dels.ctl <- Cs.ctl[rep(1:nrow(Cs.ctl), each=nrow(Cs.ctl)), ]
    weights.ctl <- 1*(Ts.ctl>0)*dels.ctl
    Z.ctl <- Z.ctl[weights.ctl==1,]
    Z <- rbind(Z.trt, Z.ctl)
  }else{
    stop("optAUCX only handles binary (type=b) and survial (type=s) outcomes. \n")
  }

  yb.sim<-rbinom(nrow(Z),1,0.3)
  data.new<-cbind(yb.sim,Z)
  Y.new<-data.new[yb.sim==1,]
  X.new<-data.new[yb.sim==0,]
  X.new[,-1]<--X.new[,-1]
  data.fin<-rbind(Y.new,X.new)
  ##
  data.fin<-cbind(data.fin,c(rep(1,nrow(Z.trt)),rep(0,nrow(Z.ctl))))
  names(data.fin)[ncol(data.fin)]<-trtvar
  ##
  data.fin <- as.data.frame(data.fin)

  if(xselect==TRUE){
    var.sel <- xselect.glmnet(data=data.fin, type=type, yvar="yb.sim", xvars=xvars, trtvar = trtvar, trtref = trtref, censorvar=censorvar, cv.iter.xselect=cv.iter.xselect, k.fold.xselect=k.fold.xselect, method.xselect=method.xselect,data.orig = data, yvar.orig = yvar)
  }else{
    var.sel <- xvars
  }

  #

  fit.formula<-as.formula(paste("yb.sim~",paste(var.sel,sep="",collapse="+"),"-1"))
  beta.proposal<-glm(fit.formula,data=data.fin,family=binomial)$coefficients
  #  beta.hat <- theta2beta(beta2theta(beta.proposal))
  beta.hat <- beta.proposal
  names(beta.hat) <- var.sel

#   if(length(var.sel)>1)
#   {
#     correct<-try(abs(HIC(beta=beta.hat, X=X.trt[,var.sel], Y=Y.trt[,var.sel], Z=Z.trt[,var.sel], weights=weights))+abs(HIC(beta=beta.hat, X=X.ctl[,var.sel], Y=Y.ctl[,var.sel], Z=Z.ctl[,var.sel], weights=weights)),silent = T)
#     if (sum(beta.hat!=0)>5) ;
#     if (class(correct)=="try-error") ; stop;
#   }

  #####
  if (length(var.sel)>1)
  {
    if(type=="b"){
      correct <- try(abs(HIC(beta=beta.hat, X=X.trt[,var.sel], Y=Y.trt[,var.sel], Z=Z.trt[,var.sel], weights=weights)) +
        abs(HIC(beta=beta.hat, X=X.ctl[,var.sel], Y=Y.ctl[,var.sel], Z=Z.ctl[,var.sel], weights=weights)),silent=T)

      if (class(correct)=="try-error") correct<-NA
    }else if(type=="s"){
      correct <- try(abs(HIC(beta=beta.hat, X=Zi.trt[,var.sel], Y=Zj.trt[,var.sel], Z=Z.trt[,var.sel], weights=weights.trt)) +
        abs(HIC(beta=beta.hat, X=Zi.ctl[,var.sel], Y=Zj.ctl[,var.sel], Z=Z.ctl[,var.sel], weights=weights.ctl)),silent=T)
      if (class(correct)=="try-error") correct<-NA
    }
    list(beta.hat=beta.hat, bias=correct)
  } else {
    list(beta.hat=beta.hat, bias=0)
  }
}


#' Optimal Combination of biomarkers for prognostic and predictive signature development Based on total AUC.
#' @description Searching for optimal linear combination of multiple biomarkers that (1) for prognostic signature development, maximizes the total area under the receiver operating characteristic (ROC) curves (for binary outcome)/C-index (for survival outcome); (2) for predictive signature development, maximizes the area between ROC curves (for binary outcome)/C-index (for survival outcome) between treatment and control arm.
#'
#' @param outcome variable (column) name for response variable if data is provided; otherwise, a vector of responses.
#' @param predictor vector of variable names for predictors (covariates) if data is provided; otherwise, a data matrix of predictors.
#' @param censorvar name of censoring variable (1: event; 0: censor), default = NULL.
#' @param treatment variable name for treatment variable if data is provided; otherwise, a vector of treatment variables. default = NULL (prognostic signature).
#' @param type type of response variable: "b" binary; "s" survival.
#' @param data data frame for training dataset.
#' @param resp.ref coding (in the column of outcome) for responder, used for binary outcome.
#' @param trt.ref coding (in the column of treatment) for treatment arm.
#' @param xselect a logical flag. TRUE: use LASSO for variable selection; FALSE: use all predictors to calculate the composite score. Default = FALSE.
#' @param scale a logical flag indicating whether standardization is to be performed to the dataset before the combination, default is TRUE.
#' @param n.boot the number of bootstrap for the MCglmnet procedure.
#' @param cv.iter.xselect the number of iterations for internal cross-validation used for xselect.
#' @param k.fold.xselect the number of folds for internal cross-validation used for xselect.
#' @param method.xselect the method used for variable selection: "CV" (default) cross validated AUC/ABC; "aAUC" approximate AUC (only for prognostic case).
#' @param pre.filter NULL, no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected.
#' @param filter.method NULL, no prefiltering, "univariate", univaraite filtering; "glmnet", glmnet filtering, "unicart": univariate rpart filtering for prognostic case.
#'
#' @return a list containing following entries:
#' \item{beta}{linear coefficient estimates.}
#' \item{data}{the dataset with composite score.}
#' \item{AUC}{the apparent AUC (C-index if survival outcome) calculated from the estimated beta (prognostic case only).}
#' \item{trt.AUC}{the apparent AUC (C-index if survival outcome) calculated from the estimated beta for the treatment arm (predictive case only).}
#' \item{ctl.AUC}{the apparent AUC (C-index if survival outcome) calculated from the estimated beta for the control arm (predictive case only).}
#' \item{ABC}{the area between the apparent ROCs (C-index if survival outcome) calculated from the estimated beta for the control arm (predictive case only).}
#'
#' \item{sd}{a numeric vector of standard deviations of the columns of the training predictors.}
#' \item{mean}{a numeric vector of mean of the columns of the training predictors.}
#' \item{aAUC}{the approximate AUC correcting the bias of the model (prognostic case only).}
#' \item{aABC}{the approximate ABC correcting the bias of the model (predictive case only).}
#'
optAUCX <- function(outcome, predictor, censorvar=NULL, treatment=NULL, type="b", data=NULL, resp.ref=1, trt.ref=1, xselect=FALSE, scale=TRUE, n.boot=10, cv.iter.xselect=20, k.fold.xselect=5, method.xselect="CV", pre.filter=NULL, filter.method=NULL){
  if(class(outcome) == "character" &  class(predictor) == "character"){
    all.data <- subset(data, select=predictor)
    resp <- subset(data, select=outcome)
    if(type=="b"){
      resp <- 1*(resp==resp.ref)
    }
    if(!is.null(treatment)){
      trt <- subset(data, select=treatment)
      trt <- 1*(trt==trt.ref)
    }else{
      trt <- NULL
    }
    if(!is.null(censorvar)){
      censor <- subset(data, select=censorvar)
    }else{
      censor <- NULL
    }
    pred.len <- length(predictor)
  }else if (is.null(data)){
    all.data <- as.matrix(predictor)
    resp <- as.matrix(outcome)
    colnames(resp) <- "outcome"
    if(!is.null(censorvar)){
      censor <- as.matrix(censorvar)
      colnames(censor) <- "censor"
      censorvar <- "censor"
      if(!is.null(treatment)){
        trt <- as.matrix(treatment)
        colnames(trt) <- "treatment"
        data <- data.frame(resp, censor, trt, all.data)
      }else{
        trt <- NULL
        data <- data.frame(resp, censor, all.data)
      }
    }else{
      if(!is.null(treatment)){
        trt <- as.matrix(treatment)
        colnames(trt) <- "treatment"
        data <- data.frame(resp, trt, all.data)
      }else{
        trt <- NULL
        data <- data.frame(resp, all.data)
      }
    }
    pred.len <- ncol(all.data)
  }else{
    stop("if data is a given data frame, input columne names for outcome, treatment and predictors; if data is NULL, input data frame/matrix for outcome, treatment, predictor. \n")
  }

  if(scale==TRUE){
    sd.data <- apply(all.data, MARGIN=2, FUN=sd)
    mean.data <- colMeans(all.data)
    all.data.orig<-all.data # original all.data will be used for rescaling betahat and scores
    all.data <- scale(all.data)
  }else{
    sd.data <- NULL
    cov.data <- NULL
    mean.data <- NULL
  }

  #### needs to be checked ###
  if (!is.null(pre.filter)){
    xvars_new=filter(data=data,type=type,yvar=outcome,xvars=predictor,censorvar=censorvar,trtvar=treatment,n.boot=n.boot,cv.iter=15,pre.filter=pre.filter,filter.method=filter.method)
    xvars=xvars_new
    all.data<-all.data[,xvars]
  }
  ###############################


  if(is.null(treatment)){
    if(is.null(censorvar)){
      data.in <- data.frame(resp, all.data)
    }else{
      data.in <- data.frame(resp, censor, all.data)
    }

    result<-betahat(yvar=colnames(resp), xvars=colnames(all.data), censorvar = censorvar, type = type, data=data.in, xselect=xselect, cv.iter.xselect=cv.iter.xselect, k.fold.xselect=k.fold.xselect, method.xselect=method.xselect)

    if (scale) # in case scale=T, betahat and score need to be rescaled
    {
      result$beta.hat<-result$beta.hat/sd.data[names(result$beta.hat)]
      score.optAUCX <- t(result$beta.hat)%*%t(all.data.orig[,names(result$beta.hat), drop=FALSE])
    } else {
      score.optAUCX <- t(result$beta.hat)%*%t(all.data[,names(result$beta.hat), drop=FALSE])
    }

    data.score <- data.frame(cbind(data, t(score.optAUCX)))
    colnames(data.score)[ncol(data.score)] <- "score"
    if(type=="b"){
      X.score <- data.score[data.score[,outcome]==0, "score"]
      Y.score <- data.score[data.score[,outcome]==1, "score"]
      AUC <- AUC.emp(X.score, Y.score)
    }else if(type=="s"){
      AUC <- C.index(yvar=colnames(resp), score="score", censorvar=censorvar, data=data.score)
    }

    list(beta=result$beta.hat,
         data=data.score, AUC=AUC, sd=sd.data, mean=mean.data, aAUC=AUC-abs(result$bias))

  }else{
    if(is.null(censorvar)){
      data.in <- data.frame(resp, trt, all.data)
    }else{
      data.in <- data.frame(resp, censor, trt, all.data)
    }
    result <- betahatX(yvar=colnames(resp), xvars=colnames(all.data), censorvar=censorvar, trtvar=colnames(trt), trtref=trt.ref, type = type, data=data.in, xselect=xselect, cv.iter.xselect=cv.iter.xselect, k.fold.xselect=k.fold.xselect, method.xselect=method.xselect)

    if (scale) # in case scale=T, betahat and score need to be rescaled
    {
      result$beta.hat<-result$beta.hat/sd.data[names(result$beta.hat)]
      score.optAUCX <- t(result$beta.hat)%*%t(all.data.orig[,names(result$beta.hat), drop=FALSE])
    } else {
      score.optAUCX <- t(result$beta.hat)%*%t(all.data[,names(result$beta.hat), drop=FALSE])
    }

    data.score <- data.frame(cbind(data, t(score.optAUCX)))
    colnames(data.score)[ncol(data.score)] <- "score"

    if(type=="b"){
      X.trt.score <- data.score[data.score[,treatment]==1 & data.score[,outcome]==0, "score"]
      Y.trt.score <- data.score[data.score[,treatment]==1 & data.score[,outcome]==1, "score"]
      X.ctl.score <- data.score[data.score[,treatment]==0 & data.score[,outcome]==0, "score"]
      Y.ctl.score <- data.score[data.score[,treatment]==0 & data.score[,outcome]==1, "score"]

      trt.AUC <- AUC.emp(X.trt.score, Y.trt.score)
      ctl.AUC <- AUC.emp(X.ctl.score, Y.ctl.score)
      ABC <- trt.AUC - ctl.AUC
    }else if(type=="s"){
      trt.AUC <- C.index(yvar=colnames(resp), score="score", censorvar=censorvar, data=data.score[data.score[,treatment]==1,])
      ctl.AUC <- C.index(yvar=colnames(resp), score="score", censorvar=censorvar, data=data.score[data.score[,treatment]==0,])
      ABC <- trt.AUC - ctl.AUC
    }

    if (scale)
    {
      result$beta.hat
    }

    list(beta=result$beta.hat,
         data=data.score, trt.AUC=trt.AUC, ctl.AUC=ctl.AUC,
         ABC=ABC, sd=sd.data, mean=mean.data, aABC=ABC-abs(result$bias))
  }
}

#' Internal function used in xselect.glmnet.
#'
#' @param m the matrix of median and se from cross-validated AUC.
#'
#' @return row number corresponding to optimum AUC/ABC according to 1se rule.
#'
opt.idx.func<-function(m)
{
  #
  m<-data.frame(m)
  names(m)<-c("med","se")
  medmax<-max(m$med)
  semax<-m[which.max(m$med),"se"]
  if (sum(m$med<=(medmax-semax))==0) # in case, there is no improvement of AUC or ABC as we increase number of predictors
  {
    opt<-1
    return(opt)
  } else {
    opt<-max(m$med[m$med<=(medmax-semax)])
    return(min(which(m$med==opt))) # min is added in case there are two optimal auc values it will choose the minimum one
  }
}


#' Internal function used in xselect.glmnet.
#'
#' @param data the original dataset.
#' @param yvar the original outcome variable name.
#' @param xvars the name of predictors.
#' @param censorvar the name of censor variable.
#' @param type "b" for binary outcome or "s" for time to event outcome.
#'
#' @return cross validated aAUC for selected predictors.
#'
aauc.cv.func<-function(data,yvar,xvars,censorvar,type)
{
  if (length(xvars)==1)
  { # for one predictor no need to run optAUCx. The predictors are composite scores themselves. AUC is approx AUC as there is no bias involved.
    if (type=="b")
    {
      X.score<-data[data[,yvar]==0,xvars]
      Y.score<-data[data[,yvar]==1,xvars]
      aauc<-AUC.emp(X.score,Y.score)
    } else if (type=="s") {
      aauc<-C.index(yvar=yvar, score=xvars, censorvar=censorvar, data=data)
    }
  } else {
    #
    res<-optAUCX(outcome=yvar, predictor=xvars,censorvar=censorvar,
                 type=type, data=data, scale=F)
    aauc<-res$aAUC
  }
  return(aauc)
}

#' Internal function used in xselect.glmnet.
#'
#' @param data the original dataset.
#' @param yvar the original outcome variable name.
#' @param xvars the name of predictors.
#' @param censorvar the name of censor variable.
#' @param type "b" for binary outcome or "s" for time to event outcome.
#' @param trtvar treatment variable name.
#'
#' @return cross validated aABC for selected predictors.
#'
aabc.cv.func<-function(data,yvar,xvars,trtvar,censorvar,type)
{
  if (length(xvars)==1)
  { # for one predictor no need to run optAUCx. The predictors are composite scores themselves. ABC is approx ABC as there is no bias involved.
    if (type=="b")
    {
      X.trt.score <- data[data[,trtvar]==1 & data[,yvar]==0, xvars]
      Y.trt.score <- data[data[,trtvar]==1 & data[,yvar]==1, xvars]
      X.ctl.score <- data[data[,trtvar]==0 & data[,yvar]==0, xvars]
      Y.ctl.score <- data[data[,trtvar]==0 & data[,yvar]==1, xvars]

      trt.AUC <- AUC.emp(X.trt.score, Y.trt.score)
      ctl.AUC <- AUC.emp(X.ctl.score, Y.ctl.score)
      aabc <- trt.AUC - ctl.AUC
    } else if (type=="s") {
      trt.AUC <- C.index(yvar=yvar, score=xvars, censorvar=censorvar, data=data[data[,trtvar]==1,])
      ctl.AUC <- C.index(yvar=yvar, score=xvars, censorvar=censorvar, data=data[data[,trtvar]==0,])
      aabc <- trt.AUC - ctl.AUC
    }
  } else {
    #
    res<-optAUCX(outcome=yvar, predictor=xvars, treatment =trtvar, censorvar=censorvar, type=type, data=data, scale=F)
    aabc<-res$aABC
  }
  return(aabc)
}



#' Function for internal CV to estimate number of candidate predictors required for composite score.
#'
#' @param data the dataset used for training.
#' @param trtref the coding (in the column of treatment) for treatment arm.
#' @param yvar the original outcome variable name.
#' @param xvars the name of predictors.
#' @param censorvar the name of censor variable.
#' @param type "b" for binary outcome or "s" for time to event outcome.
#' @param cv.iter.xselect number of iterations for internal cross-validation.
#' @param k.fold.xselect number of folds for internal cross-validation.
#'
#' @return returns median and se of cross-validated AUC/ABC and across internal cv iterations (required for 1se rule).
#'
cv.xselect<-function(data, yvar, xvars, censorvar, trtvar, type,trtref, cv.iter.xselect, k.fold.xselect, pmax, family)
{
  if (type=="b") strata=data[,yvar]
  if (type=="s") strata=data[, censorvar]

  n.data <- nrow(data)
  data.index <-  1:n.data

  if(is.null(strata)){
    strata <- rep(1, n.data)
  }

  iter.count<-1
  iter.success.count<-0

  aucabc.cv.mat<-rep(NULL,length(xvars))
  for (i in 1:cv.iter.xselect)
  {
    cat("\n\nCV iteration(xselect) ", iter.count, "(Successes: ", iter.success.count, ")\n\n")
    cv.index <- balanced.folds(strata, k.fold.xselect)
    train.index <- lapply(cv.index, function(x) setdiff(data.index, x))
    # Indices for training data for each of the k folds.
    cv.mat <- rep(NULL,length(xvars))
    fold.success.count <- 0
    # Number of successful fold evaluations.


    for (j in 1:k.fold.xselect) {
      cat("Fold(xselect) ", j)
      betahat.fold<-betahat.fold.func(yvar=yvar, xvars=xvars, censorvar=censorvar, trtvar=trtvar, trtref=trtref, data=data[train.index[[j]],], type=type, pmax=pmax, family=family)

      # Apply optAUCx to the training data to get a prediction rule.
      x=data[cv.index[[j]],xvars]
      #       score <- x%*% matrix(res$beta,ncol=1)
      score<-as.matrix(x)%*%betahat.fold
      # Get Xbeta for the left-out fold of data.
      fold.vec <- rep(j, nrow(score))
      # Record fold number for current subset of data.
      score <- cbind(score, data.frame("fold"=fold.vec))
      cv.mat <- rbind(cv.mat, score)
      # Append the output of prediction function to cv.vec.
      fold.success.count <- fold.success.count + 1
      cat("\n")
    }
    cv.mat$row.numbers <- as.numeric(rownames(cv.mat))
    cv.mat <- cv.mat[order(cv.mat$row.numbers), ] #order row numbers
    cv.mat <- subset(cv.mat, select=setdiff(names(cv.mat), c("fold","row.numbers"))) #exclude row numbers column
    aucabc.cv.vec<-apply(cv.mat,2,function(x) aucabc.func(score=x,data=data,yvar=yvar,xvars=xvars,trtvar=trtvar,censorvar=censorvar,type=type)) # get AUC or ABC for each column in cv.mat

    aucabc.cv.mat<-rbind(aucabc.cv.mat,aucabc.cv.vec)
    iter.success.count<-iter.success.count+1
    iter.count<-iter.count+1
  }

  aucabc.cv.med<-apply(aucabc.cv.mat,2,median)
  aucabc.cv.1se<-apply(aucabc.cv.mat,2,sd)
  return(cbind.data.frame(median=aucabc.cv.med,se=aucabc.cv.1se))
}


#' Function for AUC/ABC used in cv.xselect.
#'
#' @param score composite score from cv.xselect.
#' @param data the training dataset.
#' @param type type of response.
#' @param trtvar treatment variable name.
#' @param censorvar censoring variable name.
#' @param yvar response variable name.
#' @param xvars a vector of predictor names.
#'
#' @return AUC/ABC for each score from cv.xselect.
#'
aucabc.func<-function(score,data,type,trtvar,censorvar,yvar,xvars)
{
  cv.data<-cbind.data.frame(data,"score"=score)
  if (is.null(trtvar))
  {
    if (type=="b")
    {
      X.score<-cv.data[cv.data[,yvar]==0,"score"]
      Y.score<-cv.data[cv.data[,yvar]==1,"score"]
      auc.cv<-AUC.emp(X.score,Y.score)
    } else if (type=="s") {
      auc.cv<-C.index(yvar=yvar, score="score", censorvar=censorvar, data=cv.data)
    }
    return(auc.cv)
  } else {
    if (type=="b")
    {
      X.trt.score <- cv.data[cv.data[,trtvar]==1 & cv.data[,yvar]==0, "score"]
      Y.trt.score <- cv.data[cv.data[,trtvar]==1 & cv.data[,yvar]==1, "score"]
      X.ctl.score <- cv.data[cv.data[,trtvar]==0 & cv.data[,yvar]==0, "score"]
      Y.ctl.score <- cv.data[cv.data[,trtvar]==0 & cv.data[,yvar]==1, "score"]

      trt.AUC <- AUC.emp(X.trt.score, Y.trt.score)
      ctl.AUC <- AUC.emp(X.ctl.score, Y.ctl.score)
      abc.cv <- trt.AUC - ctl.AUC
    }else if(type=="s"){
      trt.AUC <- C.index(yvar=yvar, score="score", censorvar=censorvar, data=cv.data[cv.data[,trtvar]==1,])
      ctl.AUC <- C.index(yvar=yvar, score="score", censorvar=censorvar, data=cv.data[cv.data[,trtvar]==0,])
      abc.cv <- trt.AUC - ctl.AUC
    }
    return(abc.cv)
  }
}


#' beta estimates for each candidate predictor set.
#'
#' @param yvar response variable name.
#' @param xvars a vector of predictor variable names.
#' @param censorvar censoring variable name.
#' @param trtvar treatment variable name.
#' @param trtref treatment reference code.
#' @param data data frame for train data.
#' @param type "b" binary; "s" survival.
#' @param pmax maximum number of variables ever to be nonzero as required in glmnet function.
#' @param family response type.
#' @param nfolds number of CV folder used in cv.glmnet.
#'
#' @return the coefficient estimates for each fold in internal CV.
#'
#' @export
#'
betahat.fold.func<-function(yvar, xvars, censorvar, trtvar, trtref, data, type, pmax, family, nfolds = 5)
{
  if (is.null(trtvar))
  {
    if(type=="b"){
      X <- data[data[,yvar]==0, xvars, drop=FALSE]
      Y <- data[data[,yvar]==1, xvars, drop=FALSE]
      Y.idx <- rep(1:nrow(Y), each=nrow(X))
      X.idx <- rep(1:nrow(X), nrow(Y))
      Z<-(Y[Y.idx, ] - X[X.idx,])
      # set CV folds
      fold.X <- sample.int(nfolds, nrow(X), replace = TRUE)
      fold.Y <- sample.int(nfolds, nrow(Y), replace = TRUE)
      fold.Z <- ifelse(fold.X[X.idx] == fold.Y[Y.idx],
                       fold.X[X.idx], NA)
      weights <- 1
    }else if(type=="s"){
      Zi <- data[, xvars, drop=FALSE]
      Zj <- data[, xvars, drop=FALSE]
      Ti <- data[, yvar, drop=FALSE]
      Tj <- data[, yvar, drop=FALSE]
      Cs <- data[, censorvar, drop=FALSE]
      i.idx <- rep(1:nrow(Zi), each=nrow(Zj))
      j.idx <- rep(1:nrow(Zj), nrow(Zi))
      Z<-(Zi[i.idx, ] - Zj[j.idx,])
      Ts <- (Ti[i.idx, ] - Tj[j.idx,])
      dels <- Cs[i.idx, ]
      # set CV folds
      fold.ij <- sample.int(nfolds, nrow(data), replace = TRUE)
      fold.Z <- ifelse(fold.ij[i.idx] == fold.ij[j.idx],
                       fold.ij[i.idx], NA)
      weights <- 1*(Ts<0)*dels
      Z <- Z[weights==1,]
      fold.Z <- fold.Z[weights==1]
    }else{
      stop("optAUCX only handles binary (type=b) and survial (type=s) outcomes. \n")
    }
    yb.sim<-rbinom(nrow(Z),1,0.5)
    data.new<-cbind(yb.sim,Z)

    Y.new<-data.new[yb.sim==1,]
    X.new<-data.new[yb.sim==0,]
    yvar.new<-"yb.sim"
    X.new[,-1]<--X.new[,-1]
    data.fin<-rbind.data.frame(Y.new,X.new)
    data.fin <- data.frame(data.fin, cvfolds=fold.Z)

    x<-as.matrix(data.fin[,xvars])
    y<-data.fin[,yvar.new]

    # fit cv.glmnet to get all choices of nonzero estimates for different values of tuning parameter in LASSO
    cv.idx <- !is.na(data.fin$cvfolds)
    fit.cv=cv.glmnet(x[cv.idx,], y[cv.idx], family="binomial", nfolds=5, foldid=data.fin$cvfolds[cv.idx], pmax=pmax, intercept = FALSE, type.measure="auc", parallel = TRUE)
  } else {
    if(type=="b"){
      X.trt <- data[data[,trtvar]==1 & data[,yvar]==0, xvars, drop=FALSE]
      Y.trt <- data[data[,trtvar]==1 & data[,yvar]==1, xvars, drop=FALSE]
      X.ctl <- data[data[,trtvar]==0 & data[,yvar]==0, xvars, drop=FALSE]
      Y.ctl <- data[data[,trtvar]==0 & data[,yvar]==1, xvars, drop=FALSE]

      wt.trt <- (nrow(X.trt)*nrow(Y.trt)+nrow(X.ctl)*nrow(Y.ctl))/nrow(X.trt)*nrow(Y.trt)
      wt.ctl <- (nrow(X.trt)*nrow(Y.trt)+nrow(X.ctl)*nrow(Y.ctl))/nrow(X.ctl)*nrow(Y.ctl)

      X.trt.idx <- rep(1:nrow(X.trt), nrow(Y.trt))
      Y.trt.idx <- rep(1:nrow(Y.trt), each=nrow(X.trt))
      X.ctl.idx <- rep(1:nrow(X.ctl), nrow(Y.ctl))
      Y.ctl.idx <- rep(1:nrow(Y.ctl), each=nrow(X.ctl))

      Z.trt <- Y.trt[Y.trt.idx, ] - X.trt[X.trt.idx,]
      Z.ctl <- X.ctl[X.ctl.idx,] - Y.ctl[Y.ctl.idx, ]
      Z <- rbind(Z.trt, Z.ctl)
      wt.Z <- c(rep(wt.trt, nrow(Z.trt)), rep(wt.ctl, nrow(Z.ctl)))

      # set CV folds
      fold.X.trt <- sample.int(nfolds, nrow(X.trt), replace = TRUE)
      fold.Y.trt <- sample.int(nfolds, nrow(Y.trt), replace = TRUE)
      fold.Z.trt <- ifelse(fold.X.trt[X.trt.idx] == fold.Y.trt[Y.trt.idx],
                           fold.X.trt[X.trt.idx], NA)
      fold.X.ctl <- sample.int(nfolds, nrow(X.ctl), replace = TRUE)
      fold.Y.ctl <- sample.int(nfolds, nrow(Y.ctl), replace = TRUE)
      fold.Z.ctl <- ifelse(fold.X.ctl[X.ctl.idx] == fold.Y.ctl[Y.ctl.idx],
                           fold.X.ctl[X.ctl.idx], NA)
      fold.Z <- c(fold.Z.trt, fold.Z.ctl)
    }else if(type=="s"){
      Zi.trt <- data[data[,trtvar]==1, xvars, drop=FALSE]
      Zj.trt <- data[data[,trtvar]==1, xvars, drop=FALSE]
      Ti.trt <- data[data[,trtvar]==1, yvar, drop=FALSE]
      Tj.trt <- data[data[,trtvar]==1, yvar, drop=FALSE]
      Cs.trt <- data[data[,trtvar]==1, censorvar, drop=FALSE]

      i.idx.trt <- rep(1:nrow(Zi.trt), each=nrow(Zj.trt))
      j.idx.trt <- rep(1:nrow(Zj.trt), nrow(Zi.trt))

      Z.trt<-(Zi.trt[i.idx.trt, ] - Zj.trt[j.idx.trt,])
      Ts.trt <- (Ti.trt[i.idx.trt, ] - Tj.trt[j.idx.trt,])
      dels.trt <- Cs.trt[i.idx.trt, ]
      # set CV folds
      fold.ij.trt <- sample.int(nfolds, nrow(Zi.trt), replace = TRUE)
      fold.Z.trt <- ifelse(fold.ij.trt[i.idx.trt] == fold.ij.trt[j.idx.trt],
                       fold.ij.trt[i.idx.trt], NA)

      weights.trt <- 1*(Ts.trt<0)*dels.trt
      Z.trt <- Z.trt[weights.trt==1,]
      fold.Z.trt <- fold.Z.trt[weights.trt==1]

      Zi.ctl <- data[data[,trtvar]==0, xvars, drop=FALSE]
      Zj.ctl <- data[data[,trtvar]==0, xvars, drop=FALSE]
      Ti.ctl <- data[data[,trtvar]==0, yvar, drop=FALSE]
      Tj.ctl <- data[data[,trtvar]==0, yvar, drop=FALSE]
      Cs.ctl <- data[data[,trtvar]==0, censorvar, drop=FALSE]

      i.idx.ctl <- rep(1:nrow(Zi.ctl), each=nrow(Zj.ctl))
      j.idx.ctl <- rep(1:nrow(Zj.ctl), nrow(Zi.ctl))

      Z.ctl<-(Zj.ctl[j.idx.ctl,] - Zi.ctl[i.idx.ctl, ])
      Ts.ctl <- (Ti.ctl[i.idx.ctl, ] - Tj.ctl[j.idx.ctl,])
      dels.ctl <- Cs.ctl[i.idx.ctl, ]

      # set CV folds
      fold.ij.ctl <- sample.int(nfolds, nrow(Zi.ctl), replace = TRUE)
      fold.Z.ctl <- ifelse(fold.ij.ctl[i.idx.ctl] == fold.ij.ctl[j.idx.ctl],
                           fold.ij.ctl[i.idx.ctl], NA)

      weights.ctl <- 1*(Ts.ctl<0)*dels.ctl
      Z.ctl <- Z.ctl[weights.ctl==1,]
      fold.Z.ctl <- fold.Z.ctl[weights.ctl==1]
      Z <- rbind(Z.trt, Z.ctl)
      wt.Z <- rep(1, nrow(Z))
      fold.Z <- c(fold.Z.trt, fold.Z.ctl)
    }else{
      stop("optAUCX only handles binary (type=b) and survial (type=s) outcomes. \n")
    }
    yb.sim<-rbinom(nrow(Z),1,0.5)
    data.new<-cbind(yb.sim,Z)
    Y.new<-data.new[yb.sim==1,]
    X.new<-data.new[yb.sim==0,]
    yvar.new<-"yb.sim"
    X.new[,-1]<--X.new[,-1]
    data.fin<-rbind(Y.new,X.new)
    ##
    data.fin<-cbind(data.fin,c(rep(1,nrow(Z.trt)),rep(0,nrow(Z.ctl))))
    names(data.fin)[ncol(data.fin)]<-trtvar
    ##
    data.fin <- as.data.frame(data.fin)
    data.fin <- data.frame(data.fin, cvfolds=fold.Z)

    y=data.fin[,yvar.new]
    x=as.matrix(data.fin[,xvars])

    # fit cv.glmnet to get all choices of nonzero estimates for different values of tuning parameter in LASSO
    cv.idx <- !is.na(data.fin$cvfolds)
    fit.cv=cv.glmnet(x[cv.idx,], y[cv.idx], family="binomial", weights = wt.Z[cv.idx],
                     nfolds=5, foldid=data.fin$cvfolds[cv.idx], pmax=pmax,
                     intercept = FALSE, type.measure="auc", parallel = TRUE)

  }

  fit.cv.nzero.idx<-which(!duplicated(fit.cv$nzero))[-1]
  fit.cv.nzero.mat<-as.matrix(fit.cv$glmnet.fit$beta[,fit.cv.nzero.idx])

  betahat.mat <- fit.cv.nzero.mat
  rownames(betahat.mat)<-xvars
  colnames(betahat.mat)<-paste("order",1:ncol(fit.cv.nzero.mat),sep="")

  # fit.cv.nzero.list<-try(rank.list.func(fit.cv.nzero.mat,xvars),silent=T)
  # p<-length(xvars)
  # betahat.mat<-matrix(0,p,p)
  # rownames(betahat.mat)<-xvars
  # colnames(betahat.mat)<-paste("order",1:p,sep="")
  #
  # for (k in 1:p){
  #
  #   xvars.nzero.idx<-which(xvars %in% fit.cv.nzero.list[[k]])
  #   xvars.nzero<-xvars[xvars.nzero.idx]
  #   if (length(xvars.nzero)==1)
  #   {
  #     betahat.mat[xvars.nzero.idx,k]<-1
  #   } else {
  #     res <- optAUCX(outcome=yvar, predictor=xvars.nzero, censorvar=censorvar, treatment=trtvar,type=type, data=data, trt.ref=trtref, scale=F)
  #     betahat.mat[xvars.nzero.idx,k]<-res$beta
  #   }
  # }
  return(list(betahat.mat=betahat.mat, fit.cv=fit.cv))
}


#' Internal function used in betahat.fold.func
#'
#' @param betamat betahat matrix from cv.glmnet.
#' @param xvars a vector of predictor variables.
#'
#' @return order of each set of predictors as returned by cv.glmnet.
#'
rank.list.func<-function(betamat,xvars)
{
  if(ncol(betamat)==1){
    xvars.nzero<-rownames(betamat)[betamat!=0]
    xvars.nzero.list<-list(xvars.nzero)
    return(xvars.nzero.list[rep(1,length(xvars))])
  } else {
    xvars.nzero.list<-apply(betamat,2,function(x) return(names(x[which(x!=0)])))

    l<-length(ls(xvars.nzero.list))
    xvars.nzero.idx<-NULL
    for ( i in 1:l)
    {
      if (i==1)
      {
        xvars.new<-xvars.nzero.list[[i]]
      } else {
        xvars.new<-setdiff(xvars.nzero.list[[i]],xvars.nzero.list[[i-1]])
      }
      xvars.nzero.idx<-c(xvars.nzero.idx,rep(i,length(xvars.new)))
    }
    if (length(xvars.nzero.idx)<length(xvars)) xvars.nzero.idx<-c(xvars.nzero.idx,rep(xvars.nzero.idx[length(xvars.nzero.idx)],length(xvars)-length(xvars.nzero.idx)))

    return(xvars.nzero.list[xvars.nzero.idx])
  }
}
