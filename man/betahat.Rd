% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimalAUCX_functions_surv.R
\name{betahat}
\alias{betahat}
\title{optAUC estimates for prognostic case.}
\usage{
betahat(
  yvar,
  xvars,
  censorvar = NULL,
  data,
  type = "b",
  xselect = FALSE,
  cv.iter.xselect = 20,
  k.fold.xselect = 5,
  method.xselect = "CV"
)
}
\arguments{
\item{yvar}{outcome variable (column) name for response variable assuming data is provided.}

\item{xvars}{vector of variable names for predictors (covariates) assuming data is provided.}

\item{censorvar}{name of censoring variable (1: event; 0: censor), default = NULL.}

\item{data}{the data frame for training dataset.}

\item{type}{the type of response variable: "b" binary; "s" survival.}

\item{xselect}{a logical flag. TRUE: use Cross Validated AUC or aAUC for variable selection; FALSE: use all predictors to calculate the composite score. Default = FALSE}

\item{cv.iter.xselect}{the number of iterations of cross-validation used for MCglmnet.}

\item{k.fold.xselect}{the number of folds required for cross-validated AUC.}

\item{method.xselect}{the method used for variable selection: CV (default) cross validated AUC; "aAUC" approximate AUC.}
}
\value{
This function returns a list of the following:
\item{beta}{the estimated linear coefficients.}
\item{bias}{the HIC for the selected model.}
}
\description{
Fit optAUC to the data in prognostic case. Estimates along with the model bias correction are returned.
}
