% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimalAUCX_functions_surv.R
\name{opt.idx.func}
\alias{opt.idx.func}
\title{Internal function used in xselect.glmnet.}
\usage{
opt.idx.func(m)
}
\arguments{
\item{m}{the matrix of median and se from cross-validated AUC.}
}
\value{
row number corresponding to optimum AUC/ABC according to 1se rule.
}
\description{
Internal function used in xselect.glmnet.
}
