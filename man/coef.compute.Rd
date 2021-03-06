% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Coef_comp.R
\name{coef.compute}
\alias{coef.compute}
\title{Model Coefficient Generator}
\usage{
\method{coef}{compute}(y, grpvec, pooling, ...)
}
\arguments{
\item{y}{vector of data}

\item{grpvec}{pooling vector of interest}

\item{pooling}{corresponding pooling to grpvec as generated by grp.asn()}
}
\description{
This function generates the model coefficients corresponding to the model fit 
for a given partition.
}
\examples{
coef.compute()
}
