\name{optTx}
\alias{optTx}
\alias{optTx,DynTxRegime,data.frame-method}
\alias{optTx,DynTxRegime,missing-method}

\title{ Methods to Retrieve Estimated or to Predict Optimal Treatment}

\description{
Retrieves the estimated or predicted optimal treatment.
}

\usage{
\S4method{optTx}{DynTxRegime,missing}(x, newdata, ...)  
\S4method{optTx}{DynTxRegime,data.frame}(x, newdata, ...)  
\S4method{optTx}{DynTxRegime,data.frame}(x, newdata, dp, ...)  
}

\arguments{
  \item{x }{an object of class DynTxRegime.}

  \item{newdata }{an object of class data.frame containing
                  covariates and treatment histories.}

  \item{dp} {an object of class integer. For multi-decision-point
             value-search, dp is the decision point for which
             predictions are requested.}

  \item{\dots }{ignored.}

}

\value{
       For outcome regression methods Q-Learning and IQ-Learning a list
       \item{qFunctions}{ a matrix of estimated Q-Functions.}
       \item{optimalTx}{a vector of estimated optimal treatments.}
       For value-search methods optimalClass and optimalSeq a 
       a matrix of estimated optimal treatments.
}

