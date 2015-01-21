\name{fit}
\alias{fit}
\alias{fit,modelObj,data.frame,vector-method}
\alias{predictor,modelObjFit-method}
\alias{predictorArgs,modelObjFit-method}

\title{Execute Regression Step}

\description{Executes the model fitting procedure for an object of class modelObj.}

\usage{
\S4method{fit}{modelObj,data.frame,vector}(object, data, response, ...)
}

\arguments{

  \item{object  }{an object of class modelObj as returned by the 
                  `buildModelObj' function.}

  \item{data    }{an object of class data.frame 
                  containing the variables in the model.}

  \item{response}{an object of class vector containing the response variable.}

  \item{\dots   }{ignored}

}

\details{If defined by the modeling function, the 
         following methods can be applied to a modelObjFit object to obtain 
         standard fit results:
         \code{coef},
         \code{plot},
         \code{predict},
         \code{residuals}, and
         \code{summary}.}

\value{
       Returns an object of class "modelObjFit", which contains
       the object returned by the modeling function and
       the method to be used to obtain predictions.

       Slots of the modelObjFit object can be retrieved using
       \code{fitObject(object)},
       \code{predictor(object)}, and
       \code{predictorArgs(object)}.
}

\author{
  Shannon T. Holloway <sthollow@ncsu.edu>
}

\examples{
    #----------------------------------------------------#
    # Generate data
    #----------------------------------------------------#
    X <- matrix(rnorm(1000,0,1),
                ncol=4,
                dimnames=list(NULL,c("X1","X2","X3","X4")))

    Y <- X \%*\% c(0.1, 0.2, 0.3, 0.4) + rnorm(250)

    X <- data.frame(X)

    #----------------------------------------------------#
    # Create modeling object using a formula
    #----------------------------------------------------#
    mo <- buildModelObj(model=Y ~ X1 + X2 + X3 + X4,
                   solver.method='lm')

    #----------------------------------------------------#
    # Fit model
    #----------------------------------------------------#
    fit.obj <- fit(object=mo, data=X, response=Y)

    coef(fit.obj)
    head(residuals(fit.obj))
    plot(fit.obj)
    head(predict(fit.obj,X))
    summary(fit.obj)

}
