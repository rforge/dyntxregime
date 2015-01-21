\name{residuals}
\alias{residuals}
\alias{residuals,modelObjFit-method}

\title{Extract Model Residuals.}

\description{Returns the residuals from an object of class
             modelObjFit.}

\usage{ \S4method{residuals}{modelObjFit}(object, ...) }

\arguments{

  \item{object }{an object of class modelObjFit containing the 
                 regression object from which residuals are to be
                 retrieved.}

  \item{\dots  }{ignored}

}

\value{
       If \code{residuals()} is defined for the regression method
       used to obtain parameter estimate, returns the 
       residuals as defined by the regression method. 
       If method does not exist, a warning messages is printed.
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
    mo <- buildModelObj(model = Y ~ X1 + X2 + X3 + X4,
                        solver.method = 'lm')

    #----------------------------------------------------#
    # Fit model
    #----------------------------------------------------#
    fit.obj <- fit(object=mo, data=X, response=Y)

    res <- residuals(fit.obj)
    head(res)

}
