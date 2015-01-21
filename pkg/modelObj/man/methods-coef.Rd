\name{coef}
\alias{coef}
\alias{coef,modelObjFit-method}

\title{Extract Model Coefficients}

\description{Extracts model coefficients from an object of class
             modelObjFit created by a call to modelObj::fit().}

\usage{ \S4method{coef}{modelObjFit}(object, ...) }

\arguments{

  \item{object }{an object of class modelObjFit.}

  \item{\dots  }{passed through to coef() of the modelObj regression method.}

}

\value{Coefficients extracted from the modelObjFit object 'object'.
       For standard model fitting
       classes the value returned is a named vector. If no
       coef() method is defined for the regression method specified
       in the governing modelObj, NULL is returned.}

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

    coef(fit.obj)

}
