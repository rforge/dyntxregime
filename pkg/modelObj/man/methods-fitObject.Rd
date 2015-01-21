\name{fitObject}
\alias{fitObject}
\alias{fitObject,modelObjFit-method}

\title{ Retrieve Regression Object }

\description{ Retrieves the value object returned by the regression method
              used to obtain parameter estimates. }

\usage{ \S4method{fitObject}{modelObjFit}(object, ...) }

\arguments{

  \item{object }{an object of class modelObjFit.}

  \item{\dots  }{ignored.}

}

\details{This function is useful for accessing methods that are defined
         by the regression method but are not directly accessible from
         the modelObjFit object. For example, for many regression methods,
         users can retrieve the fitted values by calling fitted(object). 
         This method is not directly accessible from a modelObjFit. However,
         fitted() can be applied to the object returned by fitObject().}

\value{The Value returned by the regression method specified in
       the governing modelObj. The exact structure of the value will
       depend on the regression method. For example, if nls() is the
       regression method, a list is returned.  }

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

    obj <- fitObject(fit.obj)
    fobj <- fitted(obj)
    head(fobj)

}
