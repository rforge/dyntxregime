\name{plot}
\alias{plot}
\alias{plot,modelObjFit-method}

\title{ Plotting }

\description{ Generates the plot from an object of class modelObjFit. }

\usage{ \S4method{plot}{modelObjFit}(x, y, ...) }

\arguments{

  \item{x    }{an object of class modelObjFit.}
  \item{y    }{ignored}
  \item{\dots}{ignored}

}

\value{
       If \code{plot} is defined for the regression method
       that was used to obtain the regression parameters, 
       generates the plot as defined by the regression method. 
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

    plot(fit.obj)

}
