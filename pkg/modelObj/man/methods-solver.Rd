\name{solver}
\alias{solver,modelObj-method}
\alias{solver}

\title{Retrieve Regression Method.}

\description{Retrieve the regression method from an object of class modelObj.}

\usage{
       \S4method{solver}{modelObj}(object, ...)
}

\arguments{

  \item{object }{an object of class modelObj.}

  \item{\dots  }{ignored.}

}

\value{
       Returns an object of class "character" giving the name of the regression
       method to be used to obtain parameter estimates.
}

\author{ Shannon T. Holloway <sthollow@ncsu.edu> }

\examples{
    #----------------------------------------------------#
    # Create modeling object using a formula
    #----------------------------------------------------#
    mo <- buildModelObj(model=Y ~ X1 + X2 + X3 + X4,
                        solver.method='lm')

    solver(mo)

}
