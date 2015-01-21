\name{solverArgs}
\alias{solverArgs}
\alias{solverArgs,modelObj-method}

\title{Retrieve Arguments for Regression Method.}

\description{Retrieve the arguments to be sent to the 
             regression method.}

\usage{ \S4method{solverArgs}{modelObj}(object, ...)}

\arguments{

  \item{object}{object of class modelObj, from which the arguments sent to the
                regression method are to be retrieved.}

  \item{...   }{ignored.}

}

\author{
  Shannon T. Holloway <sthollow@ncsu.edu>
}

\examples{
    #----------------------------------------------------#
    # Create modeling object using a formula
    #----------------------------------------------------#
    mo <- buildModelObj(model = Y ~ X1 + X2 + X3 + X4,
                   solver.method = 'lm', 
                   predict.method = 'predict.lm',
                   predict.args = list(type='response'))

    solverArgs(mo)
    solverArgs(mo) <- list(model=TRUE)
    solverArgs(mo)

}
