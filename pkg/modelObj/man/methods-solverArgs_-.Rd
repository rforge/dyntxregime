\name{solverArgs<-}
\alias{solverArgs<-}
\alias{solverArgs<-,modelObj-method}

\title{Reset Arguments for Regression Method}

\description{Reset the arguments to be sent to the regression method.
             This capability is intended for package 
             developers making use of the ``model object" framework.}

\usage{ \S4method{solverArgs}{modelObj}(object) <- value}

\arguments{

  \item{object}{object of class modelObj, for which the arguments sent to the
                regression method are to be reset.}

  \item{value }{A named list containing the new arguments to be sent to the 
                regression method.}

}
\details{
The first two elements of the list contain the name for the formula object
 and the name for the data.frame. buildModelObj creates modelObj using an internal
convention, which is critical to the correct implementation of this package. 
These two elements should not be changed and will be carried over from the 
original argument list.
}
\value{modelObj is updated with new argument list.}

\author{ Shannon T. Holloway <sthollow@ncsu.edu> }

\examples{
    #----------------------------------------------------#
    # Create modeling object using a formula
    #----------------------------------------------------#
    mo <- buildModelObj(model=Y ~ X1 + X2 + X3 + X4,
                        solver.method='lm')

    solverArgs(mo)
    argList <- list("x"=TRUE)
    solverArgs(mo) <- argList
    solverArgs(mo)
}
