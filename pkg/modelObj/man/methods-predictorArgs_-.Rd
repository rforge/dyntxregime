\name{predictorArgs<-}
\alias{predictorArgs<-}
\alias{predictorArgs<--methods}
\alias{predictorArgs<-,modelObj-method}

\title{Reset Arguments for Prediction Method}

\description{Reset the arguments to be sent to the prediction method.
             This capability is intended for package developers making use of
             the ``model object" framework. }

\usage{\S4method{predictorArgs}{modelObj}(object) <- value}
\arguments{

  \item{object}{object of class modelObj, for which the arguments sent to the
                prediction method are to be reset.}

  \item{value }{A named list containing the new arguments to be sent to the 
                prediction method.}

}
\details{
The first two elements of the list contain the name for the object returned
by the regression method and the name for the data.frame for which
predictions are desired. buildModelObj creates modelObj using an internal
convention, which is critical to the correct implementation of this package. 
These two elements should not be changed and will be carried over from the 
original argument list.
}

\value{modelObj with an updated argument list to be passed to the
       prediction method.}

\author{ Shannon T. Holloway <sthollow@ncsu.edu> }

\examples{
    #----------------------------------------------------#
    # Create modeling object using a formula
    #----------------------------------------------------#
    mo <- buildModelObj(model=Y ~ X1 + X2 + X3 + X4,
                        solver.method='lm', 
                        predict.method='predict.lm',
                        predict.args=list(type='response'))

    predictorArgs(mo)
    argList <- list("type"='terms')
    predictorArgs(mo) <- argList
    predictorArgs(mo)
}
