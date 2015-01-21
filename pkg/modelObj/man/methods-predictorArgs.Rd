\name{predictorArgs}
\alias{predictorArgs}
\alias{predictorArgs,modelObj-method}

\title{Retrieve Arguments Specified for the Prediction Method}

\description{Retrieve from objects of class modelObj created by
             function buildModelObj the arguments that will be sent to the 
             prediction method. This capability is intended for package 
             developers making use of the ``model object" framework. }

\usage{
\S4method{predictorArgs}{modelObj}(object, ...)
}

\arguments{

  \item{object }{an object of class modelObj.}

  \item{\dots  }{ignored.}

}

\value{
       Returns an object of class "list" containing the arguments
       specified to be passed to the prediction method.
}

\author{
  Shannon T. Holloway <sthollow@ncsu.edu>
}

\examples{
    #----------------------------------------------------#
    # Create modeling object using a formula
    #----------------------------------------------------#
    mo <- buildModelObj(model=Y ~ X1 + X2 + X3 + X4,
                        solver.method='lm', 
                        predict.method='predict.lm',
                        predict.args=list(type='response'))

    predictorArgs(mo)

}
