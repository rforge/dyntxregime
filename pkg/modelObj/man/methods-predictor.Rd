\name{predictor}
\alias{predictor}
\alias{predictor,modelObj-method}

\title{ Retrieve Prediction Method }

\description{ Retrieve the prediction method from an object of class modelObj.
              This capability is intended for package 
              developers making use of the ``model object" framework.}

\usage{
\S4method{predictor}{modelObj}(object, ...)
}

\arguments{

  \item{object }{an object of class modelObj.}

  \item{\dots  }{ignored.}

}

\value{
       Returns an object of class "character"; the name of the prediction
       method specified in the modelObj object.
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

    predictor(mo)

}
