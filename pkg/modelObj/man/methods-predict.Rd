\name{predict}
\alias{predict}
\alias{predict,modelObjFit-method}

\title{Model Prediction.}

\description{Executes the prediction method for an object of class
             modelObjFit.}

\usage{ \S4method{predict}{modelObjFit}(object, newdata = NULL, ...) }

\arguments{

  \item{object }{an object of class modelObjFit containing the
                 regression object which which predictions are 
                 to be obtained.}

  \item{newdata }{An optional data frame of variables with which to make
                  predictions. If not provided, fitted values from the
                  regression analysis are returned.}

  \item{\dots  }{ignored}

}

\value{
       A  matrix of the predicted response for the fitted model.
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

    pred <- predict(object=fit.obj, newdata=X)

    head(pred)

}
