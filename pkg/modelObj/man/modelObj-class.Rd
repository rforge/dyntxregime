\name{modelObj-class}
\Rdversion{0.4}
\docType{class}
\alias{modelObj-class}

\title{Class \code{"modelObj"}}
\description{
A class for model objects.
}
\section{Objects from the Class}{
Objects should not be created directly. The utility function
buildModelObj() should be used.
}
\section{Slots}{
  \describe{
    \item{\code{model}:    }{Object of class \code{formula} }
    \item{\code{solver}:   }{Object of class \code{methodObj} 
                             method to obtain parameter estimates.}
    \item{\code{predictor}:}{Object of class \code{methodObj} 
                             method to obtain predicted values.}
  }
}
\section{S4Methods}{
  \describe{
    \item{fit  }{\code{signature(object = "modelObj", data="data.frame", 
                       response="vector", ...)}: 
                 Executes regression step. }
    \item{model }{\code{signature(object = "modelObj")}: 
                 Retrieve model. }
    \item{solver  }{\code{signature(object = "modelObj")}: 
                 Retrieve regression method name. }
    \item{solverArgs  }{\code{signature(object = "modelObj")}: 
                 Retrieve arguments to be sent to regression method. }
    \item{solverArgs(object)<-  }{\code{signature(object = "modelObj")}: 
                 Set arguments to be sent to regression method. }
    \item{predictor  }{\code{signature(object = "modelObj")}: 
                 Retrieve prediction method name. }
    \item{predictorArgs  }{\code{signature(object = "modelObj")}: 
                 Retrieve arguments to be sent to prediction method. }
    \item{predictorArgs(object)<-  }{\code{signature(object = "modelObj")}: 
                 Set arguments to be sent to prediction method. }
  }
}
\examples{
showClass("modelObj")
}
\keyword{classes}
