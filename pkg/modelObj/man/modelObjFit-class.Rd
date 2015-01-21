\name{modelObjFit-class}
\Rdversion{0.4}
\docType{class}
\alias{modelObjFit-class}
\alias{show,modelObjFit-method}
\alias{summary,modelObjFit-method}

\title{Class \code{"modelObjFit"}}
\description{A class containing the fit results for a model object
             defined as an object of class modelObj.
}

\section{Objects from the Class}{
  Objects should not be created directly. Objects are created by a call to
  fit().
}
\section{Slots}{
  \describe{
    \item{\code{fitObj}:   }{Value object of modeling function
                             (i.e., class(fitObj)="lm" if solver.method=lm). }
    \item{\code{func}:     }{Object of class \code{methodObj} 
                             the method to be used to obtain predicted values. }
  }
}
\section{Methods}{
  \describe{

    \item{coef     }{\code{signature(object = "modelObjFit")}: if, 
                     defined for the modeling function, extracts
                     model coefficients. Else, returns NULL.}

    \item{fitObject  }{\code{signature(object = "modelObjFit")}: 
                       Retrieve fit object. }

    \item{plot     }{\code{signature(x = "modelObjFit", y=NULL)}: if 
                     defined for the modeling function, invokes
                     plot(modelObj@fitObj). Else, returns a message that a 
                     plot(modelObj@fitObj) method is not defined.}

    \item{predict  }{\code{signature(object = "modelObjFit")}: 
                     returns predicted values obtained using fitObj and
                     the newdata data.frame.}

    \item{predictor  }{\code{signature(object = "modelObjFit")}: 
                 Retrieve prediction method name. }

    \item{predictorArgs  }{\code{signature(object = "modelObjFit")}: 
                 Retrieve arguments to be sent to prediction method. }

    \item{residuals}{\code{signature(object = "modelObjFit")}: if 
                     defined for the modeling function, returns the
                     residuals(modelObjFit@fitObj).}

    \item{show     }{\code{signature(x = "modelObjFit")}: displays
                     standard results returned by modeling function.}

    \item{summary  }{\code{signature(x = "modelObjFit")}: if, 
                     defined for the modeling function, returns 
                     summary(modelObjFit@fitObj). Else, returns a message that a 
                     summary(modelObjFit@fitObj) method is not defined.}
	 }
}
\examples{
showClass("modelObjFit")
}
\keyword{classes}
