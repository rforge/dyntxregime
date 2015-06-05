\name{buildmodelObj}
\alias{buildModelObj}

\title{
Create an Object of Class modelObj
}

\description{
A utility function to transfer user defined
models and estimation methods to an object of class modelObj.
}

\usage{ 
buildModelObj(model, solver.method=NULL, solver.args=NULL, 
predict.method=NULL, predict.args=NULL)
}
\arguments{

\item{model}{A formula object defining the regression model.}
 
\item{solver.method}{An object of class character specifying 
                     the name of the R function to be used to obtain 
                     parameter estimates. 
                     For example,  `lm', `glm', or `rpart'. 
                     The specified modeling function MUST
                     have a corresponding predict method.}
 
\item{solver.args  }{An object of class list containing 
                     additional arguments to be sent to solver.method. 
                     Arguments must be provided as a list, where the name of  
                     each element matches a formal argument of 
                     solver.method. For example, if a logistic regression 
                     using glm is desired, 
                     \deqn{solver.method = ``glm"}
                     \deqn{solver.args = list(``family"=binomial)}

                     It is assumed that the R method  
                     specified in solver.method has
                     formal arguments ``formula" and ``data".
                     If solver.method does not use 
                     ``formula" for the formula object or ``data" as the  
                     data.frame, solver.args must explicitly indicate
                     the variable names used for these inputs.
                     For example, list(``x"=``formula") if the formula 
                     object is passed to solver.method through input argument
                     ``x" or list(``df"=``data") if the data.frame object
                     is passed to solver.method through input argument ``df." }
 
\item{predict.method}{A character.
                     The name of the R function to be used to obtain 
                     predictions. 
                     For example,  `predict.lm', `predict', or `predict.glm'. 
                     If no function is explicitly 
                     given, the generic \code{predict} is assumed.
                     For many methods, the generic method is appropriate. }
 
\item{predict.args }{A list. 
                     Additional arguments to be sent to 
                     predict.method. This must be provided as a list, 
                     where the name of each element matches a formal 
                     argument of predict.method. For example, 
                     if a logistic regression using glm was used to fit the 
                     model formula object and predictions on the scale of 
                     the response are desired, 
                     \deqn{predict.method = ``predict.glm"}
                     \deqn{predict.args = list(``type"=``response").}

                       It is assumed that the 
                       predict.method has formal arguments
                       ``object" and ``newdata". If predict.method does not 
                       use these 
                       formal arguments, predict.args must explicitly indicate
                       the variable names used for these inputs. 
                       For example, list(``x"=``object") if the object 
                       returned by solver.method is passed to predict.method 
                       through input argument ``x"  or list(``ndf"=``newdata") 
                       if the data.frame object is passed to predict.method 
                       through input argument ``ndf." .}
 

 
}

\details{
Unless changed by the user in solver.args and/or
predict.args, default settings are assumed for the
specified methods.  
}

\value{
 An object of type modelObj, which
  contains a complete description of the model, 
  the method to be used for parameter estimates, and the prediction method.
}


\examples{

    #----------------------------------------------------#
    # Create modeling object using a formula
    #----------------------------------------------------#
    mo <- buildModelObj(model=Y ~ X1 + X2 + X3 + X4,
                        solver.method='lm', 
                        predict.method='predict.lm',
                        predict.args=list(type='response'))

}
