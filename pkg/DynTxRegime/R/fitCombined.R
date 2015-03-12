#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# fitCombined - Combines main effects and contrasts models to obtain           #
#               parameter estimates in a single call to fit method.            #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# moMain   : an object of class modelObj that defines the model and R methods  #
#            to be used to obtain parameter estimates and predictions for the  #
#            main effects component of an outcome regression step.             #
#            See ?modelObj for details.                                        #
#            NULL is an acceptable value for moMain if moCont is defined.      #
#            However, at least one of {moMain,moCont} must be defined.         #
#                                                                              #
# moCont   : an object of class modelObj that defines the model and R methods  #
#            to be used to obtain parameter estimates and predictions for the  #
#            contrasts component of an outcome regression step.                #
#            See ?modelObj for details.                                        #
#            NULL is an acceptable value for moMain if moCont is defined.      #
#            However, at least one of {moMain,moCont} must be defined.         #
#                                                                              #
# response : response vector                                                   #
#                                                                              #
# txName   : column header of data containing tx variable                      #
#                                                                              #
# data     : data.frame of covariates                                          #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns an object of class SimpleFit                                       =#
#=                                                                            =#
#==============================================================================#
fitCombined <- function(moMain,
                        moCont,
                        response,
                        txName,
                        data){

  #--------------------------------------------------------------------------#
  # Recast response as matrix with column header YinternalY                  #
  #--------------------------------------------------------------------------#
  YinternalY <- matrix(data = response, 
                       ncol = 1L,  
                       dimnames = list(NULL,"YinternalY"))

  if( is(data[,txName], "factor") ) data[,txName] <- as.factor(data[,txName])

  #--------------------------------------------------------------------------#
  # Combine main effects and contrast models                                 #
  #--------------------------------------------------------------------------#
  if( !is(moCont, "NULL") && !is(moMain, "NULL") ){

    type <- "Combined"

    tempCont <- terms(model(moCont))
    contPart <- paste(attr(tempCont,"term.labels"), collapse="+")

    tempMain <- terms(model(moMain))
    mainPart <- paste(attr(tempMain,"term.labels"), collapse="+")

    if( attr(tempMain,"intercept") < 0.5 ) {
      mainPart <- paste("0 + ", mainPart, sep="")
    }

    if( attr(tempCont,"intercept") > 0.5 ) {
      newForm <- paste("~", mainPart, "+ ", txName, " + ",
                       txName, ":(", contPart, ")", sep="")
    } else if( attr(tempCont,"intercept") < 0.5 ) {
      newForm <- paste("~", mainPart, "+ ",
                       txName, ":(", contPart, ")", sep="")
    } 

    newForm <- as.formula(newForm)

    s <- solver(moMain)
    sa <- solverArgs(moMain)
    sa[[1L]] <- "formula"
    sa[[2L]] <- "data"
    p <- predictor(moMain)
    pa <- predictorArgs(moMain)
    pa[[1L]] <- "object"
    pa[[2L]] <- "newdata"

    obj <- buildModelObj(model = newForm, 
                         solver.method = s,
                         solver.args = sa,
                         predict.method = p,
                         predict.args = pa)

  } else if( !is(moCont, "NULL") && is(moMain, "NULL") ) {

    type <- "moCont"

    tempCont <- terms(model(moCont))
    contPart <- paste(attr(tempCont,"term.labels"), collapse="+")

    if( attr(tempCont,"intercept") > 0.5 ) {
      newForm <- paste("~ 0 + ", txName, " + ", 
                       txName, ":(", contPart, ")", sep="")
    } else if( attr(tempCont,"intercept") < 0.5 ){
      newForm <- paste("~ 0 + ", txName, ":(", contPart, ")", sep="")
    } 
    newForm <- as.formula(newForm)

    s <- solver(moCont)
    sa <- solverArgs(moCont)
    sa[[1L]] <- "formula"
    sa[[2L]] <- "data"
    p <- predictor(moCont)
    pa <- predictorArgs(moCont)
    pa[[1L]] <- "object"
    pa[[2L]] <- "newdata"

    obj <- buildModelObj(model = newForm, 
                         solver.method = s,
                         solver.args = sa,
                         predict.method = p,
                         predict.args = pa)

  } else if( is(moCont, "NULL") && !is(moMain, "NULL")){

    type <- "moMain"

    obj <- moMain

  } else {

    DeveloperError("moMain and moCont made it in as null.", "fitCombined")

  }

  #--------------------------------------------------------------------------#
  # Obtain fit                                                               #
  # (fit method of package modelObj)                                         #
  #--------------------------------------------------------------------------#
  fitCombined <- fit(object = obj, 
                     data = data, 
                     response = YinternalY)

  fittedY <- predict(object = fitCombined,
                     newdata = data)

  n <- nrow(data)
  txVec <- data[,txName]

  if( is(data[,txName], "factor") ) {

    #----------------------------------------------------------------------#
    # Set tx to base case to obtain main effects contribution to predicted #
    # response. (predict method of package modelObj)                       #
    #----------------------------------------------------------------------#
    base <- levels(data[,txName])[1L]
    data[,txName] <- factor(rep(base,n),
                            levels = levels(data[,txName]))
  } else {

    #----------------------------------------------------------------------#
    # Set tx to zero  to obtain main effects contribution to predicted     #
    # response. (predict method of package modelObj)                       #
    #----------------------------------------------------------------------#
    base <- 0L
    data[,txName] <- base

  }

  fittedMain <- predict(object = fitCombined, 
                        newdata = data)

  fittedCont <- (fittedY - fittedMain)
  if( !is(txVec,"factor") ) {
    fittedCont <- fittedCont/txVec
    fittedCont[is.infinite(fittedCont)] <- 0.0
  }

  residuals <- as.numeric(YinternalY - fittedY)

  result <- new("SimpleFit",
                "baseLevel" = base,
                "txName" = txName,
                "fitType" = type, 
                "modelObjectFit" = fitCombined,
                "residuals" = residuals,
                "yMainHat" = as.vector(fittedMain),
                "yContHat" = as.vector(fittedCont))

  return(result)
}
