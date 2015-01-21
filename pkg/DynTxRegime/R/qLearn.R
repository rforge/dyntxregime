#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# qLearn - public function to perform a step of the Q-Learning algorithm       #
#          If an object of class QLearn is passed, it is assumed to be the     #
#          preceding step of the Q-Learning algorithm and models are fit       #
#          using the Ytilde variable of the QLearn object. If a vector         #
#          is passed, it is assumed that this is the first step in the         #
#          Q-Learning algorithm and models are fit using the response.         #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# moMain  : an object of class modelObj or a list of objects of class          #
#           modelObjSubset, which define the models and R methods to           #
#           be used to obtain parameter estimates and predictions              #
#           for the main effects component of the outcome regression.          #
#           See ?modelObj and/or ?modelObjSubset for details.                  #
#           NULL is an acceptable value if moCont is defined.                  #
#                                                                              #
# moCont  : an object of class modelObj or a list of objects of class          #
#           modelObjSubset, which define the models and R methods to           #
#           be used to obtain parameter estimates and predictions              #
#           for the contrasts component of the outcome regression.             #
#           See ?modelObj and/or ?modelObjSubset for details.                  #
#           NULL is an acceptable value if moMain is defined.                  #
#                                                                              #
# data    : data frame of covariates and treatment histories                   #
#                                                                              #
# response: response vector or object of class qLearn from a previous          #
#           Q-Learning step.                                                   #
#                                                                              #
# txName  : character string giving column header of treatment variable        #
#           in data                                                            #
#                                                                              #
# fSet    : A function.                                                        #
#           This argument allows the user to specify the subset of tx          #
#           options available to a patient.                                    #
#           The functions should accept as input either                        #
#           1) explicit covariate names as given in column names of data       #
#           2) a vector of covariates (i.e. a row of a data.frame)             #
#           and must return a vector of tx options available to the            #
#           patient                                                            #
#           Note this function is used for an INDIVIDUAL patient, matrix       #
#           results are not appropriate                                        #
#                                                                              #
# iter    : an integer                                                         #
#                                                                              #
#           >=1 if moMain and moCont are to be fitted iteratively              #
#           The value is the maximum number of iterations.                     #
#           Note the iterative algorithms is as follows:                       #
#           Y = Ymain + Ycont                                                  #
#            (1) hat(Ycont) = 0                                                #
#            (2) Ymain = Y - hat(Ycont)                                        #
#            (3) fit Ymain ~ moMain                                            #
#            (4) set Ycont = Y - hat(Ymain)                                    #
#            (5) fit Ycont ~ moCont                                            #
#            (6) Repeat steps (2) - (5) until convergence or                   #
#            a maximum of iter iterations.                                     #
#                                                                              #
#           <=0 moMain and moCont will be combined and fit as a single object. #
#                                                                              #
#           Either categorical or integer data can be provided for the tx.     #
#           If categorical, the fitted contrast and main effects are defined   #
#           relative to the base category {defined as levels()[1]}. The values #
#           may not be those returned by predict(object) if iterate fits are   #
#           used. If integer, the fitted contrast and main effects are defined #
#           relative to no tx (tx = 0).                                        #
#                                                                              #
#           Note that if iter <= 0, all non-model components of the            #
#           moMain and moCont must be identical                                #
#                                                                              #
# base   : An integer indicating the base tx or NULL (ordinal tx)              #
#                                                                              #
# ...    : ignored                                                             #
#                                                                              #
# Function returns an object of class QLearn                                   #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

qLearn <- function(...,
                   moMain,
                   moCont, 
                   data, 
                   response, 
                   txName, 
                   fSet = NULL, 
                   iter = 0L){


  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  #++++++                         Verify Input                         ++++++#
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  if( !is(data, 'data.frame') ) {

    UserError("input", 
              paste("'data' must be an object of class data.frame.", 
                    "Received an object of class ", 
                    paste(is(data), collapse=","), ".",
                    sep=""))

  }

  if( !is(txName, 'character') ) {
    msg <- "'txName', must be a character object."
    UserError("input", msg)
  }

  txVec <- try(data[,txName], silent = TRUE)

  if( is(txVec,"try-error") ) {
    UserError("input",
              paste(txName, " not found in data.", sep="") )
  }

  if( !is(txVec,"factor") ) {
    if( !isTRUE(all.equal(txVec, round(txVec,0L))) ) {
      UserError("input",
                "Treatment variable must be a factor or an integer.")
    }
    data[,txName] <- as.integer(round(data[,txName],0L))
  }

  if( is.list(moMain) ) {

    if( length(moMain) < 0.5 ) {
      UserError("input", 
                "'moMain' must have length > 0")
    }

    tst <- sapply(X = moMain, FUN = class) != 'ModelObjSubset'
    if( any(tst) ) {
      UserError("input", 
                paste("If class(moMain) == list, ",
                      "all elements must be of class modelObjSubset.\n", 
                      "Received objects of class ", 
                      paste(is(moMain), collapse=", "), ".",
                      sep=""))
    }

    if( is(fSet, "NULL") ) {
      UserError("input", 
                paste("When using objects of class modelObjSubset, ",
                      "fSet must be provided.", sep=""))
    }

    moMain <- new("ModelObjSubsetList",
                  loo = moMain)

  } else if( !is(moMain, "modelObj") && !is(moMain, "NULL") ) {

    UserError("input", 
              paste("If modeling the superset of treatment options, ",
                    "moMain must be of class modelObj.\n", 
                    "Received an object of class ", 
                    paste(is(moMain),collapse=","), ".",
                    sep=""))

  }

  if( is(moCont, "list") ){

    if( length(moCont) < 0.5 ) {
      UserError("input", 
                "'moCont' must have length > 0")
    }

    tst <- sapply(X = moCont, FUN = class) != 'ModelObjSubset'
    if( any(tst) ) {
      UserError("input", 
                paste("If class(moCont) == list, ",
                      "all elements must be of class modelObjSubset.\n", 
                      "Received a objects of class ", 
                      paste(is(moCont), collapse=","), ".",
                      sep=""))
    }

    if( is(fSet, "NULL") ) {
      UserError("input", 
                paste("When using objects of class modelObjSubset, ",
                      "fSet must be provided.", sep=""))
    }

    moCont <- new("ModelObjSubsetList",
                  loo = moCont)

  } else if( !is(moCont, "modelObj") && !is(moCont, "NULL") ) {

    UserError("input", 
              paste("If modeling the superset of treatment options, ",
                    "moCont must be of class modelObj.\n", 
                    "Received an object of class ", 
                    paste(is(moCont),collapse=","), ".",
                    sep=""))

  }

  if( is(moMain, "NULL") && is(moCont, "NULL") ){
    UserError("input", 
              "Must provide at least one of {moMain, moCont}.")
  }

  if( is(response, "QLearn")  ){
    step <- IStep(response) + 1L
    response <- YTilde(response)
  } else if( is(response, "vector") ){
    step <- 1L
  } else {
    UserError("input", 
              paste("'response' must be a vector of responses or ",
                    "an object returned by a prior call to qLearn().",
                    sep=""))
  }

  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
  #++++++                         Calculation                          ++++++#
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

  cat("\n\nStep ", step, " of Q-learning algorithm.\n\n")

  if( is(fSet, "NULL") || is(fSet, "function") ){

    txInfo <- txProcess(txVar = txName, 
                        data = data, 
                        fSet = fSet)

  } else {

    UserError("input", 
              paste("If provided, fSet must be an object of class function.\n", 
                    "Received an object of class ", 
                    paste(is(fSet),collapse=","), ".",
                    sep=""))

  }

  est <- qLearnEst(moMain = moMain, 
                   moCont = moCont,
                   data = data, 
                   response = response, 
                   iter = iter, 
                   txInfo = txInfo)

  #--------------------------------------------------------------------------#
  # Calculate Q-Functions at each tx                                         #
  #--------------------------------------------------------------------------#
  if( is(data[,txName],"factor") ) {
    nms <- levels(data[,txName])
  } else {
    nms <- SuperSet(txInfo)
  }

  qFunctions <- matrix(data = 0.0,
                       nrow = nrow(data),
                       ncol = length(nms),
                       dimnames = list(NULL, nms))

  for( i in 1L:length(nms) ) {

    if( is(data[,txName],"factor") ) {
      data[, txName] <- factor(rep(levels(data[,txName])[i],nrow(data)), 
                               levels = levels(data[,txName]))
    } else {
      data[, txName] <- as.integer(nms[i])
    }

    qFunctions[,i] <- PredictMain(object = est, newdata = data) +
                      PredictCont(object = est, newdata = data)

  }

  q2opt <- max.col(qFunctions, ties.method="first")
  if( is(data[,txName], "factor") ) {
    optTx <- factor(colnames(qFunctions)[q2opt],
                    levels = colnames(qFunctions))
  } else {
    optTx <- as.integer(colnames(qFunctions)[q2opt])
  }

  result <- new("QLearn", 
                step = step,
                qFunctions = qFunctions,
                optimalTx = optTx,
                call = match.call(),
                est,
                txInfo)

  show(result)

  return(result)
}

