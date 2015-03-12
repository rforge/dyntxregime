#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# optimalCore : obtains parameter estimates for propensity of treatment and    #
#               outcome regression models for both the DR and classification   #
#               methods.                                                       #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# moPropen: a single object of class modelObj, a list of objects of class      #
#           modelObj, or a list of objects of class modelObjSubset.            #
#           Define the models and R methods to be used to obtain               #
#           parameter estimates and predictions for the propensity for         #
#           treatment.                                                         #
#           See ?modelObj and/or ?modelObjSubset for details.                  #
#                                                                              #
#           If the prediction method specified in moPropen returns             #
#           predictions for only a subset of the categorical tx data,          #
#           it is assumed that the base level defined by levels(tx)[1] is      #
#           the missing category.                                              #
#                                                                              #
# moMain  : a single object of class modelObj, a list of objects of class      #
#           modelObj, or a list of objects of class modelObjSubset.            #
#           Define the models and R methods to be used to obtain               #
#           parameter estimates and predictions for the main effects           #
#           component of the outcome regression.                               #
#           See ?modelObj and/or ?modelObjSubset for details.                  #
#           NULL is an appropriate value.                                      #
#                                                                              #
# moCont  : a single object of class modelObj, a list of objects of class      #
#           modelObj, or a list of objects of class modelObjSubset.            #
#           Define the models and R methods to be used to obtain               #
#           parameter estimates and predictions for the contrasts              #
#           component of the outcome regression.                               #
#           See ?modelObj and/or ?modelObjSubset for details.                  #
#           NULL is an appropriate value.                                      #
#                                                                              #
# data    : a data frame of the covariates and tx histories                    #
#           tx variable will be recast as factor if not provided as such.      #
#                                                                              #
# response: response vector                                                    #
#                                                                              #
# txInfo  : an object of class txInfo or a list of objects of class            #
#           txInfo where each element of the list corresponds to a             #
#           single decision point.                                             #
#                                                                              #
# refit   : TRUE/FALSE flag indicating if conditional expectations             #
#           are to be refit for each new set of tx regime                      #
#           parameters (TRUE), or Q-learning is to be used (FALSE).            #
#           Can only be 'true' for sequential method.                          #
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
#           Note that if iter <= 0, all non-model components of the            #
#           moMain and moCont must be identical                                #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns a list                                                             =#
#=  $outcome     : an object of class QLearnEst or QLearnEstList              =#
#=                                                                            =#
#=  $propensity  : an object of class PropenFit or PropenFitList              =#
#=                                                                            =#
#==============================================================================#
optimalCore <- function(moPropen,
                        moMain = NULL,
                        moCont = NULL,
                        data = NULL,
                        response = NULL,
                        txInfo = NULL,
                        refit = FALSE,
                        iter = 0L){

  if( is(txInfo,'TxInfoList') ) {
    nDP <- length(txInfo)
  } else if( is(txInfo,'TxInfo') ) {
    nDP <- 1L
  } else {
    DeveloperError("txInfo not of appropriate class", "optimalCore")
  }

  objs <- list("moPropen" = moPropen, 
               "moMain" = moMain, 
               "moCont" = moCont)

  for( i in 1L:3L ) {

    if( is(objs[[i]], 'NULL') ) next

    if( is(objs[[i]], 'modelObj') ){

      if( nDP != 1L ) {
        UserError("input",
                  paste("Insufficient number of models specified for ",
                        names(objs)[i],".", sep=""))
      }

    } else if( is(objs[[i]],'ModelObjList') ){

      if( as.integer(round(length(objs[[i]]),0L)) != nDP ) {
        UserError("input",
                  paste("Incorrect number of models specified for ",
                        names(objs)[i],".", sep=""))
      }

    } else if( is(objs[[i]],'ModelObjSubsetList') ) {

      if( length(objs[[i]]) < nDP - 1.5e-8 ) {
          UserError("input",
                    paste("Insufficient number of models specified for ",
                          names(objs)[i],".", sep=""))
      }

    } else if( !is(objs[[i]][[1]], "NULL") ) {

      DeveloperError(paste(objs[[i]], "made it in as ",
                           paste(is(objs[[i]]),collapse=","), sep=""),
                     "optimalCore")
    }
  }


  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #                            Obtain tx probabilities                       #
  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

  propen <- list()

  for( i in 1L:nDP ) {

    mods <- pickPropenModels(moPropen, i)

    if( is(txInfo,'TxInfo') ){
      txI <- txInfo
    } else if( is(txInfo,'TxInfoList') ){
      txI <- txInfo[[i]]
    } else {
      DeveloperError("txInfo not of appropriate class", "optimalCore")
    }

    propen[[i]] <- fitPropen(moPropen = mods,
                             txInfo = txI,  
                             data = data)

  }

  if( as.integer(round(length(propen),0L)) == 1L ) {
    propen <- propen[[1]]
  } else {
    propen <- new("PropenFitList",
                  loo = propen)
  }

  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
  #                     Process conditional expectations                     #
  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

  outcome <- list()

  if( !is(moMain, "NULL") || !is(moCont, "NULL") ){

    #----------------------------------------------------------------------#
    # Fit the last decision point model                                    #
    #----------------------------------------------------------------------#

    mods <- pickOutcomeModels(moMain, moCont, nDP)

    if( is(txInfo,'TxInfo') ){
      txI <- txInfo
    } else if( is(txInfo,'TxInfoList') ){
      txI <- txInfo[[nDP]]
    } else {
      DeveloperError("txInfo not of appropriate class", "optimalCore")
    }

    outcome[[nDP]] <- qLearnEst(moMain = mods$moMain,
                                moCont = mods$moCont,
                                data = data,
                                response = response,
                                iter = iter,
                                txInfo = txI)

    j <- nDP - 1L

    while(j > 0L){

      mods <- pickOutcomeModels(moMain, moCont, j)

      outcome[[j]] <- qLearnEst(moMain = mods$moMain,
                                moCont = mods$moCont,
                                data = data,
                                response = YTilde(outcome[[j+1L]]),
                                iter = iter,
                                txInfo = txInfo[[j]])

      j <- j - 1L
    }
  }

  if( length(outcome) < 0.5 ) {
    outcome <- NULL
  } else if( as.integer(round(length(outcome),0L)) == 1L ) {
    outcome <- outcome[[1]]
  } else {
    outcome <- new("QLearnEstList",
                   loo = outcome)
  }

  return(list(   "outcome" = outcome, 
              "propensity" = propen))

}
