#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# fitChoose - Choose appropriate fitting algorithm                             #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# moMain  : an object of class modelObj or class modelObjSubset that defines   #
#           the model and R methods to be used to obtain parameter estimates   #
#           and predictions for main effects component of an outcome           #
#           regression step.                                                   #
#           See ?modelObj and/or ?modelObjSubset for details.                  #
#           NULL is an appropriate value if iter <=0 and moCont is defined.    #
#                                                                              #
# moCont  : an object of class modelObj or class modelObjSubset that defines   #
#           the model and R methods to be used to obtain parameter estimates   #
#           and predictions for contrasts component of an outcome              #
#           regression step.                                                   #
#           See ?modelObj and/or ?modelObjSubset for details.                  #
#           NULL is an appropriate value if iter <=0 and moMain is defined.    #
#                                                                              #
# data    : data.frame of covariates                                           #
#                                                                              #
# response: response vector                                                    #
#                                                                              #
# txName  : column header of data containing treatment variable                #
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
#= Returns an object of class simpleFit or iterateFit                         =#
#=                                                                            =#
#==============================================================================#
fitChoose <- function(moMain,
                      moCont,
                      data, 
                      response, 
                      txName, 
                      iter){


  if( is(moCont, "NULL") || is(moMain, "NULL") ) {

    iter <- 0L

  } else {

    if( iter <= 0L ){
      once <- isTRUE(all.equal(solver(moMain), solver(moCont)))

      if(!once){
        warning(paste("Solver method for main effects and contrasts differ.\n",
                      "Solutions will be obtained using iterative method.", 
                      sep=""))
        iter <- 100L
      }
    }
  }

  #--------------------------------------------------------------------------#
  # call appropriate fit routine                                             #
  #--------------------------------------------------------------------------#
  if( iter <= 0L ) {
    fit <- fitCombined(moMain = moMain, 
                       moCont = moCont, 
                       response = response, 
                       txName = txName, 
                       data = data)
  } else {
    fit <- fitIterate(moMain = moMain,  
                      moCont = moCont, 
                      response = response,  
                      txName = txName,  
                      data = data, 
                      max.iter = iter)
  }

  #--------------------------------------------------------------------------#
  # Return fit object                                                        #
  #--------------------------------------------------------------------------#
  return(fit)
}

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# pickOutcomeModels - extract all models needed for decision point ndp.        #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# moMain   : an object of class modelObj or class modelObjSubset that defines  #
#            the model and R methods to be used to obtain parameter estimates  #
#            and predictions for main effects component of an outcome          #
#            regression step.                                                  #
#            See ?modelObj and/or ?modelObjSubset for details.                 #
#                                                                              #
# moCont   : an object of class modelObj or class modelObjSubset that defines  #
#            the model and R methods to be used to obtain parameter estimates  #
#            and predictions for contrasts component of an outcome             #
#            regression step.                                                  #
#            See ?modelObj and/or ?modelObjSubset for details.                 #
#            NULL is an appropriate value.                                     #
#                                                                              #
# ndp      : decision point                                                    #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns a list of models for moCont and moMain                             =#
#=                                                                            =#
#==============================================================================#

pickOutcomeModels <- function(moMain, 
                              moCont, 
                              ndp){

  if( !is(moMain, "NULL") && is(moCont, "NULL") ){

    moContTemp <- NULL

    if( is(moMain,'modelObj') ){
      if(ndp != 1L) {
        UserError("input","Method requires more than 1 decision point.")
      } else {
        moMainTemp <- moMain
      }

    } else if( is(moMain,'ModelObjList') ){

      if( length(moMain) < ndp ) {
        UserError("input","Insufficient number of models provided in moMain.")
      }

      moMainTemp <- moMain[[ndp]]

    } else if( is(moMain,'ModelObjSubsetList') ) {

      moMainTemp <- list()

      for( i in 1L:length(moMain) ){
        if( DecisionPoint(moMain[[i]]) != ndp ) next

        moMainTemp <- c(moMainTemp, moMain[[i]])
      }

      if( length(moMainTemp) < 0.5 ) {
        UserError("input", 
                  paste("Could not find models for decision point.", 
                        ndp, sep=""))
      }

      moMainTemp <- new("ModelObjSubsetList",
                        loo = moMainTemp)

    } else {

      DeveloperError("Could not identify class of moMain.","pickOutcomeModels")

    }

  } else if( is(moMain, "NULL") && !is(moCont, "NULL") ){

    moMainTemp <- NULL

    if( is(moCont,'modelObj') ){

      if( ndp != 1L ) {
        UserError("input","Method requires more than 1 decision point.")
      } else {
        moContTemp <- moCont
      }

    } else if( is(moCont,'ModelObjList') ){

      if(length(moCont) < ndp ) {
        UserError("input", "Insufficient number of models provided in moCont.")
      }

      moContTemp <- moCont[[ndp]]

    } else if( is(moCont,'ModelObjSubsetList') ) {

      moContTemp <- list()

      for( i in 1L:length(moCont) ){
        if( DecisionPoint(moCont[[i]]) != ndp ) next

        moContTemp <- c(moContTemp, moCont[[i]])
      }

      if( length(moContTemp) < 0.5 ) {
        UserError("input", 
                  paste("Could not find models for decision point.", 
                        ndp, sep=""))
      }

      moContTemp <- new("ModelObjSubsetList",
                        loo = moContTemp)

    } else {
      DeveloperError("Could not identify class of moCont.", "pickOutcomeModels")
    }

  } else if( !is(moMain, "NULL") && !is(moCont, "NULL") ){

    if( is(moMain,'modelObj') && is(moCont,'modelObj') ) {
      if( ndp != 1L ) {
        UserError("input", "Method requires more than 1 decision point.")
      } else {
        moMainTemp <- moMain
        moContTemp <- moCont
      }

    } else if(  is(moMain,'ModelObjList') && is(moCont,'ModelObjList')  ){

      if( (length(moMain) < ndp) || (length(moCont) < ndp) ) {
        UserError("input", "Insufficient number of models provided.")
      }

      if( !isTRUE(all.equal(length(moMain),length(moCont))) ) {
        UserError("input",
                  paste("moMain and moCont object must contain ",
                        "the same number of models.", sep=""))
      }

      moMainTemp <- moMain[[ndp]]
      moContTemp <- moCont[[ndp]]

    } else if(  is(moMain,'ModelObjSubsetList') && 
                is(moCont,'ModelObjSubsetList')  ){

      moMainTemp <- list()
      moContTemp <- list()

      for( i in 1L:length(moMain) ){

        if( DecisionPoint(moMain[[i]]) != ndp ) next

        moMainTemp <- c(moMainTemp, moMain[[i]])
        txg <- Subset(moMain[[i]])

        found <- FALSE

        for( k in 1L:length(moCont)  ) {

          txgC <- Subset(moCont[[k]])

          if( all(txg %in% txgC) && all(txgC %in% txg) ){
            moContTemp <- c(moContTemp, moCont[[k]]) 
            found <- TRUE
          }
        }

        if( !found ){
          UserError("input",
                    "Unable to match subset of moMain to a subset in moCont.")
        }
      }

      if( length(moMainTemp) < 0.5 ) {
        UserError("input",
                  "Unable to identify moMain for specified decision point.")
      }


      moMainTemp <- new("ModelObjSubsetList",
                        loo = moMainTemp)

      moContTemp <- new("ModelObjSubsetList",
                        loo = moContTemp)

    } else {
      DeveloperError("Could not identify class of moMain and moCont objects.",
                     "pickOutcomeModels")
    }
  }

  return(list("moMain" = moMainTemp,
              "moCont" = moContTemp))
}

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# pickPropenModels - extract all models needed for decision point ndp.         #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# moPropen : an object of class modelObj or class modelObjSubset that defines  #
#            the model and R methods to be used to obtain parameter estimates  #
#            and predictions for propensity of treatment.                      #
#            See ?modelObj and/or ?modelObjSubset for details.                 #
#                                                                              #
# ndp      : decision point                                                    #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns a list of models for moPropen                                      =#
#=                                                                            =#
#==============================================================================#
pickPropenModels <- function(moPropen, 
                             ndp){

  if( is(moPropen,'modelObj') ) {

    if( ndp != 1L ) {
      UserError("input", "Method requires more than 1 decision point.")
    } else {
      moPropenTemp <- moPropen
    }

  } else if( is(moPropen,'ModelObjList') ){

    if(length(moPropen) < ndp) {
      UserError("input", "Insufficient number of models in moPropen")
    }

    moPropenTemp <- moPropen[[ndp]]

  } else if( is(moPropen,'ModelObjSubsetList') ) {

    moPropenTemp <- list()
    for( i in 1L:length(moPropen) ){
      if( DecisionPoint(moPropen[[i]]) != ndp) next
      moPropenTemp <- c(moPropenTemp, moPropen[[i]])
    }

    moPropenTemp <- new("ModelObjSubsetList",
                        loo = moPropenTemp)

  } else {

    DeveloperError("Could not identify class of moPropen object.",
                   "pickPropenModels")
  }

  return( moPropenTemp )
}
