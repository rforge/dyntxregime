#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# fitPropen : Fit propensity score models.                                     #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# moPropen      : an object of class modelObj or modelObjSubsetList, which     #
#                 define the models and R method to be used to obtain          #
#                 parameter estimates and predictions for the propensity of    #
#                 treatment.                                                   #
#                 See ?modelObj and/or ?modelObjSubset for details.            #
#                                                                              #
# txInfo        : an object of class txInfo                                    #
#                                                                              #
# data          : full data frame of covariates.                               #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns an object of class propenFit                                       =#
#=                                                                            =#
#==============================================================================#
fitPropen <- function(moPropen,
                      txInfo,
                      data){

  subsets <- Subsets(txInfo)
  ptsSubset <- PtsSubset(txInfo)
  txName <- TxName(txInfo)

  if( is(moPropen, "ModelObjSubsetList") ){

    propenFitObj <- list()

    for( k in 1L:length(moPropen) ){

      #------------------------------------------------------------------#
      # Extract the subset name(s) for the kth model                     #
      #------------------------------------------------------------------#
      modelSubset <- Subset(moPropen[[k]])

      #------------------------------------------------------------------#
      # Identify which fSet subsets match the model subsets              #
      #------------------------------------------------------------------#
      indSubset <- which( names(subsets) %in% modelSubset )

      if( length(indSubset) < 0.5 ) {
        UserError("input", 
                  paste("Unable to match moPropen subset ", modelSubset, 
                        " to a subset defined by the fSet.", sep=""))
      }

      #------------------------------------------------------------------#
      # Determine which pts are associated with the model subsets        #
      #------------------------------------------------------------------#
      use4fit <- ptsSubset %in% modelSubset

      if( sum(use4fit) < 0.5 ) {
        UserError("input", 
                  paste("No observations match moPropen model subset:",
                        paste(modelSubset,collapse=","),sep=" "))
      }

      predArgs <- predictorArgs(moPropen[[k]]@modelObject)

      if( is(predArgs[["propen.missing"]], "NULL") ) {
        small <- TRUE
      } else if( tolower(predArgs[["propen.missing"]]) == "smallest" ) {
        small <- TRUE
        predArgs[["propen.missing"]] <- NULL
      } else if( tolower(predArgs[["propen.missing"]]) == "largest" ) {
        small <- FALSE
        predArgs[["propen.missing"]] <- NULL
      } 
      predictorArgs(moPropen[[k]]@modelObject) <- predArgs
      
      #------------------------------------------------------------------#
      # Subset modeling causes problems with identifiability of tx levels#
      # The following tries to track what is being estimated. If tx is   #
      # a factor, the tx vector for the subset of data is recast as      #
      # another factor variable to get the correct levels. If tx is an   #
      # integer, the unique tx values for the subset of data are sorted  #
      # It is thus assumed that the prediction method with return the    #
      # probability matrix in the default order of factors or in sorted  #
      # order of the integer values.                                     #
      #------------------------------------------------------------------#
      tData <- data[use4fit,,drop=FALSE]

      if( is.factor(data[,txName]) ) {
        tempTx <- factor(data[use4fit,txName])
        levs <- levels(tempTx)
      } else {
        tempTx <- data[use4fit,txName]
        levs <- sort(unique(round(tempTx,0L)))
        levs <- as.character(levs)
      }
      
      tData[,txInfo@txName] <- tempTx
      
      fitResult <- Fit(object = moPropen[[k]], 
                       data = tData, 
                       response = tData[,txName])

      propenFitObj[[k]] <- new("PropenSubsetFit",
                               subset = modelSubset,
                               levels = levs,
                               small = small,
                               modelObjectFit = fitResult)

    }

    propenFitObj <- new("PropenSubsetFitList",
                        loo = propenFitObj,
                        txInfo = txInfo)

  } else if( is(moPropen, "modelObj") ){

    #----------------------------------------------------------------------#
    # Eliminate patients with only 1 tx option from dataset for fit        #
    #----------------------------------------------------------------------#
    useGrps <- sapply(X = subsets, FUN = length) > 1.5
    use4fit <- ptsSubset %in% names(subsets)[useGrps]

    if( sum(use4fit) < 0.5 ) {
      UserError("input",
                "No patients have more than 1 treatment option.")
    }
    
    #----------------------------------------------------------------------#
    # feasibility rules also causes problems w/ identifiability of levels  #
    # The following tries to track what is being estimated. If tx is       #
    # a factor, the tx vector for the subset of data is recast as          #
    # another factor variable to get the correct levels. If tx is an       #
    # integer, the unique tx values for the subset of data are sorted      #
    # It is thus assumed that the prediction method with return the        #
    # probability matrix in the default order of factors or in sorted      #
    # order of the integer values.                                         #
    #----------------------------------------------------------------------#
    if( is.factor(data[,txName]) ) {
      levs <- factor(data[use4fit,txName])
      levs <- levels(levs)
    } else {
      levs <- sort(unique(round(data[use4fit,txName],0L)))
      levs <- as.character(levs)
    }

    predArgs <- predictorArgs(moPropen)
    if( is(predArgs[["propen.missing"]], "NULL") ) {
      small <- TRUE
    } else if( tolower(predArgs[["propen.missing"]]) == "smallest" ) {
      small <- TRUE
      predArgs[["propen.missing"]] <- NULL
    } else if( tolower(predArgs[["propen.missing"]]) == "largest" ) {
      small <- FALSE
      predArgs[["propen.missing"]] <- NULL
    } 

    predictorArgs(moPropen) <- predArgs

    #----------------------------------------------------------------------#
    # fit propen model using provided solver                               #
    #----------------------------------------------------------------------#
    propenFitObj <- fit(object = moPropen, 
                        data = data[use4fit,,drop=FALSE], 
                        response = data[use4fit,txName])

    propenFitObj <- new("PropenModelObjFit",
                        modelObjectFit = propenFitObj,
                        levels = levs,
                        small = small)
  } else {

    DeveloperError("unrecognized class for moPropen", "fitPropen")

  }

  result <- new("PropenFit",
                fits = propenFitObj)

  return(result)
}


#------------------------------------------------------------------------------#
# prob_matrix : determine propensity for treatment matrix.                     #
#------------------------------------------------------------------------------#
#                                                                              #
# fitObj : an object of class modelObjFit                                      #
#                                                                              #
# data   : full data frame of covariates.                                      #
#                                                                              #
# txName : a character indicating the column of data containing tx             #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= Returns a matrix of probabilities for all treatment categories of txg      =#
#=                                                                            =#
#==============================================================================#
prob_matrix <- function(fitObj, 
                        data, 
                        txName,
                        sset,
                        msg = FALSE) {

  tol <- 1.5e-8
  
  #--------------------------------------------------------------------------#
  # Obtain predictions for entire dataset                                    #
  #--------------------------------------------------------------------------#
  mm <- PredictPropen(object = fitObj, 
                      newdata = data, 
                      txName = txName, 
                      sset = sset)

  #--------------------------------------------------------------------------#
  # Verify that treatment probabilities are not negative                     #
  #--------------------------------------------------------------------------#
  if( any(mm < -tol) ) {
    UserError("input",
              "Treatment probabilities are negative. Verify moPropen inputs.")
  }

  #--------------------------------------------------------------------------#
  # Verify that the total probability for each patient is <= 1               #
  #--------------------------------------------------------------------------#
  if( sum(rowSums(mm) > {1.0 + tol}) > 0.5 ) {
    UserError("input",
              "Sum of treatment probabilities >1. Verify moPropen inputs.")
  }

  return(mm)
}
