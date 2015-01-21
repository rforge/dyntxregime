#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# iqLearn_pm : Calculate value functions for tx = +1/-1 for new patient        #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# object : object of class iqLearnSS, iqLearnFS_ME, iqLearnFS_C, or            #
#          iqLearnFS_VHet                                                      #
#                                                                              #
# newdata : data.frame of covariate information for new patient                #
#                                                                              #
#==============================================================================#
#=                                                                            =#
#= returns a list :                                                           =#
#=  $pos : value function if positive tx option given                         =#
#=  $neg : value function if negative tx option given                         =#
#=                                                                            =#
#==============================================================================#
iqLearn_pm <- function(object, 
                       newdata){

  n <- nrow(newdata)

  data <- rbind(newdata, newdata)

  tx <- c(rep(-1L,n),rep(1L,n))

  if( TxName(object) %in% colnames(newdata) ){
    data[,TxName(object)] <- tx
  } else {
    data <- cbind(data,tx)
    colnames(data) <- c(colnames(newdata),TxName(object))
  }
  #--------------------------------------------------------------------------#
  # Calculate value function for positive treatment                          #
  # Note that predict method has to be called because combined fits cannot   #
  # be broken down into "main" and "contrast"; "contrast" will be sent back  #
  # as zero.                                                                 #
  # For iterate fits, the formula has been modified to include the treatment #
  # variable, therefore the value of treatment does not need to be considered#
  #--------------------------------------------------------------------------#
  mn <- PredictMain(object = object, newdata = data)
  cn <- PredictCont(object = object, newdata = data)

  #--------------------------------------------------------------------------#
  # Return positive and negative value functions                             #
  #--------------------------------------------------------------------------#
  return( matrix(data = mn + cn, ncol = 2L) )
}

