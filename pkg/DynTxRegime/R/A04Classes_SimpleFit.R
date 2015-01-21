#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                               CLASS SimpleFit                                #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# contains results obtained by a call to modelObjFit with a combined model     #
#                                                                              #
#   fitType        : character indicating ME, C, or ME+C models                #
#                                                                              #
#   modelObjectFit : An object of class modelObjFit                            #
#                                                                              #
#   residuals      : residuals of the combined fit                             #
#                                                                              #
#   yContHat       : fitted contrast                                           #
#                                                                              #
#   yMainHat       : fitted main effect                                        #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
if(!isClass("character or numeric")){
  setClassUnion("character or numeric", 
                members = c("character","numeric"))
}

setClass(Class = "SimpleFit", 
         slots = c(     baseLevel = "character or numeric",
                           txName = "character",
                          fitType = "character",
                   modelObjectFit = "modelObjFit",
                        residuals = "numeric",
                         yContHat = "numeric",
                         yMainHat = "numeric"))

setClass(Class = "SimpleFitList",
         contains = "dpList")

setMethod(f = "Base",    
          signature = c(object = "SimpleFit"), 
          definition = function(object, ...){ 
                         return( object@baseLevel ) 
                       } )

setMethod(f = "Coef", 
          signature = c(object = "SimpleFit"), 
          definition = function(object, ...){
                         res <- list()
                         nms <- object@fitType
                         res[[ nms ]] <- coef(object@modelObjectFit, ...)
                         return( res )
                       } )

setMethod(f = "FitObject", 
          signature = c(object = "SimpleFit"), 
          definition = function(object, ...){
                         res <- list()
                         nms <- object@fitType
                         res[[ nms ]] <- fitObject(object@modelObjectFit, ...)
                         return( res )
                       } )

setMethod(f = "FittedCont", 
          signature = c(object = "SimpleFit"), 
          definition = function(object, ...){ return( object@yContHat ) } )

setMethod(f = "FittedMain", 
          signature = c(object = "SimpleFit"), 
          definition = function(object, ...){ return( object@yMainHat ) } )

setMethod(f = "FitType",
          signature = c(object = "SimpleFit"),
          definition = function(object, ...){ return( object@fitType ) } )

setMethod(f = "ModelObjectFit", 
          signature = c(object = "SimpleFit"), 
          definition = function(object, ...){ return( object@modelObjectFit ) } )

setMethod(f = "MySummary", 
          signature = c(object = "SimpleFit"), 
          definition = function(object, ...){
                         res <- list()
                         nms <- object@fitType
                         res[[ nms ]] <- summary(object@modelObjectFit, ...)
                         return( res )
                       } )

setMethod(f = "Plot", 
          signature = c(x = "SimpleFit"), 
          definition = function(x, suppress=FALSE, ...){

                         argList <- list(...)

                         nms <- x@fitType

                         if( !suppress ) {
                           if( is(argList[[ "main" ]], "NULL") ) {
                             argList[[ "main" ]] <- nms
                           } else if( is(argList[[ "sub" ]], "NULL") ) {
                             argList[[ "sub" ]] <- nms
                           } else {
                             argList[[ "sub" ]] <- paste(argList[[ "sub" ]], 
                                                         " (", nms, ")", sep="")
                           }
                         }
                         argList[[ "x" ]] <- x@modelObjectFit
                         do.call(what = plot, args = argList)
                       } )

setMethod(f = "PredictCont", 
          signature = c(object = "SimpleFit",  
                        newdata = "data.frame"), 
          definition = function(object, newdata, ...){

                         fittedY <- predict(object = object@modelObjectFit,  
                                            newdata = newdata)

                         n <- nrow(newdata)

                         if( is(object@baseLevel, "character") ) {
                           newdata[,object@txName] <- 
                             factor(rep(object@baseLevel,n), 
                                    levels = levels(newdata[,object@txName]))
                         } else {
                           newdata[,object@txName] <- object@baseLevel
                         }

                         baseY <- predict(object = object@modelObjectFit,  
                                          newdata = newdata)

                         fittedCont <- fittedY - baseY

                         return( fittedCont )

                       } )

setMethod(f = "PredictMain", 
          signature = c(object = "SimpleFit", 
                        newdata = "data.frame"), 
          definition = function(object, newdata, ...){

                         n <- nrow(newdata)

                         if( is(object@baseLevel, "character") ) {
                           newdata[,object@txName] <- 
                             factor(rep(object@baseLevel,n), 
                                    levels = levels(newdata[,object@txName]))
                         } else {
                           newdata[,object@txName] <- object@baseLevel
                         }

                         baseY <- predict(object = object@modelObjectFit,  
                                          newdata = newdata)

                         return( baseY )
                       } )

setMethod(f = "Print",
          signature = c(x = "SimpleFit"),
          definition = function(x, ...){
                         cat("\n *** ", x@fitType, " Fit ***\n", sep="")
                         print(x@modelObjectFit)
                        } )

setMethod(f = "Residuals", 
          signature = c(object = "SimpleFit"), 
          definition = function(object, ...){ 
                         return( object@residuals ) 
                        } )

setMethod(f = "Show", 
          signature = c(object = "SimpleFit"), 
          definition = function(object, ...){
                         cat("\n *** ", object@fitType, " Fit ***\n", sep="")
                         show(object@modelObjectFit)
                       } )


