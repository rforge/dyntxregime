#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                              CLASS IterateFit                                #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# contains results obtained by a call to modelObjFit with separate models      #
#                                                                              #
#   modelObjectFitMain : An object of class modelObjFit for main effects model #
#                                                                              #
#   modelObjectFitCont : An object of class modelObjFit for contrast model     #
#                                                                              #
#   residuals          : residuals of the combined fit                         #
#                                                                              #
#   yMainHat           : fitted main effects                                   #
#                                                                              #
#   yContHat           : fitted contrast                                       #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass(Class = "IterateFit", 
         slots = c(baseLevel = "character or numeric",
                   txName = "character",
                   modelObjectFitMain = "modelObjFit",
                   modelObjectFitCont = "modelObjFit",
                   residuals = "numeric",
                   yMainHat = "numeric", 
                   yContHat = "numeric"))

setClass(Class = "IterateFitList",
         contains = "dpList")

setMethod(f = "Base",    
          signature = c(object = "IterateFit"), 
          definition = function(object, ...){ 
                         return( object@baseLevel ) 
                       } )

setMethod(f = "Coef",
          signature = c(object = "IterateFit"),
          definition = function(object, ...) {
                         me <- coef(object@modelObjectFitMain, ...)
                         cn <- coef(object@modelObjectFitCont, ...)
                         res <- list("MainEffect" = me,
                                     "Contrast" = cn)
                         return( res )
                       } )

setMethod(f = "FitObject", 
          signature = c(object = "IterateFit"), 
          definition = function(object, ...){
                         me <- fitObject(object@modelObjectFitMain, ...)
                         cn <- fitObject(object@modelObjectFitCont, ...)
                         res <- list("MainEffect" = me,
                                     "Contrast" = cn)
                         return( res )
                       } )

setMethod(f = "FittedCont", 
          signature = c(object = "IterateFit"), 
          definition = function(object, ...){ return( object@yContHat ) } )

setMethod(f = "FittedMain", 
          signature = c(object = "IterateFit"), 
          definition = function(object, ...){ return( object@yMainHat ) } )

setMethod(f = "ModelObjectFit", 
          signature = c(object = "IterateFit"), 
          definition = function(object, ...){
                         me <- object@modelObjectFitMain
                         cn <- object@modelObjectFitCont
                         res <- list("MainEffect" = me,
                                     "Contrast" = cn)
                         return( res )
                       } )

setMethod(f = "MySummary",
          signature = c(object = "IterateFit"),
          definition = function(object, ...) {
                         me <- MySummary(object@modelObjectFitMain, ...)
                         cn <- MySummary(object@modelObjectFitCont, ...)
                         res <- list("MainEffect" = me,
                                     "Contrast" = cn)
                         return( res )
                       } )

setMethod(f = "Plot", 
          signature = c(x = "IterateFit"), 
          definition = function(x, suppress=FALSE, ...){

                         argList <- list(...)
                         if( !suppress ) {
                           if( is(argList[[ "main" ]], "NULL") ) {
                             argList[[ "main" ]] <- "moMain"
                           } else if( is(argList[[ "sub" ]], "NULL") ) {
                             argList[[ "sub" ]] <- "moMain"
                           } else {
                             argList[[ "sub" ]] <- paste(argList[[ "sub" ]], 
                                                         " (moMain)", sep="")
                           }
                         }
                         argList[[ "x" ]] <- x@modelObjectFitMain
                         do.call(what = plot, args = argList)

                         argList <- list(...)
                         if( !suppress ) {
                           if( is(argList[[ "main" ]], "NULL") ) {
                             argList[[ "main" ]] <- "moCont"
                           } else if( is(argList[[ "sub" ]], "NULL") ) {
                             argList[[ "sub" ]] <- "moCont"
                           } else {
                             argList[[ "sub" ]] <- paste(argList[[ "sub" ]], 
                                                         " (moCont)", sep="")
                           }
                         }
                         argList[[ "x" ]] <- x@modelObjectFitCont
                         do.call(what = plot, args = argList)
                       } )


setMethod(f = "PredictCont", 
          signature = c(object = "IterateFit", 
                        newdata = "data.frame"), 
          definition = function(object, newdata, ...){

                         if( is(object@baseLevel, "character") ) {

                           newdata[,object@txName] <- 
                             as.factor(newdata[,object@txName])

                         } else {
                           newdata[,object@txName] <- as.integer(newdata[,object@txName])
                         }

                         main <- predict(object = object@modelObjectFitMain,  
                                         newdata = newdata)

                         contrast <- predict(object = object@modelObjectFitCont, 
                                             newdata = newdata)

                         fittedY <- main + contrast

                         n <- nrow(newdata)

                         if( is(object@baseLevel, "character") ) {

                           newdata[,object@txName] <- 
                             as.factor(rep(object@baseLevel,n))

                         } else {
                           newdata[,object@txName] <- object@baseLevel
                         }

                         mainBase <- predict(object = object@modelObjectFitMain,  
                                             newdata = newdata)

                         contrastBase <- predict(object = object@modelObjectFitCont, 
                                                 newdata = newdata)

                         fittedBase <- mainBase + contrastBase

                         fittedCont <- fittedY - fittedBase

                         return( fittedCont )
                       } )

setMethod(f = "PredictMain", 
          signature = c(object = "IterateFit", 
                        newdata = "data.frame"), 
          definition = function(object, newdata, ...){

                         n <- nrow(newdata)

                         if( is(object@baseLevel, "character") ) {
                           newdata[,object@txName] <- 
                             as.factor(rep(object@baseLevel,n))
                         } else {
                           newdata[,object@txName] <- object@baseLevel
                         }

                         mainBase <- predict(object = object@modelObjectFitMain,  
                                             newdata = newdata)

                         contrastBase <- predict(object = object@modelObjectFitCont, 
                                                 newdata = newdata)

                         fittedBase <- mainBase + contrastBase

                         return( fittedBase )
                       } )

setMethod(f = "Print", 
          signature = c(x = "IterateFit"), 
          definition = function(x, ...){
                         cat("\n *** moMain Fit *** \n")
                         show(x@modelObjectFitMain)
                         cat("\n *** moCont Fit *** \n")
                         show(x@modelObjectFitCont)
                       } )

setMethod(f = "Residuals", 
          signature = c(object = "IterateFit"), 
          definition = function(object, ...){ 
                         return( object@residuals ) 
                       } )

setMethod(f = "Show", 
          signature = c(object = "IterateFit"), 
          definition = function(object, ...){
                         cat("\n *** moMain Fit *** \n")
                         show(object@modelObjectFitMain)
                         cat("\n *** moCont Fit *** \n")
                         show(object@modelObjectFitCont)
                       } )

if(!isClass("SimpleFit or IterateFit")){
  setClassUnion("SimpleFit or IterateFit", 
                members = c("SimpleFit","IterateFit"))
}


