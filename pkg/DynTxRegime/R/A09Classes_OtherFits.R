#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                            CLASS PropenSubsetFit                             #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# subset         : subset of treatment for which the model is to be used       #
#                                                                              #
# small          : an object of class logical indicating missing tx            #
#                                                                              #
# modelObjectFit : an object of class modelObjFit                              #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass("PropenSubsetFit",
         slots = c(        subset = "character",
                            small = "logical",
                   modelObjectFit = "modelObjFit"))

setClass(Class = "PropenSubsetFitList",
         slots = c(txInfo = "TxInfo"),
         contains = "List" )

setMethod(f = "Coef",
          signature = c(object = "PropenSubsetFit"), 
          definition = function(object, ...){
                         return( Coef(object = object@modelObjectFit, ...) ) 
                       } )

setMethod(f = "FitObject",
          signature = c(object = "PropenSubsetFit"), 
          definition = function(object, ...){
                         return( FitObject(object = object@modelObjectFit) )
                       } )

setMethod(f = "Fitted",
          signature = c(object = "PropenSubsetFit"), 
          definition = function(object, ...){
                         return( predict(object = object@modelObjectFit, ...) )
                       } )

setMethod(f = "ModelObjectFit",
          signature = c(object = "PropenSubsetFit"), 
          definition = function(object, ...){ 
                         return( object@modelObjectFit ) 
                       } )

setMethod(f = "MySummary",
          signature = c(object = "PropenSubsetFit"), 
          definition = function(object, ...){
                         return( summary(object = object@modelObjectFit, ...) ) 
                       } )

setMethod(f = "Plot",
          signature = c(x = "PropenSubsetFit"), 
          definition = function(x, ...){
                         plot(x = x@modelObjectFit, ...) 
                       } )

setMethod(f = "Predict",
          signature = c(object = "PropenSubsetFit",
                        newdata = "data.frame"), 
          definition = function(object, newdata, ...){

                         mm <- predict(object = object@modelObjectFit, 
                                       newdata = newdata, ...) 

                         return( mm )
                       } )

setMethod(f = "PredictPropen", 
          signature = c(object = "PropenSubsetFit",
                        newdata = "data.frame"),
          definition = function(object, newdata, txName, ...){

                         mm <- predict(object = object@modelObjectFit, 
                                       newdata=newdata, ...)

                         nfv <- ncol(mm)

                         if( is(newdata[,txName], "factor") ) {
                           nlevs <- length(levels(newdata[,txName]))
                           nms <- levels(newdata[,txName])
                         } else {
                           nlevs <- length(unique(newdata[,txName]))
                           nms <- as.character(sort(unique(newdata[,txName])))
                         }

                         if( nfv == {nlevs-1L} ) {

                           correction <- 1.0 - rowSums(mm)

                           if( object@small ) {                           
                             mm <- cbind(correction, mm)
                           } else {                           
                             mm <- cbind(mm, correction)
                           }

                         } else if( nfv != nlevs ) {
                           msg <- paste("Number of tx options in data ",
                                        "and/or fSet do not match ",
                                        "model predictions from moPropen ",
                                        "object.",sep="")
                           UserError("input", msg)
                           
                         }

                         colnames(mm) <- nms

                         return( mm )
                       } )

setMethod(f = "Residuals",
          signature = c(object = "PropenSubsetFit"), 
          definition = function(object, ...){
                         return(Residuals(object = object@modelObjectFit, ...))
                       } )

setMethod(f = "Show", 
          signature = c(object = "PropenSubsetFit"), 
          definition = function(object, ...){
                         Show(object@modelObjectFit)
                       } )

setMethod(f = "Subset",
          signature = c(object = "PropenSubsetFit"), 
          definition = function(object, ...){ 
                         return( object@subset ) 
                       } )

setMethod(f = "Coef",
          signature = c(object = "PropenSubsetFitList"),
          definition = function(object, ...){

                         nFits <- length(object)

                         res <- list()

                         for( i in 1L:nFits ) {
                           res[[i]] <- Coef(object = object[[i]], ...) 
                         }

                         return(res)

                       } )

setMethod(f = "Fitted",
          signature = c(object = "PropenSubsetFitList"),
          definition = function(object, ...){

                         nFits <- length(object)

                         ptsSubset <- PtsSubset(object@txInfo)

                         res <- numeric(length(ptsSubset))

                         for( i in 1L:nFits ) {
                           u4f <- ptsSubset %in% Subset(object[[i]])
                           res[u4f] <- Fitted(object = object[[i]], ...) 
                         }

                         return(res)
                       } )

setMethod(f = "Plot",
          signature = c(x = "PropenSubsetFitList"),
          definition = function(x, ...){
                         nFits <- length(x)
                         for( i in 1L:nFits ) {
                           Plot(x = x[[i]], ...)
                         }
                       } )

setMethod(f = "Predict",
          signature = c(object = "PropenSubsetFitList",
                        newdata = "data.frame"),
          definition = function(object, newdata, ...){

                         txInfo <- object@txInfo

                         fs <- feasibility(superSet = SuperSet(txInfo), 
                                           fSet = SubsetRule(txInfo), 
                                           txName = TxName(txInfo),
                                           data = newdata)

                         newPtSubset <- fs$ptsSubset

                         nFits <- length(object)

                         res <- numeric(length(newPtSubset))

                         res <- NULL

                         for( i in 1L:nFits ) {
                           u4f <- newPtSubset %in% Subset(object[[i]])

                           if( sum(u4f) == 0L ) next

                           temp <- Predict(object = object[[i]],
                                           newdata = newdata[u4f,,drop=FALSE])

                           if( is(res, "NULL") ){
                             res <- matrix(data = 0.0, 
                                           nrow = nrow(newdata),
                                           ncol = ncol(temp) )
                           }
                           res[u4f,] <- temp
                         }

                         return(res)

                       } )


setMethod(f = "PredictPropen",
          signature = c(object = "PropenSubsetFitList",
                        newdata = "data.frame"),
          definition = function(object, newdata, txName, ...){

                         txInfo <- object@txInfo

                         fs <- feasibility(superSet = SuperSet(txInfo), 
                                           fSet = SubsetRule(txInfo), 
                                           txName = TxName(txInfo),
                                           data = newdata)

                         newPtSubset <- fs$ptsSubset

                         nFits <- length(object)

                         if( is(newdata[,txName], "factor") ) {
                           nlevs <- length(levels(newdata[,txName]))
                           nms <- levels(newdata[,txName])
                         } else {
                           nlevs <- length(unique(newdata[,txName]))
                           nms <- as.character(sort(unique(newdata[,txName])))
                         }

                         res <- matrix(data = 0.0, 
                                       nrow = nrow(newdata),
                                       ncol = nlevs,
                                       dimnames = list(NULL,nms) )

                         for( i in 1L:nFits ) {
                           u4f <- newPtSubset %in% Subset(object[[i]])

                           if( sum(u4f) == 0L ) next

                           temp <- PredictPropen(object = object[[i]],
                                                 newdata = newdata[u4f,,drop=FALSE],
                                                 txName = txName)
                           nms <- colnames(temp)
                           res[u4f,nms] <- temp
                         }

                         return(res)

                       } )

setMethod(f = "Residuals",
          signature = c(object = "PropenSubsetFitList"),
          definition = function(object, ...){

                         nFits <- length(object)

                         ptsSubset <- PtsSubset(object@txInfo)

                         res <- numeric(length(ptsSubset))

                         for( i in 1L:nFits ) {
                           u4f <- ptsSubset %in% Subset(object[[i]])
                           res[u4f] <- Residuals(object = object[[i]], ...) 
                         }

                         return(res)
                       } )

setMethod(f = "Show",
          signature = c(object = "PropenSubsetFitList"),
          definition = function(object, ...){
                         nFits <- length(object)
                         for( i in 1L:nFits ) {
                           cat("\nSubset: ", object[[i]]@subset, "\n", sep="")
                           Show(object[[i]], ...)
                         }
                       } )

setMethod(f = "Print",
          signature = c(x = "PropenSubsetFitList"),
          definition = function(x, ...){
                         nFits <- length(x)
                         for( i in 1L:nFits ) {
                           cat("\nSubset: ", x[[i]]@subset, "\n", sep="")
                           Print(x[[i]], ...)
                         }
                       } )

setMethod(f = "MySummary",
          signature = c(object = "PropenSubsetFitList"),
          definition = function(object, ...){

                         nFits <- length(object)

                         res <- list()

                         for( i in 1L:nFits ) {
                           res[[i]] <- MySummary(object = object[[i]], ...) 
                         }

                         return(res)

                       } )

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                           CLASS PropenModelObjFit                            #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# subset         : subset of treatment for which the model is to be used       #
#                                                                              #
# small          : an object of class logical indicating missing tx            #
#                                                                              #
# modelObjectFit : an object of class modelObjFit                              #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass("PropenModelObjFit",
         slots = c(         small = "logical",
                   modelObjectFit = "modelObjFit"))

setMethod(f = "Coef",
          signature = c(object = "PropenModelObjFit"), 
          definition = function(object, ...){
                         return( Coef(object = object@modelObjectFit, ...) ) 
                       } )

setMethod(f = "FitObject",
          signature = c(object = "PropenModelObjFit"), 
          definition = function(object, ...){
                         return( FitObject(object = object@modelObjectFit) )
                       } )

setMethod(f = "Fitted",
          signature = c(object = "PropenModelObjFit"), 
          definition = function(object, ...){
                         return( predict(object = object@modelObjectFit, ...) )
                       } )

setMethod(f = "ModelObjectFit",
          signature = c(object = "PropenModelObjFit"), 
          definition = function(object, ...){ 
                         return( object@modelObjectFit ) 
                       } )

setMethod(f = "MySummary",
          signature = c(object = "PropenModelObjFit"), 
          definition = function(object, ...){
                         return( summary(object = object@modelObjectFit, ...) ) 
                       } )

setMethod(f = "Plot",
          signature = c(x = "PropenModelObjFit"), 
          definition = function(x, ...){
                         plot(x = x@modelObjectFit, ...) 
                       } )

setMethod(f = "Predict",
          signature = c(object = "PropenModelObjFit",
                        newdata = "data.frame"), 
          definition = function(object, newdata, ...){

                         mm <- predict(object = object@modelObjectFit, 
                                       newdata = newdata, ...) 

                         return( mm )
                       } )

setMethod(f = "PredictPropen", 
          signature = c(object = "PropenModelObjFit",
                        newdata = "data.frame"),
          definition = function(object, newdata, txName, sset, ...){

                         mm <- predict(object = object@modelObjectFit, 
                                       newdata = newdata, ...)

                         nfv <- ncol(mm)

                         if( is(newdata[,txName], "factor") ) {
                           nlevs <- length(levels(newdata[,txName]))
                           nms <- levels(newdata[,txName])
                         } else {
                           nlevs <- length(sset)
                           nms <- sset
                         }

                         if( nfv == {nlevs-1L} ) {

                           correction <- 1.0 - rowSums(mm)

                           if( object@small ) {                           
                             mm <- cbind(correction, mm)
                           } else {                           
                             mm <- cbind(mm, correction)
                           }

                         } else if( nfv != nlevs ) {
                           msg <- paste("Number of tx options in data ",
                                        "and/or fSet do not match ",
                                        "model predictions from moPropen ",
                                        "object.",nfv,txName,nlevs,sep=" ")
                           UserError("input", msg)
                           
                         }

                         colnames(mm) <- nms

                         return( mm )
                       } )

setMethod(f = "Print", 
          signature = c(x = "PropenModelObjFit"), 
          definition = function(x, ...){
                         Print(x@modelObjectFit)
                       } )

setMethod(f = "Residuals",
          signature = c(object = "PropenModelObjFit"), 
          definition = function(object, ...){
                         return(Residuals(object = object@modelObjectFit, ...))
                       } )

setMethod(f = "Show", 
          signature = c(object = "PropenModelObjFit"), 
          definition = function(object, ...){
                         Show(object@modelObjectFit)
                       } )

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                                CLASS PropenFit                               #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# contains results for the propensity for treatment for a single dp            #
#                                                                              #
#   fits  : If only one the superset of treatment is modeled, an object of     #
#           class modelObjFit. If treatment subsets are modeled individually,  #
#           an object of class PropenSubsetFitList                           #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
if(!isClass("PropenModelObjFit or PropenSubsetFitList")){
  setClassUnion("PropenModelObjFit or PropenSubsetFitList", 
                members = c("PropenModelObjFit","PropenSubsetFitList"))
}

setClass("PropenFit",
         slots = c( fits = "PropenModelObjFit or PropenSubsetFitList") )

setClass(Class = "PropenFitList",
         contains = "dpList" )

setMethod(f = "Coef", 
          signature = c(object="PropenFit"), 
          definition = function(object, ...){ 
                         return( Coef(object = object@fits, ...)) } )

setMethod(f = "FitObject", 
          signature = c(object="PropenFit"), 
          definition = function(object, ...){
                         return( FitObject(object = object@fits, ...) )
                       } )
          
setMethod(f = "ModelObjectFit", 
          signature = c(object = "PropenFit"), 
          definition = function(object, ...){
                         return( ModelObjectFit(object = object@fits, ...) )
                       } )

setMethod(f = "Plot", 
          signature = c(x="PropenFit"), 
          definition = function(x, suppress=FALSE, ...){
                         Plot(x = x@fits, suppress = suppress, ...) 
                       } )

setMethod(f = "Predict", 
          signature = c(object="PropenFit"), 
          definition = function(object, newdata, ...){ 
                         return( Predict(object = object@fits, newdata=newdata, ...) )
                       } )

setMethod(f = "PredictPropen", 
          signature = c(object="PropenFit"), 
          definition = function(object, newdata, txName, ...){ 
                         return( PredictPropen(object = object@fits, 
                                               newdata = newdata, 
                                               txName = txName, ...) )
                       } )

setMethod(f = "Print", 
          signature = c(x = "PropenFit"), 
          definition = function(x, ...){
                         Print(x = x@fits, ...)} )

setMethod(f = "Residuals", 
          signature = c(object="PropenFit"), 
          definition = function(object, ...){ 
                         return( Residuals(object = object@fits, ...) )
                       } )

setMethod(f = "Show", 
          signature = c(object = "PropenFit"), 
          definition = function(object, ...){
                         Show(object = object@fits, ...)} )

setMethod(f = "MySummary", 
          signature = c(object="PropenFit"), 
          definition = function(object, ...){
                         return( MySummary(object = object@fits, ...) )
                       } )

if(!isClass("PropenFit or PropenFitList")){
  setClassUnion("PropenFit or PropenFitList", 
                members = c("PropenFit", "PropenFitList")
  )
}


