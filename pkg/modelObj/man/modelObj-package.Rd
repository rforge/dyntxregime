\name{modelObj-package}
\alias{modelObj-package}
\docType{package}
\title{A Model Object Framework for Regression Analysis}
\description{
A utility package to facilitate the generalization of methods built on a regression framework.
}
\details{
\tabular{ll}{
Package: \tab modelObj\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-02-23\cr
License: \tab GPL-2\cr
Depends: \tab methods\cr
}
Often, new statistical methods are developed on the framework of
traditional regression or classification methods. To simplify the
creation of new R implementations of these methods, 
researchers and software developers
often make choices regarding the types of models that can be
used by the user; hard-coding the regression method into the library and
limiting or eliminating the ability of the user to modify regression
control parameters. These choices artificially limit the general application of
new methods. In addition, if a new method requires multiple models, 
a developer is
often forced to artificially break the method into multiple function calls, each
for a specific regression/classification step, or is forced to provide a
cumbersome and/or confusing interface for the user. \pkg{modelObj} is
an R package developed to facilitate the use of existing and
future R regression and classification libraries that simplifies 
the development of general, non-model-specific implementations of new 
statistical methods. 
}
\author{
Shannon T. Holloway \cr
Maintainer: Shannon T. Holloway <sthollow@ncsu.edu>
}
\keyword{ package }

