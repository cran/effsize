\name{effsize-package}
\alias{effsize-package}
\alias{effsize}
\docType{package}
\title{
Effect size computation
}
\description{
This packages contains functions to compute effect sizes both based on means difference (Cohen's d and Hedges g) and dominance matrices (Cliff's Delta).
}
\details{
\tabular{ll}{
Package: \tab effsize\cr
Type: \tab Package\cr
Version: \tab 0.3.1\cr
Date: \tab 2014-02-25\cr
License: \tab GPL-2\cr
URL:  \tab \url{http://softeng.polito.it/software/effsize/} \cr
      \tab    \url{http://softeng.polito.it/software/effsize/effsize_0.3.1.pdf} \cr
      \tab    \url{http://softeng.polito.it/software/effsize/effsize_0.3.1.tar.gz}
}
The main functions are:

\code{\link{cliff.delta}}.

\code{\link{cohen.d}}.


\bold{Change history}

\describe{
\item{0.3.1}{Fixed a bug in \code{cohen.d} when \code{PAIRED=TRUE}, now the \code{PAIRED} parameter has no effect, it is left just for compatibility. In a future code clean-up it may be removed}
}
}
\author{
Marco Torchiano \url{http://softeng.polito.it/torchiano/}
}
\keyword{ package }
\keyword{ effect size }