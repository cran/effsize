\name{effsize-package}
\alias{effsize-package}
\alias{effsize}
\docType{package}
\title{
Effect size computation
}
\description{
This packages contains functions to compute effect sizes both based on means difference (Cohen's d and Hedges g) and dominance matrices (Cliff's Delta).

The computation (especially for Cliff's Delta) is carried on with higly efficient algorithms.
}
\details{
\tabular{ll}{
Package: \tab effsize\cr
Type: \tab Package\cr
Version: \tab 0.4.1\cr
Date: \tab 2014-02-25\cr
License: \tab GPL-2\cr
URL:  \tab \url{http://softeng.polito.it/software/effsize/} \cr
      \tab    \url{http://softeng.polito.it/software/effsize/effsize_0.4.1.pdf} \cr
      \tab    \url{http://softeng.polito.it/software/effsize/effsize_0.4.1.tar.gz}
}
The main functions are:

\code{\link{cliff.delta}}.

\code{\link{cohen.d}}.


\bold{Change history}

\describe{
\item{0.3.1}{Fixed a bug in \code{cohen.d} when \code{PAIRED=TRUE}, now the \code{PAIRED} parameter has no effect, it is left just for compatibility. In a future code clean-up it may be removed}
\item{0.4}{Implemented a new algorithm with improved memory and time complexity. In particular new time complexity is T = O(n1*log(n2)) vs. the previous T = O(n1*n2), and new memory complexity M = O( n1 + n2 ) vs. the previous M = O( n1 * n2). In practice now the computation becomes feasible in a "\emph{reasonable}" time.}
\item{0.4.1}{Code clean-up and optimization using vectorized binary partioning.}
}
}
\author{
Marco Torchiano \url{http://softeng.polito.it/torchiano/}
}
\keyword{ package }
\keyword{ effect size }