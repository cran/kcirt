\name{kcirt-package}
\alias{kcirt-package}
\alias{kcirt}
\docType{package}
\title{
k-Cube Thurstonian IRT Models
}
\description{
Create, Simulate, Fit, Solve k-Cube Thurstonian IRT Models.
}
\details{
\tabular{ll}{
Package: \tab kcirt\cr
Type: \tab Package\cr
Version: \tab 0.5.2\cr
Date: \tab 2013-10-10\cr
License: \tab GPL (>= 2)\cr
}
Use \code{\link{kcirt.model}} to define a k-Cube Thurstonian IRT model. The function \code{\link{kcirt.sim}} generates a random realization.  The function \code{\link{kcirt.fitEE}} uses an expectation-expectation volley to approximately locate \code{mu} and \code{Lambda} and predict the states, \code{Eta}.  The function \code{\link{kcirt.fitMSS}} makes use of metaheuristic stochastic search to further refine the predictions/estimates.


The system of interest is defined as


\bold{\eqn{y}}\eqn{_i}* =  \bold{\eqn{\Delta}}  \bold{\eqn{\mu}}  + \bold{\eqn{\Delta}} \bold{\eqn{\Lambda}} \bold{\eqn{S}} \bold{\eqn{\eta}}\eqn{_i} + \bold{\eqn{\Delta}} \bold{\eqn{\epsilon}}\eqn{_i}

\bold{\eqn{y}} = 1, if \bold{\eqn{y}}* > 0

\bold{\eqn{y}} = 0, otherwise


\bold{\eqn{Y}} = (\eqn{y_1}, \eqn{y_2}, ..., \eqn{y_N})

where

\bold{\eqn{y}}\eqn{_i} is the (column) response vector for observation \eqn{i}.

\bold{\eqn{\Delta}} is the Delta matrix.

\bold{\eqn{\mu}} is the column vector of item means (aka, 'utilities').

\bold{\eqn{\Lambda}} is the hyperparameter matrix (aka, 'loadings').

\bold{\eqn{S}} is the Slot matrix.

\bold{\eqn{\eta}}\eqn{_i} is the row vector of latent states (aka, 'constructs', or 'scales') for observation \eqn{i}.

\bold{\eqn{\epsilon}}\eqn{_i} ~ \eqn{N[0, \Sigma_s]} is a column vector of system shocks for observation \eqn{i}.


}
\author{
Dave Zes, Jimmy Lewis, Dana Landis @ Korn/Ferry International

<zesdave@gmail.com>
}
\references{
Brown, A., & Maydeu-Olivares, A. (2012, November 12). How IRT Can Solve Problems of Ipsative Data in Forced-Choice Questionnaires. Psychological Methods. Advance online publication. doi: 10.1037/a0030641
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ IRT }
\keyword{ Item Response Theory }
\keyword{ Thurstonian }
\keyword{ psychometric }
\keyword{ forced choice  }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%~~ simple examples of the most important functions ~~
%}
