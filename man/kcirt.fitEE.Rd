\name{kcirt.fitEE}
\alias{kcirt.fitEE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Least Squares k-Cube Thurstonian IRT Fitting
}
\description{
k-Cube Thurstonian IRT Fitting using a least-squares expectation-expectation algorithm.
}
\usage{
kcirt.fitEE(model, mxHatLambda, maxIter = 40)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
	A kcirt model.  A named list of \code{class} 'kcube.irt.model'.
}
  \item{mxHatLambda}{
An initial guess for the Hyperparameters.
}
  \item{maxIter}{
Maximum number of iterations.
}
}
\details{
This function can be thought of as an expectation-expectation procedure.  The starting Hyperparameters, \code{mxHatLambda}, are used to predict \code{mxEta} (this prediction is commonly called \code{mxHatEta} is this package), and so on, back and forth.  The procedure stops when either the L2 cost first bottoms out, or \code{maxIter} is met.

In many cases, this function alone produces excellent-performing estimates/predictions.  The user may pass the returned model to \code{\link{kcirt.fitMSS}} for further refinement.
}
\value{
	A kcirt model.  A named list of \code{class} 'kcube.irt.model'.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Dave Zes, Korn/Ferry International
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also \code{\link{kcirt.fitMSS}}.
}
\examples{
constructMap.ls <- list(
c(1,1,2,2),
c(1,1,3,3),
c(2,2,3,3),
c(1,1,2,2),
c(1,1,3,3),
c(2,2,3,3)
)

qTypes <- rep("R", length(constructMap.ls))

mod <- kcirt.model(constructMap.ls=constructMap.ls, qTypes=qTypes, mxLambda=NULL)

N <- 300
set.seed(99999)
mod <- kcirt.sim(model=mod, N=N)

ikcirt.df1(mod, "self")

mxHatLambda <- mod$mxLambda - matrix( rnorm( sum(mod$ns)^2, 0, 0.3 ), sum(mod$ns), sum(mod$ns) )

mod2 <- kcirt.fitEE(model=mod, mxHatLambda=mxHatLambda, maxIter=40)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
