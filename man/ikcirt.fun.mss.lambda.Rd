\name{ikcirt.fun.mss.lambda}
\alias{ikcirt.fun.mss.lambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSS Locate Lambda
}
\description{
Used by \code{\link{kcirt.fitMSS}}; not to be called by user.
}
\usage{
ikcirt.fun.mss.lambda(jj, iimss, jjmss, rndTrys, mxHatLambda, penalty, usetruesigma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{jj}{
Index for \code{\link{snowfall}} cluster call.
}
  \item{iimss}{
Row index of \code{mxHatLambda} to presently search.
}
  \item{jjmss}{
Column index of \code{mxHatLambda} to presently search.
}
  \item{rndTrys}{
Values to try.
}
  \item{mxHatLambda}{
Estimated Hyperparameter matrix.
}
  \item{penalty}{
Scalar string.  How to measure quality of fit?  Currently either 'logit', 'L2', 'L2c', 'miscat'.
}
  \item{usetruesigma}{
Scalar boolean.  Use actual (assumed) system variance?
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Fit cost.  A scalar real-valued number.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Dave Zes, Korn/Ferry International
}
\note{
I cannot imagine a situation where the user could profit by calling this function.  It is not namespaced for its possible didactic value.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also \code{\link{kcirt.fitMSS}}, \code{\link{ikcirt.fun.mss.mu}}, \code{\link{ikcirt.fun.mss.eta}}.
}

%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
