\name{ikcirt.fun.mss.mu}
\alias{ikcirt.fun.mss.mu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MSS Locate mu
}
\description{
Used by \code{\link{kcirt.fitMSS}}; not to be called by user.
}
\usage{
ikcirt.fun.mss.mu(jj, iimss, rndTrys, hatMu, useSysCov, penalty)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{jj}{
Index for \code{\link{snowfall}} cluster call.
}
  \item{iimss}{
Index of \code{hatMu} to presently search.
}
  \item{rndTrys}{
Values to try.
}
  \item{hatMu}{
Estimated item means.
}
  \item{useSysCov}{
Either the true or estimated System Variance matrix.
}
  \item{penalty}{
Scalar string.  How to measure quality of fit?  Currently either 'logit', 'L2', 'L2c', 'miscat'.
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
See Also \code{\link{kcirt.fitMSS}}, \code{\link{ikcirt.fun.mss.lambda}}, \code{\link{ikcirt.fun.mss.eta}}.
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
