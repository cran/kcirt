\name{ikcirt.data2Y}
\alias{ikcirt.data2Y}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert Rank Data to Y
}
\description{
Converts raw rank data to the dichotomous response matrix, Y.  Probably no need for user to call; called by \code{\link{kcirt.model}} and \code{\link{kcirt.sim}}.
}
\usage{
ikcirt.data2Y(mxData, mxDelta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mxData}{
Ranked Data.  A matrix of positive integers (rankings).
}
  \item{mxDelta}{
      System Delta matrix, elements in \code{\{-1,0,1\}}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
    Response matrix with elements in \code{\{NA,0,1\}}.
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
For the inverse of this function, see \code{\link{ikcirt.Y2data}}.
}
\examples{

##### here's an itty-bitty example
	
constructMap.ls <- list(
c(1,1,2,2),
c(2,2,3,3),
c(1,1,3,3)
)

qTypes <- rep("M", length(constructMap.ls))

model <- kcirt.model(constructMap.ls=constructMap.ls, qTypes=qTypes, mxLambda=NULL)

#### fake data, two cases (cases run column-wise)
mxData <- t(
rbind(
c( c(1, NA, NA, 4), c(4, NA, NA, 1), c(NA, 4, NA, 1) ),
c( c(1, 4, NA, NA), c(NA, NA, 1, 4), c(NA, 4, 1, NA) )
)
)

ikcirt.data2Y(mxData, mxDelta=model$mxDelta)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
