\name{ikcirt.Y2data}
\alias{ikcirt.Y2data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert Y to Rank Data
}
\description{
Converts a dichotomous response matrix, Y, to rank data.  Probably no need for user to call; called by \code{\link{kcirt.model}}.
}
\usage{
ikcirt.Y2data(Y, mxDelta, ns)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
Response matrix.
}
  \item{mxDelta}{
      System Delta matrix.  A \eqn{q} x \eqn{m} matrix with elements in \code{\{-1,0,1\}}.
}
  \item{ns}{
Response dimensions. A vector of length equal to the number of questions, each element giving the number of items in the question.
}

}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Ranked Data matrix.
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
For the inverse of this function, see \code{\link{ikcirt.data2Y}}.
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

#### convert data to Y
Y <- ikcirt.data2Y(mxData=mxData, mxDelta=model$mxDelta)

#### convert Y back to data
ikcirt.Y2data(Y=Y, mxDelta=model$mxDelta, ns=model$ns)




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
