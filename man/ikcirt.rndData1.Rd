\name{ikcirt.rndData1}
\alias{ikcirt.rndData1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random Data
}
\description{
Generate a random realization of a minimally defined k-Cube Thurstonian IRT Model.  Probably no need for user to call; called by \code{\link{kcirt.sim}}.
}
\usage{
ikcirt.rndData1(N, qTypes, mxDelta, ns)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
Number of observational units.
}
  \item{qTypes}{
Is a question to be fully ranked or most/least format. A character vector whose length is the number of blocks. Each element in {'R', 'M'}.
}
  \item{mxDelta}{
Delta matrix.
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
constructMap.ls <- list(
c(1,2),
c(2,3),
c(1,3)
)

qTypes <- rep("R", length(constructMap.ls))

mod <- kcirt.model(constructMap.ls=constructMap.ls, qTypes=qTypes, mxLambda=NULL)

N <- 50
set.seed(99999)
mod <- kcirt.sim(model=mod, N=N)

ikcirt.rndData1(N=N, qTypes=qTypes, mxDelta=mod$mxDelta, ns=mod$ns)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
