\name{kcirt.model}
\alias{kcirt.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create k-Cube Thurstonian IRT Model
}
\description{
Create a k-Cube Thurstonian IRT Model skeleton.
}
\usage{
kcirt.model(constructMap.ls, qTypes, data = NULL, Y = NULL, mu = 0, mxLambda = NULL, 
covEta = 1, covShocks = 1, deltaType=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{constructMap.ls}{
How the constructs map to the items.  A list of vectors, each vector representing a block (question in the instrument); each element of a vector giving the construct index.
}
  \item{qTypes}{
Is a question to be fully ranked or most/least format.  A character vector whose length is the number of blocks.  Each element in \code{\{'R', 'M'\}}.  Note that only accepts 'R' when corresponding block has 3 or 2 items.
}
  \item{data}{
Data matrix.  The number of columns is the number of observational units, the number of rows is the total number of items.
}
  \item{Y}{
Response cohort matrix, each element in \code{\{NA,0,1\}}.
}
  \item{mu}{
True item means.  A real-valued vector whose length is the total number of items.
}
  \item{mxLambda}{
True hyperparameters.  A real-valued square matrix, or a vector.  If vector, \code{mxLambda} is assumed to be diagonal (no item crosstalk) -- vector is recycled
}
  \item{covEta}{
True covariance of states.  A square, symmetric real-valued matrix.
}
  \item{covShocks}{
True covariance of shocks.  A square, symmetric real-valued matrix.
}
  \item{deltaType}{
Tell function in what pattern items are compared.  Scalar.  Directly affects structure of the Delta matrix.  Currently one of two integer values.  1 is default.  2 is that implied by Brown and Maydeau-Oliveras.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
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

## view contents
mod
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
