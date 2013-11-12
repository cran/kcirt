\name{kcirt.sim}
\alias{kcirt.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a k-Cube Thurstonian IRT Model
}
\description{
Given model parameters, create a random realization of a k-Cube Thurstonian IRT Model skeleton.
}
\usage{
kcirt.sim(model, N, type = "Eta")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
A kcirt model.  A named list of \code{class} 'kcube.irt.model'.
}
  \item{N}{
Number of observations.
}
  \item{type}{
Currently ignored.
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
## Please see example in, e.g., kcirt.fitMSS.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
