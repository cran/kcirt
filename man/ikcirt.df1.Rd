\name{ikcirt.df1}
\alias{ikcirt.df1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conservative Degrees of Freedom
}
\description{
A simple, conservative calculation of a model's degrees of freedom.  Each question-respondent counts as 1 \eqn{df}.
}
\usage{
ikcirt.df1(model, lambdaConstraint = "self")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
	A kcirt model.  A named list of \code{class} 'kcube.irt.model'.  Must possess an object from which number of respondents can be gleaned, e.g., \code{Y}.
}
  \item{lambdaConstraint}{
A scalar string.  Please see \code{\link{kcirt.fitMSS}}.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Scalar integer.
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
## Please see example in kcirt.fitMSS.


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
