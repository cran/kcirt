\name{kcirt.ystarinfo}
\alias{kcirt.ystarinfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate State Information
}
\description{
Calculate System Information imparted to States through \bold{\eqn{Y}}*.
}
\usage{
kcirt.ystarinfo(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
		A kcirt model.  A named list of \code{class} 'kcube.irt.model'.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A numeric, square, symmetric matrix, whose number of rows equals the number of model states (scales).
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
c(1,1,2,2),
c(1,1,3,3),
c(2,2,3,3),
c(1,1,2,2),
c(1,1,3,3),
c(2,2,3,3),
c(1,2),
c(2,3),
c(1,3)
)
    
    
qTypes <- rep("R", length(constructMap.ls))
    
mod <- kcirt.model(constructMap.ls=constructMap.ls, qTypes=qTypes, mxLambda=NULL)

kcirt.ystarinfo(mod)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
