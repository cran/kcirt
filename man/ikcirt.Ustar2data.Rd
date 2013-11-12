\name{ikcirt.Ustar2data}
\alias{ikcirt.Ustar2data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert \code{U}* to Data
}
\description{
Convert the non-differenced latent responses to Data.  Probably no need for user to call; called by \code{\link{kcirt.sim}}.
}
\usage{
ikcirt.Ustar2data(Ustar, qTypes, mxDelta, ns)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ustar}{
Non-differenced latent responses matrix.
}
  \item{qTypes}{
Is a question to be fully ranked or most/least format. A character vector whose length is the number of blocks. Each element in \code{\{'R', 'M'\}}.
}
  \item{mxDelta}{
System Delta matrix, elements in \code{\{-1,0,1\}}.
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

%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
