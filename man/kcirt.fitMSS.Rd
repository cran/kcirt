\name{kcirt.fitMSS}
\alias{kcirt.fitMSS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Metaheuristic k-Cube Thurstonian IRT Fitting
}
\description{
Use metaheuristic stochastic search to locate k-Cube Thurstonian IRT hyperparameters and states.
}
\usage{
kcirt.fitMSS(model, lambdaConstraint = "self", kcpus = 2, penalty = "logit", 
usetruesigma = TRUE, mss.sd = 0.2, nsearch = 19, l2zvarpow = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
A kcirt model; a named list of \code{class} "kcube.irt.model".
}
  \item{lambdaConstraint}{
Scalar string specifying how to constrain Lambda during fitting.  See Details.
}
  \item{kcpus}{
Scalar positive integer telling \code{snowfall} how many threads to initialize.  Set to 1 for no parallel processing.
}
  \item{penalty}{
Scalar string.  How to measure quality of fit?  Currently either 'logit', 'L2', 'L2c', 'miscat'.
}
  \item{usetruesigma}{
Scalar boolean.  Use actual (assumed) system variance?
}
  \item{mss.sd}{
Scalar positive number or vector of length 3.  Size of search function standard deviation for, in order, \code{mu}, \code{Lambda}, and \code{Eta} searches.
}
  \item{nsearch}{
Scalar positive integer.  How many candidate values to draw from under the MSS function.  Typically 14-20.
}
  \item{l2zvarpow}{
Scalar positive integer.  Only used when \code{penalty} is 'L2c'.  Power by which to raise \code{var(2*Y-1)}; the amount of 'smoothing' to impose on \code{2*Y-1}.  Note that when zero, setting \code{penalty} to 'L2c' produces identical results to 'L2' -- though slightly slower computationally.
}
}
\details{
\bold{lambdaConstraint} defines the model k-cube, i.e., the item crosstalk space.  Currently accepts \code{self}, only diagonal elements of \code{Lambda}; \code{withinx}, additionally includes within block items pointing to e\bold{x}ogenous constructs; \code{withini}, additionally includes items within block pointing to same construct.  \code{betweenx}, diagonal elements and all other items that point to different constructs; \code{betweeni}, additionally includes items that point to same construct.  \code{priorx} is similar to \code{betweenx} except that only previous items crosstalk with item, i.e., \code{Lambda} is lower-block-diagonal; \code{priori} additionally permits crosstalk between prior items pointing to same construct.

\bold{penalty} defines the ojective function to be minimized.  \code{logit} minimizes the logistic deviance with shrinkage ||eta||.  \code{L2} minimizes sum of squares between \code{2*Y-1} and \code{2*Yhat-1}.  \code{L2c} is similar to \code{L2}, except that it pre-projects (up to a constant) \code{2*Y-1} and \code{2*Yhat-1} using \code{var(2*Y-1)}.  Note that \code{NA} values in \code{Y} are replaced with 0.5.  Finally, \code{miscat}, minimizes the misclassification rate using a confusion table between the observations, \code{Y}, and the predicted states (scales), \code{hatEta}, split at zero.

Excepting very adventurous explorations with this model, the defaults for \bold{lambdaConstraint} and \bold{penalty} will likely best serve the user.
}
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

\seealso{
See Also \code{\link{kcirt.fitEE}}.
}
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
    
N <- 200
set.seed(99999)
mod <- kcirt.sim(model=mod, N=N)
    
    
ikcirt.df1(mod, "self")
    
####### create initial guess for hyperparameters (aka loadings)
mod$mxHatLambda <- mod$mxLambda - matrix( rnorm( sum(mod$ns)^2, 0, 0.3 ), sum(mod$ns), sum(mod$ns) )
    
####### need to assign hat states and utilities
mod$mxHatEta <- matrix(0, N, sum(mod$nuc))
mod$hatMu <- rep(0, sum(mod$ns))

\dontrun{
###### run MSS fit -- performance is R^2 btwn true and pred states (aka scales)
mod <- kcirt.fitMSS(model=mod, lambdaConstraint="self", kcpus=2, penalty="L2", 
usetruesigma=TRUE, mss.sd=1)
mod$performance
    
###### run again ...
mod <- kcirt.fitMSS(model=mod, lambdaConstraint="self", kcpus=2, penalty="L2", 
usetruesigma=TRUE, mss.sd=1)
mod$performance
    
###### run a few more times ... ...
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
