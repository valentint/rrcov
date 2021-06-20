\name{CovMrcd}
\alias{CovMrcd}
%
\title{ Robust Location and Scatter Estimation via Minimum Regularized Covariance Determonant (MRCD)}
\description{
     Computes a robust multivariate location and scatter estimate with a high
     breakdown point, using the Minimum Regularized Covariance Determonant (MRCD) estimator.
}
\usage{
CovMrcd(x,
       alpha=control@alpha, 
       h=control@h,
       maxcsteps=control@maxcsteps,
       initHsets=NULL, save.hsets=FALSE,
       rho=control@rho,
       target=control@target,
       maxcond=control@maxcond,
       trace=control@trace,
       control=CovControlMrcd())
}
\arguments{
  \item{x}{a matrix or data frame. }
  \item{alpha}{numeric parameter controlling the size of the subsets
    over which the determinant is minimized, i.e., \code{alpha*n}
    observations are used for computing the determinant.  Allowed values
    are between 0.5 and 1 and the default is 0.5.}
  \item{h}{the size of the subset (can be between ceiling(n/2) and n). 
    Normally NULL and then it \code{h} will be calculated as 
    \code{h=ceiling(alpha*n)}. If \code{h} is provided, \code{alpha} 
    will be calculated as \code{alpha=h/n}.}
  \item{maxcsteps}{maximal number of concentration steps in the
    deterministic MCD; should not be reached.}
  \item{initHsets}{NULL or a \eqn{K x h} integer matrix of initial
    subsets of observations of size \eqn{h} (specified by the indices in
    \code{1:n}).}
  \item{save.hsets}{(for deterministic MCD) logical indicating if the
    initial subsets should be returned as \code{initHsets}.}
  \item{rho}{regularization parameter. Normally NULL and will be estimated from the data.}
  \item{target}{structure of the robust positive definite target matrix: 
  a) "identity": target matrix is diagonal matrix with robustly estimated 
  univariate scales on the diagonal or b) "equicorrelation": non-diagonal 
  target matrix that incorporates an equicorrelation structure 
  (see (17) in paper). Default is \code{target="identity"}}
  \item{maxcond}{maximum condition number allowed 
  (see step 3.4 in algorithm 1). Default is \code{maxcond=50}}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
  \item{control}{ a control object (S4) of class \code{\link{CovControlMrcd-class}}
    containing estimation options - same as these provided in the function
    specification. If the control object is supplied, the parameters from it
    will be used. If parameters are passed also in the invocation statement, they will
    override the corresponding elements of the control object.}
}
\details{
  This function computes the minimum regularized covariance determinant estimator (MRCD)
  of location and scatter and returns an S4 object of class
  \code{\link{CovMrcd-class}} containing the estimates.
  Similarly like the MCD method, MRCD looks for the \eqn{h (> n/2)}
  observations (out of \eqn{n}) whose classical
  covariance matrix has the lowest possible determinant, but 
  replaces the subset-based covariance by a regularized 
  covariance estimate, defined as a weighted average of the 
  sample covariance of the h-subset and a predetermined 
  positive definite target matrix. The Minimum Regularized Covariance 
  Determinant (MRCD) estimator is then the regularized covariance 
  based on the h-subset which makes the overall determinant the smallest.
  A data-driven procedure sets the weight of the target matrix (\code{rho}), so that 
  the regularization is only used when needed.
}
\value{
  An S4 object of class \code{\link{CovMrcd-class}} which is a subclass of the
  virtual class \code{\link{CovRobust-class}}.
}
\references{
  Kris Boudt, Peter Rousseeuw, Steven Vanduffel and Tim Verdonck (2018)
  The Minimum Regularized Covariance Determinant estimator.
  submitted, available at \url{https://arxiv.org/abs/1701.07086}.

  Mia Hubert, Peter Rousseeuw and Tim Verdonck (2012) A deterministic algorithm
  for robust location and scatter.
  \emph{Journal of Computational and Graphical Statistics} \bold{21}(3), 618--637.

  Todorov V & Filzmoser P (2009),
  An Object Oriented Framework for Robust Multivariate Analysis.
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  URL \url{http://www.jstatsoft.org/v32/i03/}.
}
\author{ 
Kris Boudt, Peter Rousseeuw, Steven Vanduffel and Tim Verdonk. Improved by Joachim Schreurs and Iwein Vranckx.
Adapted for \pkg{rrcov} by Valentin Todorov \email{valentin.todorov@chello.at}
}
\seealso{
  \code{\link{CovMcd}}
}
\examples{
## The result will be (almost) identical to the raw MCD
##  (since we do not do reweighting of MRCD)
##
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
c0 <- CovMcd(hbk.x, alpha=0.75, use.correction=FALSE)
cc <- CovMrcd(hbk.x, alpha=0.75)
cc$rho
all.equal(c0$best, cc$best)
all.equal(c0$raw.center, cc$center)
all.equal(c0$raw.cov/c0$raw.cnp2[1], cc$cov/cc$cnp2)

summary(cc)

## the following three statements are equivalent
c1 <- CovMrcd(hbk.x, alpha = 0.75)
c2 <- CovMrcd(hbk.x, control = CovControlMrcd(alpha = 0.75))
## direct specification overrides control one:
c3 <- CovMrcd(hbk.x, alpha = 0.75,
             control = CovControlMrcd(alpha=0.95))
c1

\dontrun{
data(octane)
n <- nrow(octane)
p <- ncol(octane)
out <- CovMrcd(octane, h=33) 
robpca = PcaHubert(octane, k=2, alpha=0.75, mcd=FALSE)
(outl.robpca = which(robpca@flag==FALSE))

# Observations flagged as outliers by ROBPCA:
# 25, 26, 36, 37, 38, 39

# Plot the orthogonal distances versus the score distances:
pch = rep(20,n); pch[robpca@flag==FALSE] = 17
col = rep('black',n); col[robpca@flag==FALSE] = 'red'
plot(robpca, pch=pch, col=col, id.n.sd=6, id.n.od=6)

## Plot now the MRCD mahalanobis distances
pch = rep(20,n); pch[!getFlag(out)] = 17
col = rep('black',n); col[!getFlag(out)] = 'red'
plot(out, pch=pch, col=col, id.n=6)
}
}
\keyword{robust}
\keyword{multivariate}
