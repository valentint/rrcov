######
##  VT::03.08.2019
##
##
##  roxygen2::roxygenise("C:/projects/statproj/R/rrcov", clean=TRUE)
##
#'
#'
#'
#' Computer Hardware
#'
#' A data set containing relative CPU performance data of 209 machines on 8 variables.
#;  The \code{rownames} are the vendor and model descriptions. Six of the variables
#'  are predictive, one (\code{PRP}) is the goal field and one (\code{ERP}) is the
#'  linear regression's guess. The estimated relative performance values were
#'  estimated by the authors using a linear regression method.  See their article
#'  (Ein-Dor and Feldmesser, CACM 4/87, pp 308-317) for more details on how the
#'  relative performance values were set.
#'
#' @name machines
#' @docType data
#' @usage data(machines)
#' @format A data frame with 209 rows and 8 variables
#' The variables are as follows:
#'
#' \itemize{
#'   \item MMIN: minimum main memory in kilobytes (integer)
#'   \item MMAX: maximum main memory in kilobytes (integer)
#'   \item CACH: cache memory in kilobytes (integer)
#'   \item CHMIN: minimum channels in units (integer)
#'   \item CHMAX: maximum channels in units (integer)
#'   \item PRP: published relative performance (integer)
#'   \item ERP: estimated relative performance from the original article (integer)
#' }
#'
#' @source \href{http://archive.ics.uci.edu/ml/datasets/Computer+Hardware?ref=datanews.io}{UCI Archive}
#'
#' @references
#'  Phillip Ein-Dor and Jacob Feldmesser (1987), Attributes of the performance
#'      of central processing units: A relative performance prediction model,
#'      \emph{Communications of the ACM}, \bold{30}, 4, pp 308-317.
#'
#'  Kibler, D., Aha, D.W. and Albert, M. (1989).  Instance-based prediction
#'      of real-valued attributes.  \emph{Computational Intelligence}, Vo.l 5,
#'      51-57.
#'
#' @examples
#'
#'  data(machines)
#'
#'  ## Compute the medcouple of each variable of the Computer hardware data
#'      data.frame(MC=round(apply(machines, 2, mc),2))
#'
#'  ## Plot a pairwise scaterplot matrix
#'      pairs(machines[,1:6])
#'
#'      mcd <- CovMcd(machines[,1:6])
#'      plot(mcd, which="pairs")
#'
#'  ##  Remove the rownames (too long)
#'      rownames(machines) <- NULL
#'
#'  ## Start with robust PCA based on MCD (P << n)
#'      (pca1 <- PcaHubert(machines, k=3))
#'      plot(pca1, main="ROBPCA-MCD", off=0.03)
#'
#'  ## PCA with the projection algoritm of Hubert
#'      (pca2 <- PcaHubert(machines, k=3, mcd=FALSE))
#'      plot(pca2, main="ROBPCA-SD", off=0.03)
#'
#'  ## PCA with the adjusted for skewness algorithm of Hubert et al (2009)
#'      (pca3 <- PcaHubert(machines, k=3, mcd=FALSE, skew=TRUE))
#'      plot(pca3, main="ROBPCA-AO", off=0.03)
#'
#' @keywords datasets
NULL

#'
#' Skull dimensions of the wolf \emph{Canis lupus} L.
#'
#' A data set containing skull morphometric measurements on Rocky Mountain
#'  and Arctic wolves (\emph{Canis Lupus L.}). The tdata are published in Morrison (1990),
#'  originally from Jolicoeur (1959).
#'
#' @name wolves
#' @docType data
#' @usage data(wolves)
#' @format A data frame with 25 rows and 12 variables.
#' The variables are as follows (all measurements are in milimeters):
#'
#' \itemize{
#'      \item \code{class}: a factor presenting the combinations of \code{location}
#'          and \code{sex}. The levels are \code{arf} \code{arm} \code{rmf} and \code{rmm}
#'      \item \code{location}: a factor with levels \code{ar}=Arctic, \code{rm}=Rocky Mountain
#'      \item \code{sex}: a factor with levels \code{f}=female, \code{m}=male
#'      \item \code{x1}: palatal length
#'      \item \code{x2}: postpalatal length
#'      \item \code{x3}: zygomatic width
#'      \item \code{x4}: palatal width outside first upper molars
#'      \item \code{x5}: palatal width inside second upper molars
#'      \item \code{x6}: postglenoid foramina width
#'      \item \code{x7}: interorbital width
#'      \item \code{x8}: braincase width
#'      \item \code{x9}: crown length
#' }
#'
#' @source
#'  Jolicoeur, P. Multivariate geographical variation in the wolf \emph{Canis lupis L.},
#'  \emph{Evolution}, XIII, 283--299.
#'
#'  Morrison, D. F.  \emph{Multivariate Statistical Methods},  (3rd ed.), 1990.
#'  New York: McGraw-Hill, p. 288--289.
#'
#' @examples
#'
#'  data(wolves)
#'
#'  ## Remove the factors location and sex which we will not use for now
#'  x <- wolves[,-c(2:3)]
#'
#'  ## Plot a pairwise scaterplot matrix
#'  pairs(x[,2:10])
#'
#'  mcd <- CovMcd(x[, 2:10])
#'  plot(mcd, which="pairs")
#'
#'  lda <- LdaClassic(class~., data=x)
#'  lda@center
#'  lda@cov
#'
#'  predict(lda)
#'
#' @keywords datasets
NULL
