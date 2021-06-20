% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("C:/projects/statproj/R/rrcov/inst/NEWS.Rd"))
\name{NEWS}
\title{News for \R Package \pkg{rrcov}}
\encoding{UTF-8}

\section{CHANGES in rrcov VERSION 1.5-5 (2020-07-31)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item PcaHubert will crash if X is 1-dimensional and mcd=FALSE: fixed.
      \item Fixed 'noLD' issues in tlda.R and okg4.R
      \item Fixed a problem with wrong scores in PcaProj() (reported by Matthieu Lesnoff <matthieu.lesnoff@gmail.com>)
      \item Fixed a problem with nsamp="exact" or nsamp="best" in CovMve(), CovSest() (reported by Claudio Agostinelli) - 
        these functions, differently from CovMcd, should not take non-numeric 'nsamp'
    }
  }
}  
\section{CHANGES in rrcov VERSION 1.5-2 (2020-01-14)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item  Added parameter 'control' to Linda - to select the robust location and 
        covariance estimator to use in LDA. Now any estimator derived from class 'Cov'
        can be used, even such that are not in 'rrcov'. Return this parameter in the 
        returned S4 object.
      \item Linda returns now the location and covariance estimator used as 'covobj'. This
        is useful for controlling cross-validation, for example.
      \item Linda and LdaClassic use generalized inverse if the common covariance matrix is singular.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed an issue the 'predict' function.
      \item Removed the dependence on packages 'cluster' and 'ellipse'.
    }
  }
}
\section{CHANGES in rrcov VERSION 1.4-9 (2019-11-25)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item  Added data set \code{diabetes}; data set \code{soil} from package \code{rrcovHD} moved here.
      \item  Linear and quadratic discriminant analysis can use the MRCD estimates.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed an issue with CovControlMcd(nsamp="deterministic") - this would not work, 
        because nsamp was defined in the class definition as "numeric". Now it is "Cnumeric" - union of "character" and "numeric'.
      \item Corrected the code for Minimum Regularized Covariance Determinant estimator (MRCD) - 
        \code{CovMrcd()} - the step of adjusting the eignevalues in r6pack() is excluded 
        now because it has no effect when n > p.
    }
  }
}
\section{CHANGES in rrcov VERSION 1.4-7 (2018-11-14)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Added Minimum Regularized Covariance Determinant estimator (MRCD) 
        (Boudt et al. 2018) - 
        functions \code{CovMrcd()} and \code{CovControlMrcd()}
        
      \item  Added data set \code{octane}; data set \code{olitos} from package rrcovHD moved here.
      
      \item The 'pairs' plot is now available for classical covariance matrix
    }
  }
}
