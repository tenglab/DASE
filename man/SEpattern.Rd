% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SEpattern.R
\name{SEpattern}
\alias{SEpattern}
\title{Pattern of SE}
\usage{
SEpattern(se_fit_df, pattern_plot_df, cutoff_vector)
}
\arguments{
\item{se_fit_df}{fitted output file from SEfitspline}

\item{pattern_plot_df}{spline fitted line data from SEfitspline}

\item{cutoff_vector}{cutoff vector including lower and upper cutoff from SEfitspline
(default is c(-1.5,1.5))}
}
\value{
plots: list contain constitute enhancer fold change plot of each SE
se_segment_percent: SE segment percentage dataset
}
\description{
A function to get the pattern of SEs based on enhancer fitted line
}
\details{
This function will generate plot patterns of each SE and based on the intersection of
cutoff and fitted line, patterns are seperated into several segments. The percentage of
each segments are calculated by basemean of enhancer/sum of basemean within each SE
}
\examples{
seg_list <- SEpattern(SEfitspline_out,SEfitspline_cutoff)

}
