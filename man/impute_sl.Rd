% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impute_sl.R
\name{impute_sl}
\alias{impute_sl}
\title{Impute missing values in shotgun lipidomics data matrix}
\usage{
impute_sl(m)
}
\arguments{
\item{m}{a matrix of shotgun lipidomics data where rownames are lipid compositions, structures, or any other unique identifier, columns are samples or replicates, and all entries of matrix are coercible to numeric.}
}
\description{
This function performs a sample minimum imputation of a shotgun lipidomics data matrix where missing values are represented by NA
}
\examples{
m <- df \%>\% dplyr::select(c("rel.intensity","sample","comp")) \%>\% 
            tidyr::pivot_wider(names_from = sample,
                               values_from = rel.intensity,
                               values_fn = mean))

m.norm <- normalize_sl(m)

m.final <- impute_sl(m.norm)
}
\keyword{normalization}
