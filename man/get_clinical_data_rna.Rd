% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_clinical_data.R
\name{get_clinical_data_rna}
\alias{get_clinical_data_rna}
\title{Get clinical data for RNA samples}
\usage{
get_clinical_data_rna(filter = NULL, wes = FALSE)
}
\arguments{
\item{filter}{A character indicating the filter to apply to the data. Options are:
NULL, "normals", "purity", "tissues", "breast"}

\item{wes}{A logical indicating if clinical data for samples with RNA and WES should be returned}
}
\value{
A data frame with sample level clinical information
}
\description{
Get clinical data for RNA samples
}
\examples{
get_clinical_data_rna()
get_clinical_data_rna("normals")
}
