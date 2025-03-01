% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MetaAnalysis.R
\name{performMetaAnalysis}
\alias{performMetaAnalysis}
\title{Title}
\usage{
performMetaAnalysis(resultsFolder, maxCores = 8)
}
\arguments{
\item{resultsFolder}{Folder containing the per-database results, and where
The meta-analytic estimates will be written. Each DB
is assumed to have its own sub-folder.}

\item{maxCores}{Maximum number of cores used for parallel processing.}
}
\value{
Does not return anything. Produces CSV files in the results folder.
}
\description{
Title
}
