% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_to_fish_scoring.R
\name{convert_to_fish_scoring}
\alias{convert_to_fish_scoring}
\title{convert_to_fish_scoring}
\usage{
convert_to_fish_scoring(genelist)
}
\arguments{
\item{genelist}{list of human genes as character vector.}
}
\value{
a character vector of the fish orthologs of the given fish genes.
}
\description{
convert a human genelist (character vector) to fish, using DIOPT scores to remove duplicated genes. Intended to be used in cases where you don't have corresponding expression data for the genelist.
}
