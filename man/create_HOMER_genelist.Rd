% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_HOMER_genelist.R
\name{create_HOMER_genelist}
\alias{create_HOMER_genelist}
\title{create_HOMER_genelist}
\usage{
create_HOMER_genelist(genelist, save_as, remove_ribosomal_genes = T)
}
\arguments{
\item{genelist}{List of genes you want to perform HOMER on. Must include "gene" and "avg_logFC" columns.}

\item{save_as}{name that HOMER txt file will be saved as, don't have to include txt at end.}

\item{remove_ribosomal_genes}{default = T, remove ribosomal genes from FindMarkers results?}
}
\value{
doesn't return anything into R, just saves txt file in the working directory.
}
\description{
Create txt file of genes with ENSEMBL IDs to be used for HOMER.
}
