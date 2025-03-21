% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apps.R
\name{homeology.wrapper}
\alias{homeology.wrapper}
\title{homeology.wrapper}
\usage{
homeology.wrapper(
  junctions,
  width = 50,
  pad = 0,
  thresh = 0,
  stride = 0,
  genome,
  cores,
  flip = FALSE,
  bidirectional = TRUE,
  annotate = TRUE,
  savegMatrix = TRUE,
  outdir = "./"
)
}
\arguments{
\item{junctions}{Path to .vcf or bedpe or rds file of junctions or gGraph from which alt edges will be taken}

\item{width}{width around junction breakpoint around which to search for homeology}

\item{pad}{number of bases of padding around each sequence position (bin) to use when computing homeology, i.e. we then will be comparing 1 + 2*pad -mer sequences for edit distance}

\item{thresh}{string distance threshold for calling homeology in a bin}

\item{stride}{distance in bases between consecutive bins in which we will be measuring homeology}

\item{genome}{Path to .2bit or ffTrack .rds containing genome sequence}

\item{cores}{How many cores to use}

\item{flip}{if flip = FALSE, homeology search for -/- and +/+ junctions is done between a sequence and its reverse complement}

\item{bidirectional}{adding padding on both sides of each breakpoint (TRUE) or only in the direction of the fused side (FALSE)}

\item{annotate}{annotate edges in gGraph object and save it in working directory}

\item{savegMatrix}{save gMatrix object of edit distances}

\item{outdir}{output directory}
}
\description{
homeology.wrapper
}
