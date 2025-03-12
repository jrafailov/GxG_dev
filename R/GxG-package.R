#' GxG
#'
#' @import methods
#' @import R6
#' @import data.table
#' @import Matrix
#' @import GenomicRanges
#' @import igraph
#' @importFrom reshape2 melt
#' @import gUtils
#' @import gTrack
#' @import strawr
#' @import data.table
#' @import gGnome
#' @import ggplot2 
#' @import parallel
#' @import httr
#' @import futile.logger
#' @importFrom rtracklayer TwoBitFile
#' @importFrom gUtils gr2dt dt2gr %Q% gr.val %$% %*% gr.chr gr.nochr gr.findoverlaps
#' @importFrom magrittr %>%
#' @importFrom parallel mclapply
## @importFrom khtools conform_si et copy2 printerr parse.gr2
#' @import khtools 
#' @importFrom GenomicRanges GRanges GRangesList
#' @importFrom GenomeInfoDb Seqinfo seqnames seqnames<- seqinfo seqinfo<- seqlengths seqlengths<- seqlevels seqlevels<- seqlevelsStyle seqlevelsStyle<-
#' @importFrom BiocGenerics width
#' @importMethodsFrom BiocGenerics width sort
#' @importMethodsFrom IRanges trim start end
#' @importMethodsFrom GenomeInfoDb seqnames seqnames<- seqinfo seqinfo<- seqlengths seqlengths<- seqlevels seqlevels<- seqlevelsStyle seqlevelsStyle<- sortSeqlevels
#' @importMethodsFrom GenomicRanges start end
#' @importMethodsFrom MatrixGenerics rowRanges
#' @importMethodsFrom S4Vectors split mcols mcols<- values values<-
#' @importFrom Rcpp sourceCpp
#' @useDynLib GxG
registerS3method(genname = "merge", class = "data.table", method = data.table::merge.data.table)
"_PACKAGE"
