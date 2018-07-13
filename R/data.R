#' Barcode dataset from Illumina.
#'
#' 48 barcodes from Illumina TruSeq Small RNA kits
#'
#' @format A data frame with 48 rows and 2 variables:
#' \describe{
#' \item{V1}{barcode identifier}
#' \item{V2}{DNA sequence}
#' ...
#' }
#'
"IlluminaIndexesRaw"

#' Barcode dataset from Illumina with features.
#'
#' 48 barcodes from Illumina TruSeq Small RNA kits along 
#' with percentage in CG content and presence of homopolymers of size >= 3
#'
#' @format A data frame with 48 rows and 4 variables:
#' \describe{
#' \item{Id}{barcode identifier}
#' \item{sequence}{DNA sequence}
#' \item{GC_content}{percentage of G and C relative to A and T}
#' \item{homopolymer}{presence of nucleotide repetitions (>= 3)}
#' ...
#' }
#'
"IlluminaIndexes"

