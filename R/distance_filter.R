#' @title 
#' Select barcode combinations with error correction properties
#'
#' @description 
#' Filter barcode combinations for a given distance metric (hamming or seqlev) and threshold.
#'
#' @usage 
#' distance_filter(index_df, combinations_m, metric, d)
#'
#' @param index_df A dataframe containg barcode identifiers, corresponding DNA sequences along with GC content and presence of homopolymers.
#' @param combinations_m A matrix of compatible barcode combinations.
#' @param metric The type of distance (hamming or seqlev).
#' @param d The minimum value of the distance.
#'
#' @details 
#' The Hamming distance is suitable for correcting substitution errors. 
#' The SeqLev distance is suitable for correcting insertion/deletion errors.
#'
#' @return 
#' Returns a filtered matrix containing the identifiers of compatible combinations of barcodes.
#'
#' @examples
#' write.table(DNABarcodeCompatibility::IlluminaIndexes,
#'  txtfile <- tempfile(), row.names = FALSE, col.names = FALSE, quote=FALSE)
#' barcodes <- file_loading_and_checking(txtfile)
#' m <- get_all_combinations(barcodes, 3, 4)
#' distance_filter(barcodes, m, "hamming", 3)
#'
#' @seealso 
#' \code{\link{get_all_combinations}}, 
#' \code{\link{get_random_combinations}}
#' 
#' @references Buschmann, T. 2015. The Systematic Design and Application of Robust DNA Barcodes. 
#' @references Buschmann, T. 2017. DNABarcodes: an R package for the systematic construction of DNA sample tags. Bioinformatics 33, 920–922.
#' 
#' @export
#' 

# \link{http://www.qucosa.de/recherche/frontdoor/?tx_slubopus4frontend%5bid%5d=urn:nbn:de:bsz:15-qucosa-209812}

distance_filter = function(index_df, combinations_m, metric, d){
  index_distance_df =  index_distance(index_df)
  if (metric == "hamming"){
    hamming_rejection_table = low_hamming_distance(index_df, index_distance_df, d)
    filtered_combinations_m = filter_combinations(combinations_m, hamming_rejection_table)
  }else if (metric == "seqlev") {
    seqlev_rejection_table = low_seqlev_distance(index_df, index_distance_df, d) 
    filtered_combinations_m = filter_combinations(combinations_m, seqlev_rejection_table)
  }
  return (filtered_combinations_m)
}
