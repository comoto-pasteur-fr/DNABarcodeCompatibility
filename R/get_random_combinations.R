#' @title 
#' Get a large set of compatible combinations.
#'
#' @description 
#' Finds a randomly generated set of at most 1000 combinations of compatible barcodes.
#'
#' @usage 
#' get_random_combinations(index_df, mplex_level, chemistry)
#'
#' @param index_df A dataframe containing barcodes identifiers, corresponding DNA sequences along with GC content and presence of homopolymers.
#' @param mplex_level The number at which the barcodes will be multiplexed.
#' @param chemistry An integer representing the number of channels (1, 2, 4) of the Illumina plateform to be used.
#'
#' @details 
#' This function is suited if the total number of possible combinations is too high for an exhaustive search to be possible in a reasonable amount of time.
#'
#' @return 
#' A matrix containing the identifiers of compatible barcode combinations.
#'
#' @examples
#' write.table(DNABarcodeCompatibility::IlluminaIndexes,
#'  txtfile <- tempfile(), row.names = FALSE, col.names = FALSE, quote=FALSE)
#' barcodes <- file_loading_and_checking(txtfile)
#' get_random_combinations(barcodes, 3, 4)
#' 
#'
#' @seealso 
#' \code{\link{get_all_combinations}}
#' 
#' @export
#' 

get_random_combinations = function(index_df, mplex_level, chemistry){
  if (chemistry == 4 ){
    combinations_m = get_random_combinations_4_channel(index_df, mplex_level)
  } else if (chemistry == 2){
    combinations_m = get_random_combinations_2_channel(index_df, mplex_level)
  } else if (chemistry == 1){
    combinations_m = get_random_combinations_1_channel(index_df, mplex_level)
  } else {
    display_message("Please choose a correct chemistry for your experiment ")
  }
  return (combinations_m)
}
