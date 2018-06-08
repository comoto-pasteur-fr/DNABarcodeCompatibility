#' @title 
#' Get all compatible combinations.
#'
#' @description 
#' Finds the exhaustive set of compatible barcode combinations.
#'
#' @usage 
#' get_all_combinations(index_df, mplex_level, chemistry)
#'
#' @param index_df A dataframe containing barcodes identifiers, corresponding DNA sequences along with GC content and presence of homopolymers.
#' @param mplex_level The number at which the barcodes will be multiplexed.
#' @param chemistry An integer representing the number of channels (1, 2, 4) of the Illumina plateform to be used.
#'
#' @details 
#' Be aware that the total number of combinations may become prohibitively large for large barcode sets and large multiplexing numbers.
#'
#' @return 
#' A matrix containing the identifiers of compatible barcode combinations.
#'
#' @examples
#' write.table(DNABarcodeCompatibility::IlluminaIndexes,
#'  txtfile <- tempfile(), row.names = FALSE, col.names = FALSE, quote=FALSE)
#' barcodes <- file_loading_and_checking(txtfile)
#' get_all_combinations(barcodes, 3, 4)
#' 
#'
#' @seealso 
#' \code{\link{get_random_combinations}}
#' 
#' @export
#' 

get_all_combinations = function(index_df, mplex_level, chemistry){
  if (is.numeric(mplex_level)){
    if (chemistry == 4 ){
      combinations_m = get_all_combinations_4_channel(index_df, mplex_level)
    } else if (chemistry == 2){
      combinations_m = get_all_combinations_2_channel(index_df, mplex_level)
    } else if (chemistry == 1){
      combinations_m = get_all_combinations_1_channel(index_df, mplex_level)
    } else {
      display_message("Please choose a correct chemistry for your experiment ")
    }
    if(nrow(as.matrix(combinations_m)) == 0){
      display_message("The programm did'nt find any compatible combination, please check your index list")
    }else{
      return (combinations_m)
    }
  }else{
    display_message("please enter a number as mplex_level")
  }
  
}
