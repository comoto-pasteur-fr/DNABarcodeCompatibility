#' @title 
#' Find a set of barcode combinations with least redundant barcodes
#'
#' @description 
#' This function uses the Shannon Entropy to identify a set of compatible barcode combinations with least redundancy between DNA barcodes.
#'
#' @usage 
#' optimize_combinations(index_df, sample_number, multiplexing_level, chemistry, metric = NULL, d = 3)
#'
#' @param 
#' index_df A dataframe containg barcodes identifiers, corresponding DNA sequences along with GC content and presence of homopolymers.
#' sample_number The number of libraries to be sequenced.
#' multiplexing_level The number at which the barcodes will be multiplexed.
#' chemistry An integer representing the number of channels (1, 2, 4) of the desired Illumina plateform.
#' metric The type of distance (hamming or seqlev).
#' d The minimum value of the distance.
#'
#' @details 
#' 
#'
#' @return 
#' Returns a matrix containing an optimized set of compatible combinations of barcodes.
#'
#' @examples
#' write.table(DNABarcodeCompatibility::illumina, file <- tempfile(), row.names = FALSE, col.names = FALSE, quote=FALSE)
#' barcodes <- file_loading_and_checking(file)
#' optimize_combinations(barcodes, 9, 3, 4, "hamming", 3)
#' 
#'
#' @seealso 
#' \code{\link{get_all_combinations}}
#' \code{\link{get_random_combinations}}
#' 
#' @export
#' 

optimize_combinations = function (index_df, sample_number, multiplexing_level, chemistry, metric = NULL, d = 3){
  combinations_m = get_combinations(index_df, multiplexing_level, chemistry)
  if(!is.null(metric)){
    combinations_m = distance_constraints_filter (index_df, combinations_m, metric, d)
  }
  nb_lane = sample_number %>% as.numeric() / multiplexing_level %>% as.numeric()
  cb = comb_result(combinations_m, sample_number, index_number, nb_lane) %>% as.data.frame()
  result = data.frame(Id = as.vector(cb%>% t() %>% as.vector),
                      Lane = (rep(1:nb_lane, length.out = sample_number, each = multiplexing_level)))
  result$Id = as.character(result$Id)
  return(result)
}
