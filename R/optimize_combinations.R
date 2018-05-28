#' @title 
#' Find a set of barcode combinations with least redundant barcodes
#'
#' @description 
#' This function uses the Shannon Entropy to identify a set of compatible barcode combinations with least redundancy between DNA barcodes.
#'
#' @usage 
#' optimize_combinations(combination_m, nb_lane, index_number)
#'
#' @param 
#' combination_m A matrix of compatible barcode combinations.
#' nb_lane The number of lanes to be use for sequencing (i.e. the number of libraries divided by the multiplex level).
#' index_number The total number of distinct DNA barcodes in the dataset.
#'
#' @details 
#' 
#'
#' @return 
#' Returns a matrix containing an optimized set of combinations of compatible barcodes.
#'
#' @examples
#' write.table(DNABarcodeCompatibility::illumina, file <- tempfile(), row.names = FALSE, col.names = FALSE, quote=FALSE)
#' barcodes <- file_loading_and_checking(file)
#' m <- get_random_combinations(barcodes, 3, 4)
#' optimize_combinations(m, 12, 48)
#' 
#'
#' @seealso 
#' \code{\link{get_all_combinations}}
#' \code{\link{get_random_combinations}}
#' 
#' @export
#' 

optimize_combinations = function (combination_m, nb_lane, index_number){
  max = entropy_max(index_number, length(combination_m) * nb_lane)
  sample_combination = combination_m[sample(1:nrow(combination_m),100)]
  if(nb_lane < nrow(combination_m)){
    if(combination_m > 100){
      possible_combination = 
        for (i in 1:10){
          sample_combination = combination_m[sample(1:nrow(combination_m),100)]
          one_combination = recursive_entropy()
        }
    }
    
  }
}
