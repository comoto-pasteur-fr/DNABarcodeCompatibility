#' @title 
#' Find a set of barcode combinations with least redundant barcodes
#'
#' @description 
#' This function uses the Shannon Entropy to identify a set of compatible barcode combinations with least redundancy between DNA barcodes.
#'
#' @usage 
#' optimize_combinations(combination_m, nb_lane, index_number)
#'
#' @param combination_m A matrix of compatible barcode combinations.
#' @param nb_lane The number of lanes to be use for sequencing (i.e. the number of libraries divided by the multiplex level).
#' @param index_number The total number of distinct DNA barcodes in the dataset.
#'
#' @details 
#' 
#'
#' @return 
#' Returns a matrix containing an optimized set of combinations of compatible barcodes.
#'
#' @examples
#' write.table(DNABarcodeCompatibility::illumina,
#'  txtfile <- tempfile(), row.names = FALSE, col.names = FALSE, quote=FALSE)
#' barcodes <- file_loading_and_checking(txtfile)
#' m <- get_random_combinations(barcodes, 3, 4)
#' optimize_combinations(m, 12, 48)
#' 
#'
#' @seealso 
#' \code{\link{get_all_combinations}}, 
#' \code{\link{get_random_combinations}}, 
#' \code{\link{experiment_design}}
#' 
#' @export
#' 

optimize_combinations = function (combination_m, nb_lane, index_number){
  max = entropy_max(index_number, ncol(combination_m) * nb_lane)
  print(max)
  if(nb_lane < nrow(combination_m)){
    if (nrow(combination_m > 80)){ 
      a_combination = recursive_entropy(combination_m[sample(1:nrow(combination_m),80),],nb_lane)
      i = 0
      while ((i < 10) && (entropy_result(a_combination) < max) ){
        temp_combination = recursive_entropy(combination_m[sample(1:nrow(combination_m),80),], nb_lane)
        if (entropy_result(temp_combination) > entropy_result(a_combination) ){
          a_combination = temp_combination
        }
        i = i+1
      }
    } else {
      a_combination = recursive_entropy(combination_m,nb_lane)
      i = 0
      while ((i < 10) && (entropy_result(a_combination) < max) ){
        temp_combination = recursive_entropy(combination_m, nb_lane)
        if (entropy_result(temp_combination) > entropy_result(a_combination) ){
          a_combination = temp_combination
        }
      }
    }
  } else {
    n = nb_lane - nrow(combination_m)
    combination_m = combination_m[sample(1:nrow(combination_m),nrow(combination_m)),]
    part_combination = combination_m[sample(1:nrow(combination_m),n, replace = TRUE),]
    a_combination = rbind(combination_m, part_combination)
    i = 0
    while ((i < 30) && (entropy_result(a_combination) < max) && n > 0){
      part_combination = combination_m[sample(1:nrow(combination_m),n, replace = TRUE),]
      temp_combination = rbind(combination_m, part_combination)
      if (entropy_result(temp_combination) > entropy_result(a_combination) ){
        a_combination = temp_combination
      }
      i = i + 1
    }
    a_combination = a_combination[sample(1:nrow(a_combination),nrow(a_combination)),]
  }
  print(entropy_result(a_combination))
  return(a_combination)
}

