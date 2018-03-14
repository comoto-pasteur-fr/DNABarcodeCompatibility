#' @title 
#' Get compatible combinations of DNA barcodes
#'
#' @description 
#' Performs either an exhaustive or a random search of compatible DNA-barcode combinations depending on the size of the DNA-barcode population and of the number of samples to be multiplexed.
#'
#' @usage 
#' get_combinations(index, multiplexing_level)
#'
#' @param 
#' index is a dataframe containing the barcode names, the DNA sequence and the corresponding binary sequence.
#' @param 
#' multiplexing_level is the number of samples to be multiplexed for sequencing.
#'
#' @return 
#' Returns compatible DNA-barcode combinations as a matrix.
#'
#' @examples
#' write.table(DNABarcodeCompatibility::illumina, file <- tempfile(), row.names = FALSE, col.names = FALSE, quote=FALSE)
#' barcodes <- file_loading_and_checking(file) 
#' index <- transform(barcodes, binary = sequence_binary_conversion(sequence))
#' comb <- get_combinations(index,multiplexing_level=3)
#' head(comb); tail(comb)
#'
#' @export
#' 
get_combinations = function (index, multiplexing_level){
  if (choose(nrow(index),multiplexing_level) <= 2024){
    return(get_all_combinations(index, multiplexing_level))
  }else {
    return(get_random_combinations(index, multiplexing_level))
  }
}


# get all compatible combinations of an index
#' @describeIn get_combinations force to get all compatible combinations of DNA barcodes
#' @examples
#' allComb <- get_all_combinations(index,multiplexing_level=3)
#' head(allComb); tail(allComb)
#'
#' @export
get_all_combinations = function(index, multiplexing_level){
  index = index %>% arrange(Id)
  matrix_id = combn(x = index$Id, m = multiplexing_level)
  matrix_binary_sequence = matrix_id_to_binary_sequence(matrix_id = matrix_id, index = index)#matches Ids to binary sequences
  logical_rigth_combination = as.logical(x = list_of_good_combinations(m = matrix_binary_sequence))
  list_of_all_combinations = matrix_id[, logical_rigth_combination] %>% t()
  return(list_of_all_combinations)
}


# For a random search
#' @describeIn get_combinations force to perform random search to get some compatible combinations of DNA barcodes.
#' @details 
#' The \code{\link{get_random_combinations}} function generates k-combinations of a set S composed of n elements by random sampling, and it filters out incompatible combinations on the fly. It stops when 700 k-combinations have been tested for compatibility.  HERE MOTIVATE nrow = 700.
#' @examples
#' ranComb <- get_random_combinations(index,multiplexing_level=3)
#' head(ranComb); tail(ranComb)
#'
#' @export

get_random_combinations = function (index, multiplexing_level){
  list_of_good_combs = matrix(nrow = 700, ncol = multiplexing_level)
  j = 0
  for (i in 1:700){
    combination =index[sample(x = 1:nrow(index), size = (multiplexing_level), replace = FALSE),]
    if(is_a_good_combination(combination$binary)){
      j = j+1
      combination = arrange(combination, Id) #facilitates the distance filter
      list_of_good_combs [j,] = combination$Id
    }
  }
  M = list_of_good_combs %>% unique() %>% na.omit()
  M = M[order(M[,1]),]#facilitates the distance filter
  return (M)
}






