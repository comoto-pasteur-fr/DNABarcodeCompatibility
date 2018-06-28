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
#' N/k compatible combinations are then selected using a Shannon entropy maximization approach.
#' It can be shown that the maximum value of the entropy that can be attained for a selection of N barcodes among n, with possible repetitions, reads: 
#'  \deqn{S_{max}=-(n-r)\frac{\lfloor N/n\rfloor}{N}\log(\frac{\lfloor N/n\rfloor}{N})-r\frac{\lceil N/n\rceil}{N}\log(\frac{\lceil N/n\rceil}{N})}
#'  
#'  where r denotes the rest of the division of N by n, while
#'  \deqn{\lfloor N/n\rfloor} and \deqn{\lceil N/n\rceil} denote
#'  the lower and upper integer parts of N/n, respectively.
#'
#'     + Case 1: number of lanes < number of compatible DNA-barcode combinations
#'     
#' This function seeks for compatible DNA-barcode combinations of highest entropy.
#' In brief this function uses a randomized greedy descent algorithm to find an optimized selection. 
#' Note that the resulting optimized selection may not be globally optimal.
#' It is actually close to optimal and much improved in terms of non-redundancy of DNA barcodes used, compared to a randomly choosen set of combinations of compatible barcodes.
#'
#'     + Case 2: number of lanes >= number of compatible DNA-barcode combinations
#'     
#' In such a case, there are not enough compatible DNA-barcode combinations and redundancy is inevitable.
#'
#' @md
#' 
#' @return 
#' A matrix containing an optimized set of combinations of compatible barcodes.
#'
#' @examples
#' m <- get_random_combinations(DNABarcodeCompatibility::IlluminaIndexes, 3, 4)
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

optimize_combinations = function (combination_m, nb_lane, index_number, thrs_size_comb=120, max_iteration=15, method="greedy_exchange"){
  # browser()
  if (nrow(as.matrix(combination_m)) == 0){
    display_message("No combinations have been found")
  } else {
    if(is.numeric(index_number)){
      if (is.numeric(nb_lane)){
        
        max = entropy_max(index_number, ncol(combination_m) * nb_lane)
        print(paste("Theoretical max entropy:",round(max, 3)))
        
        if(nb_lane < nrow(combination_m)){
          if (nrow(combination_m) > thrs_size_comb){ 
            a_combination = recursive_entropy(combination_m[sample(1:nrow(combination_m),thrs_size_comb),], nb_lane, method=method)
            i = 0
            while ((i < max_iteration) && (entropy_result(a_combination) < max) ){
              temp_combination = recursive_entropy(combination_m[sample(1:nrow(combination_m),thrs_size_comb),], nb_lane, method=method)
              if (entropy_result(temp_combination) > entropy_result(a_combination) ){
                a_combination = temp_combination
              }
              i = i+1
            }
          } else {
            if (nrow(combination_m) > 25) {
              a_combination = recursive_entropy(combination_m, nb_lane, method=method)
              i = 0
              while ( (i < max_iteration) && (entropy_result(a_combination) < max) ){
                temp_combination = recursive_entropy(combination_m, nb_lane, method=method)
                if (entropy_result(temp_combination) > entropy_result(a_combination) ){
                  a_combination = temp_combination
                }
                i = i+1
              }
            } else {
              a_combination = recursive_entropy(combination_m,nb_lane, method="greedy_descent")
              i = 0
              while ((i < max_iteration) && (entropy_result(a_combination) < max) ){
                temp_combination = recursive_entropy(combination_m, nb_lane, method="greedy_descent")
                if ( entropy_result(temp_combination) > entropy_result(a_combination) ){
                  a_combination = temp_combination
                }
                i = i+1
              }
            }
          }
        } else {
          n = nb_lane - nrow(combination_m)
          combination_m = combination_m[sample(1:nrow(combination_m),nrow(combination_m)),,drop = F]
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
        
        print(paste("Entropy of the optimized set:", round(entropy_result(a_combination),3)))
        return(a_combination)
        
      } else {
        display_message("Please enter a number as nb_lane")
      }
    } else {
      display_message("Please enter a number as index_number")
    }
  }
}

