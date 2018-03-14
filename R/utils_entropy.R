#' Group of utility functions
#' 
#' Group of functions for Shannon's entropy calculations
#' 
#' These functions allow the user to calculate the entropy and max entropy of barcode combinations
#'
#' @param frequence probability of occurrence of a given barcode in a set of compatible barcode combinations.
#' @param index_combination matrix of compatible DNA-barcode combinations.
#' @param index_number total number of DNA barcodes provided by the user.
#' @param sample_number number of libraries to be sequenced.
#' @param nb_lane an integer, given the sample_number and multiplex level.
#' 
#' @details 
#' In the context of DNA barcodes, the maximum entropy is reached when the selected barcodes are the least redundant. It corresponds to an ideal case for which consumables are used as much homogeneously as possible in order to extend the life time of expensive library-preparation kits.
#'
#' Case 1: index_number > sample_number
#' In such a case, it would be ideal to assign distinct barcodes to the various libraries to be sequenced (regardless multiplexing), and the resulting Shannon entropy would correspond to the maximum entropy (all barcodes are different): \code{Hmax = log(sample_number)}.
#'
#' Case 2: index_number < sample_number
#' Since there are not as many distinct barcodes as libraries to be sequenced, it is inevitable to use the same barcodes for different libraries. The maximum entropy expression is non-trivial and it depends on both the available number of barcodes and the number of libraries.
#'
#' 
#' @return The return value of each function is briefly documented thereafter for each function.
#' 
#' @examples 
#' # For a set of 48 barcodes and 96 samples (libraries)
#' entropy_max(index_number = 48, sample_number = 96)
#' 
#' @name Utilities_Entropy
NULL

# Shannon's entropy
#' @rdname Utilities_Entropy 
#' @md
#' @section Functions:
#' 
#' * shannon_entropy(): Shannon's entropy
#'  
#' @export
shannon_entropy = function(frequence){
  return(-1 * sum(frequence*log(frequence)))
}

#for a matrix of combination
#' @rdname Utilities_Entropy 
#' @md
#' @section Functions:
#' 
#' * entropy_result(): Shannon's entropy for a matrix of combinations
#'  
#' @export
entropy_result = function (index_combination){
  d = index_combination  %>% table()
  d = d/sum(d)
  entropy = shannon_entropy(d)
  return (entropy)
}

#' @rdname Utilities_Entropy 
#' @md
#' @section Functions:
#' 
#' * entropy_max(): calculte the max entropy given the number of barcodes and libraries
#'  
#' @export
entropy_max = function (index_number,sample_number){
  if(sample_number > index_number){
    k = index_number
    n = sample_number
    entropy =
      - (k - (n %% k)) * (floor(n/k)/n) * log(floor(n/k)/n) -
      (n %% k) * (ceiling(n/k)/n) * log(ceiling(n/k)/n)
    return (entropy)
  }
  else {
    return(log(sample_number))
  }
}


#' @rdname Utilities_Entropy 
#' @md
#' @section Functions:
#' 
#' * comb_result(): Integrate the compatible barcodes combinations with user-deduced experimental constaints: number of libraries to be sequenced, the total number of DNA barcodes provided by the user and the number of lanes in the experiment.
#'  
#' @export
comb_result = function (index_combination, sample_number, index_number, nb_lane){
  # maximal entropy determination for the experiment
  max = entropy_max(index_number,sample_number)
  print(max)
  if(nb_lane < nrow(index_combination)){
    # best combination determination
    one_combination = index_combination [sample(x = 1:nrow(index_combination),nb_lane),]
    an_entropy = entropy_result(one_combination)
    if (an_entropy == max){
      return(one_combination)# we stop if we reach the maximal value of entropy
    }
    else{
      for (i in 1 : 1000){# to debate
        another_combination = index_combination [sample(x = 1:nrow(index_combination),nb_lane),]
        another_entropy = entropy_result(another_combination)
        if(another_entropy == max){
          one_combination = another_combination
          i = 1000 # we stop if we reach the maximal value of entropy
        }
        else if(another_entropy > an_entropy){
          one_combination = another_combination # we keep the best one
          an_entropy = another_entropy # we keep the higher one
        }
      }
    }
  }
  else if (nb_lane >=  nrow(index_combination)){#calcul entopie ????
    {
      one_combination = index_combination[sample(x = 1:nrow(index_combination), nrow(index_combination)),]
      for(i in 2:(floor(nb_lane/nrow(index_combination)))){
        one_combination = rbind(one_combination,
                                index_combination[sample(x = 1:nrow(index_combination), nrow(index_combination)),])
      }
      one_combination = rbind(one_combination,
                              index_combination[sample(x = 1:nrow(index_combination), nb_lane %% nrow(index_combination)),] )
    }
  }
  print(entropy_result(one_combination))
  return (one_combination)
}


# gets the result
#' @rdname Utilities_Entropy 
#' @md
#' @section Functions:
#' 
#' * get_result(): get a compatible set of DNA barcodes for 4-channel sequencing plateforms. DNA-barcode redundancy is reduced by using a Shannon-entropy based optimization.
#'  
#' @export
get_result = function (index, sample_number, multiplexing_level, filter = NULL, metric = 3){
  index = index %>% index_binary_conversion()
  index_combinations = get_combinations(index, multiplexing_level)
  if(!is.null(filter)){
    if(filter == "hamming"){
      hamming_rejection_table = low_hamming_distance(index_distances, metric) 
      index_combinations <<- filter_combinations(index_combinations, hamming_rejection_table)
    }else {
      seqlev_rejection_table = low_seqlev_distance(index_distances, metric) 
      index_combinations <<- filter_combinations(index_combinations, seqlev_rejection_table)
    }
  }
  nb_lane = sample_number %>% as.numeric() / multiplexing_level %>% as.numeric()
  cb = comb_result(index_combinations, sample_number, index_number, nb_lane) %>% as.data.frame()
  result = data.frame(Id = as.vector(cb%>% t() %>% as.vector),
                      Lane = (rep(1:nb_lane, length.out = sample_number, each = multiplexing_level)))
  result$Id = as.character(result$Id)
  return(result)
}

