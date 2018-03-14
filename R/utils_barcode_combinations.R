#' Group of utility functions
#' 
#' Group of functions for loading and checking DNA barcodes
#' 
#' These functions allow the user to read the input file and check its content 
#'
#' 
#' @param sequence the name of the file which the data are to be read from. Each row of the table appears as one line of the file. If it does not contain an absolute path, the file name is relative to the current working directory, \code{\link[base]{getwd}}.
#' @param index a dataframe containing a representation of the DNA-barcode dataset.
#' 
#' @return A data frame containing a representation of the loaded data.
#' @name Utilities_DNAbarcodeCombinations
NULL

#' @rdname Utilities_DNAbarcodeCombinations 
#' @md
#' @section Functions:
#' 
#' * binary_word_into_numeric(): conversion of the string sequence into a vector of numeric
#' 
#' @examples
#' 
#' @export

#conversion of the string sequence into a vector of numeric
binary_word_into_numeric = function (binary_word){
  as.numeric(unlist(strsplit(as.character(binary_word),"")))
}

# conversion of each index sequence combination into a matrix
#' @rdname Utilities_DNAbarcodeCombinations 
#' @md
#' @section Functions:
#' 
#' * vectors_into_matrix(): conversion of each index sequence combination into a matrix
#' 
#' @examples
#' 
#' @export
vectors_into_matrix = function (binary_word){
  m =  mapply(binary_word_into_numeric,binary_word)
  return(m)
}


#test if a column/line of a index combination is correct
#' @rdname Utilities_DNAbarcodeCombinations 
#' @md
#' @section Functions:
#' 
#' * any_different(): test if a column/line of a index combination is correct
#' 
#' @examples
#' 
#' @export
any_different = function(binary_sequence){
  if (length(unique(binary_sequence)) > 1){
    return (TRUE)
  }
  else{
    return (FALSE)
  }
}

# check if a combination of indexes is correct
#' @rdname Utilities_DNAbarcodeCombinations 
#' @md
#' @section Functions:
#' 
#' * is_a_good_combination(): check if a combination of indexes is correct
#' 
#' @examples
#' 
#' @export
is_a_good_combination = function (combination_matrix){
  all_combinations = vectors_into_matrix(combination_matrix)
  results = prod(apply(all_combinations,1,any_different))
  return(results)
}

# keeps only the good ones :
#' @rdname Utilities_DNAbarcodeCombinations 
#' @md
#' @section Functions:
#' 
#' * list_of_good_combinations(): keep only good combinations
#' 
#' @examples
#' 
#' @export
list_of_good_combinations = function (matrix_id){
  list = apply(matrix_id,2,is_a_good_combination)
  return(list)
}

##super fast and furious
#matches an id to its binary_sequence
#' @rdname Utilities_DNAbarcodeCombinations 
#' @md
#' @section Functions:
#' 
#' * id_into_binary_sequence(): matches an id to its binary_sequence
#' 
#' @examples
#' 
#' @export
id_into_binary_sequence = function (index_id_combination, index){
  index_rows = subset(index, Id == index_id_combination)
  index_binary_sequence_combination = as.character(index_rows$binary)
  return (index_binary_sequence_combination)
}


#matches the matrix_id to the matrix_binary_sequence
#' @rdname Utilities_DNAbarcodeCombinations 
#' @md
#' @section Functions:
#' 
#' * matrix_id_to_binary_sequence(): matches the matrix_id to the matrix_binary_sequence
#' 
#' @examples
#' 
#' @export
matrix_id_to_binary_sequence = function(matrix_id, index){
  m = matrix_id
  n = nrow(m)
  m[1:n,] = sapply(X = m[1:n,], FUN = id_into_binary_sequence, index)
  return (m)
}

# low distance tab = hamming rejection table or seq lev rejection table
#' @rdname Utilities_DNAbarcodeCombinations 
#' @md
#' @section Functions:
#' 
#' * filter_combinations(): low distance tab = hamming rejection table or seq lev rejection table
#' 
#' @examples
#' 
#' @export
filter_combinations = function(combinations, low_distance_tab){
  
  id1 = combinations %>% as.data.frame() %>% select_all() == unlist(low_distance_tab$Id1)
  id1 = apply(id1, 1, any)
  id2 = combinations%>% as.data.frame() %>% select_all() == unlist(low_distance_tab$Id2)
  id2 = apply(id2, 1, any)
  
  to_remove = id1 * id2
  return(combinations [!as.logical(to_remove),])
  
}
