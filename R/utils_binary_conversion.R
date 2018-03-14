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
#' @return The return value of each function is briefly documented thereafter for each function.
#' @name Utilities_DNAsequenceBinaryConversion
NULL

#' @rdname Utilities_DNAsequenceBinaryConversion 
#' @md
#' @section Functions:
#' 
#' * sequence_binary_conversion(): converts a DNA sequence into a binary sequence for which A or C becomes 0 and G or T becomes 1
#'  
#' @export

#A function is used for the traduction of a sequence :
sequence_binary_conversion = function(sequence){
  binary_sequence = 
    toupper(sequence) %>% 
    gsub("A|C", "0", .) %>% 
    gsub("G|T", "1", .)
  return(binary_sequence)
}

#' @rdname Utilities_DNAsequenceBinaryConversion 
#' @md
#' @section Functions:
#' 
#' * index_binary_conversion(): binarize each DNA sequence in a vector
#' 
#' @export
#' 
#creates and add to the index the corresponding binary sequence to each row
index_binary_conversion = function(index){
  index = index %>% mutate (binary = sequence_binary_conversion(sequence))
  return (index)
}