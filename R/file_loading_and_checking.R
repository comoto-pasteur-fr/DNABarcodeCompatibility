#' @title
#' Loading and checking DNA barcodes.
#'
#' @description
#' Loads the file containing DNA barcodes and analyze barcode content.
#'
#' @usage
#' file_loading_and_checking(file)
#'
#' @param file The input data file that contains 2 columns separated by a space
#' or a tabulation, namely the sequence identifiers and
#' corresponding DNA sequence.
#'
#' @details
#' This function loads the DNA barcodes from the input file and checks barcodes
#' for unicity (identifier and sequence), DNA content, and equal size.
#' It also calculates the fraction of G and C relative to A and T, as referred
#' to as "GC content", and it detects the presence of
#' homopolymers of length >= 3.
#'
#' @return
#' A dataframe containing sequence identifiers, nucleotide sequence, GC content,
#' presence of homopolymers.
#'
#' @examples
#' write.table(DNABarcodeCompatibility::IlluminaIndexesRaw,
#' txtfile <- tempfile(), row.names = FALSE, col.names = FALSE, quote=FALSE)
#' file_loading_and_checking(txtfile)
#'
#' @export
#' 

file_loading_and_checking = function(file){
    index = read_index(file) 
    if (!is.null(index) && index_check(index)){#  if no issue
        index_number <- nrow(index)
        index  = index %>% mutate (GC_content = get_index_GC_content(index), 
                                    homopolymer = get_index_homopolymer(index))
        return (index)
    }else{
        return(NULL)
    }
}