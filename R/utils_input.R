# Inputs ------------------------------------------------------------------

#' Group of utility functions
#' 
#' Group of functions for loading and checking DNA barcodes
#' 
#' These functions allow the user to read the input file and check its content 
#'
#' 
#' @param file the name of the file which the data are to be read from. Each row of the table appears as one line of the file. If it does not contain an absolute path, the file name is relative to the current working directory, \code{\link[base]{getwd}}.
#' 
#' @return A data frame containing a representation of the loaded data.
#' @name Utilities_Input

NULL



#' @rdname Utilities_Input 
#' @md
#' @section Functions:
#' 
#' * file_loading_and_checking(): Load and check DNA barcodes
#' 
#' @examples
#' write.table(DNABarcodeCompatibility::illumina, file <- tempfile(), row.names = FALSE, col.names = FALSE, quote=FALSE)
#' index <- file_loading_and_checking(file)
#' index
#' 
#' @export

file_loading_and_checking = function(file){
  index = read_index(file) 
  if (!is.null(index) && index_check(index)){#  if no issue
    index_number <<- nrow(index)
    index  = index %>% mutate (GC_content = get_index_GC_content(index), 
                               homopolymer = get_index_homopolymer(index))
    index_distances <<- index_distance(index)
    return (index)
  }else{
    return(NULL)
  }
}


#' @rdname Utilities_Input 
#' @md
#' @section Functions:
#' 
#' * read_index(): load DNA barcodes only
#' 
#' @export

read_index = function(file) {
  if(!file.exists(file)){
    display_message("Your file doesn't exist, please check the path", call. = FALSE )
    index <<- NULL
  }else{
    index <<- try(as.data.frame(read.table(file, 
                                           header = FALSE, 
                                           sep = "", 
                                           col.names = c("Id", "sequence"),
                                           colClasses = c("character", "character"), 
                                           stringsAsFactors = FALSE)), silent = TRUE)
    if (exists("index")){
      if(class(index) == "try-error"){
        index <<- NULL
        display_message("An error occurred, please check the content of your file")
      }
    }
  }
  return(index)
}


#' @rdname Utilities_Input 
#' @md
#' @section Functions:
#' 
#' * display_message(): display a message on the screen
#' 
#' @export
display_message = function (a_message){
  error_messsage <<- a_message
  print(a_message)
}







