#' Group of utility functions
#' 
#' Group of functions for finding the possible multiplex levels given the number of libraries
#' 
#' These functions help the user to identify the possible values of multiplex levels given the number of libraries.
#'
#' @param sample_number the number of libraries to be sequenced.
#' @param multiplexing_level the number librairies to be pooled for sequencing
#' 
#' @return The return value of each function is briefly documented thereafter for each function
#' @name Utilities_Multiplexing
NULL


#' @rdname Utilities_Multiplexing 
#' @md
#' @section Functions:
#' 
#' * sample_and_multiplexing_level_check(): check the adequacy between the number of libraries and the multiplexing level; Returns a boolean
#' 
#' @export

sample_and_multiplexing_level_check = function(sample_number,multiplexing_level){
  if (sample_number_check(sample_number)){
    possible_multiplexing_level = multiplexing_level_set(sample_number)
    if (!(multiplexing_level %in% possible_multiplexing_level)){
      display_message( "The sample number isn't a multiple of the multiplexing level. Here are the possible multiplexing levels :")
      display_message(possible_multiplexing_level)
      return(FALSE)
    }else {
      return (TRUE)
    }
  }else {return(FALSE)}}


#' @rdname Utilities_Multiplexing 
#' @md
#' @section Functions:
#' 
#' * sample_number_check(): check if the number of libraries is realistic; Returns a boolean
#' 
#' @export
sample_number_check = function (sample_number){
  # if (!try(x = sample_number, silent = TRUE) ){
  #   print("you have to enter an integer value")
  #   return(FALSE)
  # }
  # if (!exists(deparse(substitute(sample_number)))) {
  #   display_message("you have to enter an integer value1")
  #   return(FALSE)
  # }
  # else 
  if (sample_number!= floor(sample_number)) {
    display_message("you have to enter an integer value2")
    return(FALSE)
  }
  else if (is.na(sample_number)) {
    display_message("you have to enter an integer value3")
    return(FALSE)
  }
  else if (isPrime(sample_number)) {
    display_message("this a prime number, you can't use multiplexing")
    return(FALSE)
  }
  else if (sample_number < 2) {
    message("You need at least 2 samples in order to use multiplexing !")
    return(FALSE)
  }else if (sample_number > 1000){
    display_message("The sample number number is to high, please enter a value under 1000")
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

#' @rdname Utilities_Multiplexing 
#' @md
#' @section Functions:
#' 
#' * multiplexing_level_set(): calculate the multiplex level choices given the number of libraires; Returns a vector of multiplex levels 
#' 
#' @export
multiplexing_level_set = function (sample_number){
  v  = 2 : (sample_number-1)
  multiplexing_level_choices = v[sample_number %% v == 0]
  return (multiplexing_level_choices)
}

