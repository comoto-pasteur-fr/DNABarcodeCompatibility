#' @title 
#' Display message.
#'
#' @description 
#' Display a message for the user and export it in the global environment.
#'
#' @usage 
#' display_message(a_message)
#'
#' @param a_message A character string containg the message to be displayed.
#'
#' @details 
#' The error_messsage variable that contains the message is exported in the global environment. This variable is read by the Java GUI.
#'
#' @return 
#' Prints the message and returns the error_messsage character string
#'
#' @examples
#' display_message("Hello world")
#' error_messsage
#'
#' 
#' @export
#' 

error_messsage = ""
display_message = function (a_message){
  error_messsage <<- a_message
  print(a_message)
}