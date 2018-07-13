#' DNABarcodeCompatibility: 
#' to find least-redundant sets of compatible barcodes for multiplex
#' experiments performed on next generation sequencing platforms.
#'
#' 
#' The DNABarcodeCompatibility package provides six functions to load DNA 
#' barcodes, and to generate, filter and optimize sets of barcode combinations 
#' for multiplex sequecing experiments.
#' In particular, barcode combinations are selected to be compatible with 
#' respect to Illumina chemistry constraints, and can be filtered to keep those
#' that are robust against substitution and insertion/deletion errors 
#' thereby facilitating the demultiplexing step.
#' In addtion, the package provides an optimizer function to further favor
#' the selection of compatible barcode combinations with 
#' least redundancy of DNA barcodes.  
#'    
#' @docType package
#' @name DNABarcodeCompatibility-package
#' @import dplyr
#' @import tidyr
#' @import numbers
#' @import purrr
#' @import stringr
#' @import DNABarcodes
NULL