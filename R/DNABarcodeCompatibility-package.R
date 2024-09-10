#' DNABarcodeCompatibility: 
#' to find optimised sets of compatible barcodes with least heterogeneity in 
#' barcode usage for multiplex
#' experiments performed on next generation sequencing platforms.
#'
#' 
#' The DNABarcodeCompatibility package provides six functions to load DNA 
#' barcodes, and to generate, filter and optimise sets of barcode combinations 
#' for multiplex sequencing experiments.
#' In particular, barcode combinations are selected to be compatible with 
#' respect to Illumina chemistry constraints, and can be filtered to keep those
#' that are robust against substitution and insertion/deletion errors 
#' thereby facilitating the demultiplexing step.
#' In addition, the package provides an optimiser function to further favor
#' the selection of compatible barcode combinations with least 
#' heterogeneity in barcode usage.
#'    
#' @docType package
#' @name DNABarcodeCompatibility-package
#' @import dplyr
#' @import tidyr
#' @import numbers
#' @import purrr
#' @import stringr
#' @import methods
NULL