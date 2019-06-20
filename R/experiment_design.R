#' @title
#' Find a set of barcode combinations with least heterogeneity in 
#' barcode usage for single and dual indexing
#'
#' @description
#' This function uses the Shannon Entropy to identify a set of compatible
#' barcode combinations with least heterogeneity in barcode usage,
#' in the context of single and dual indexing.
#' It performs either an exhaustive or a random-greedy search of compatible
#' DNA-barcode combinations depending on the size of the
#' DNA-barcode population and of the number of samples to be multiplexed.
#'
#' @usage
#' experiment_design(file1, sample_number, mplex_level, platform = 4,
#' file2 = NULL, export = NULL, metric = NULL, d = 3, thrs_size_comb,
#' max_iteration, method)
#'
#' @param file1 The input data file that contains 2 columns separated by a
#' space or a tabulation, namely the sequence identifiers and corresponding
#' DNA sequence
#' @param sample_number Number of libraries to be sequenced.
#' @param mplex_level The number at which the barcodes will be multiplexed.
#' @param file2 The input data file that contains 2 columns separated by
#' a space or a tabulation, namely the sequence identifiers and
#' corresponding DNA sequence; used for dual-indexing, see details below.
#' @param platform An integer representing the number of channels (1, 2, 4)
#' of the desired Illumina platform: 1 for iSeq; 2 for NextSeq, NovaSeq,
#' MiniSeq; 4 for HiSeq and MiSeq. 0 represents any other platform than
#' Illumina.
#' @param export If not NULL, results are saved in a csv file at the
#' specified path.
#' @param metric The type of distance (hamming or seqlev).
#' @param d The minimum value of the distance.
#' @param thrs_size_comb The maximum size of the set of compatible
#' combinations to be used for the greedy optimization.
#' @param max_iteration The maximum number of iterations during the
#' optimizing step.
#' @param method The choice of the greedy search: 'greedy_exchange' or
#' 'greedy_descent'.
#'
#' @details
#' By specifying the total number of libraries and the number of libraries
#' to be multiplexed, this function returns an optimal set of DNA-barcode
#' combinations to be used for sequencing.
#' 
#' In the case of **single indexing**, one must only provide one input file 
#' containing a list of DNA barcodes (file1 argument). The file2 argument being
#' optional, the program runs the optimisation for single indexing if this 
#' argument is left empty. The output shows the sample ID with its respective
#' single barcode.
#' 
#' In the case of **dual indexing**, one must provide two input files
#' containing DNA barcodes as two separate sets of barcodes. The program will 
#' detect these two files and automatically switch to the 'dual indexing' mode. 
#' The program then runs the optimisation for each barcode set separately. 
#' The output shows the sample ID with its respective pair of barcodes.
#' 
#' The inputs of the algorithm are a list of n distinct barcodes,
#' the number N of required libraries, and the multiplex level k; N = ak,
#' where a is the number of lanes of the flow cells to be used for the
#' experiment.
#'
#'
#' * Step 1:
#'
#' This step consists of identifying a set of compatible barcode
#' combinations. Given the number of barcodes and the multiplex level,
#' the total number of barcode combinations (compatible or not) reads:
#' \deqn{{n}\choose{k}}
#' If this number is not too large, the algorithm will perform an
#' exhaustive search and output all compatible combinations of k barcodes.
#' Otherwise, it will proceed by picking up combinations at random, in
#' order to identify a large enough set of compatible barcode combinations.
#'
#' * Step 2:
#'
#' Finds an optimized set of barcode combinations in which barcode
#' redundancy is minimized (see details in
#' \code{\link{optimize_combinations}})
#'
#'
#' @md
#'
#' @return
#' A dataframe containing compatible DNA-barcode combinations organized
#' by lanes of the flow cell.
#'
#' @examples
#' write.table(DNABarcodeCompatibility::IlluminaIndexesRaw,
#' txtfile <- tempfile(), row.names = FALSE, col.names = FALSE, quote=FALSE)
#' experiment_design(file1=txtfile, sample_number=18, mplex_level=3,
#' platform=4)
#'
#'
#' @importFrom stats na.omit
#' @importFrom utils combn read.table write.csv2
#'
#' @export
#'
#'

experiment_design = function (
    file1,
    sample_number,
    mplex_level,
    platform = 4,
    file2 = NULL,
    export = NULL,
    metric = NULL,
    d = 3,
    thrs_size_comb = 120,
    max_iteration = 50,
    method = "greedy_exchange") {
    
    if (is.null(file2)) {
        index_df_1  = file_loading_and_checking(file1)
        if (!is.null(index_df_1)) {
            if (sample_and_multiplexing_level_check(
                sample_number,
                mplex_level) == TRUE) {
                # print("mlx and sample ok")
                result1 = final_result( index_df = index_df_1, 
                                        sample_number = sample_number, 
                                        mplex_level = mplex_level ,
                                        platform = platform,
                                        metric = metric,
                                        d = d,
                                        thrs_size_comb = thrs_size_comb,
                                        max_iteration = max_iteration,
                                        method = method)
                if (!is.null(export)) {
                    write.csv2(result1, file = export)
                }
                return(result1)
            } else{
                stop("The multiplexing level or the sample number is wrong")
            }
        } else{
            stop("An error occured on the file")
        }
    } else{
        index_df_1 = file_loading_and_checking(file1)
            if (!is.null(index_df_1)) {
                index_df_2 = file_loading_and_checking(file2)
                if (!is.null(index_df_2)) {
                    if (sample_and_multiplexing_level_check(
                        sample_number,
                        mplex_level)) {
                result = final_result_dual(index_df_1 = index_df_1,
                                            index_df_2 = index_df_2,
                                            sample_number = sample_number,
                                            mplex_level =   mplex_level,
                                            platform = platform,
                                            metric = metric,
                                            d = d,
                                            thrs_size_comb = thrs_size_comb,
                                            max_iteration = max_iteration,
                                            method = method)
                    
                    if (!is.null(export)) {
                        write.csv2(result, file = export)
                    }
                    return(result)
                } else{
                    stop(paste(
                        "The multiplexing level",
                        "or the sample number is wrong"))
                }
            } else{
                stop("An error occured on the second file")
            }
        } else{
            stop("An error occured on the first file")
        }
    }
}
