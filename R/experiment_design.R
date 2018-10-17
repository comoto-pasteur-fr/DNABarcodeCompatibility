#' @title
#' Find a set of barcode combinations with least redundant barcodes
#' for single and dual indexing
#'
#' @description
#' This function uses the Shannon Entropy to identify a set of compatible
#' barcode combinations with least redundancy between DNA barcodes,
#' in the context of single and dual indexing.
#' It performs either an exhaustive or a random-greedy search of compatible
#' DNA-barcode combinations depending on the size of the
#' DNA-barcode population and of the number of samples to be multiplexed.
#'
#' @usage
#' experiment_design(file1, sample_number, mplex_level, chemistry = 4,
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
#' corresponding DNA sequence; used for dual-indexing.
#' @param chemistry An integer representing the number of channels (1, 2, 4)
#' of the desired Illumina platform.
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
#' chemistry=4)
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
    chemistry = 4,
    file2 = NULL,
    export = NULL,
    metric = NULL,
    d = 3,
    thrs_size_comb = 120,
    max_iteration = 50,
    method = "greedy_exchange") {
    
    if (is.null(file2)) {
        file1  = file_loading_and_checking(file1)
        if (!is.null(file1)) {
            if (sample_and_multiplexing_level_check(
                sample_number,
                mplex_level) == TRUE) {
                # print("mlx and sample ok")
                result1 = final_result(file1, sample_number, 
                                        mplex_level, chemistry, metric, d)
                if (!is.null(export)) {
                    write.csv2(result1, file = export)
                }
                return(result1)
            } else{
                stop("The multiplexing level or the sample number is wrong")
            }
        } else{
            stop("An error occured on the first file")
        }
    } else{
        file1 = file_loading_and_checking(file1)
        if (!is.null(file1)) {
            file2 = file_loading_and_checking(file2)
            if (!is.null(file2)) {
                if (sample_and_multiplexing_level_check(
                    sample_number,
                    mplex_level)) {
                    result1 = get_result(file1,
                                        sample_number,
                                        mplex_level,
                                        chemistry,
                                        metric,
                                        d)
                    result2 = get_result(file2,
                                        sample_number,
                                        mplex_level,
                                        chemistry,
                                        metric,
                                        d)
                    result2 = check_for_duplicate(result1, result2)

                    result1 = left_join(result1,
                                        select(file1, Id, sequence),
                                        by = "Id")
                    # print(result1)
                    result2 = left_join(result2,
                                        select(file2, Id, sequence),
                                        by = "Id")
                    # print(result2)
                    result = data.frame(
                        sample = seq(1,sample_number) %>% as.character(),
                        Lane = result1$Lane %>% as.character(),
                        Id1 = result1$Id %>% as.character(),
                        sequence1 = result1$sequence %>% as.character(),
                        Id2 = result2$Id %>% as.character(),
                        sequence2 = result2$sequence %>% as.character()
                    ) %>%
                        arrange(Lane)
                    
                    if (!is.null(export)) {
                        write.csv2(result1, file = export)
                    }
                    return(result)
                } else{
                    stop("The multiplexing level
                        or the sample number is wrong")
                }
            } else{
                stop("An error occured on the second file")
            }
        } else{
            stop("An error occured on the first file")
        }
    }
}
