#' @title 
#' Find a set of barcode combinations with least redundant barcodes for single and dual indexing
#'
#' @description 
#' This function uses the Shannon Entropy to identify a set of compatible barcode combinations with least redundancy between DNA barcodes, in the contexts of single and dual indexing.
#' It performs either an exhaustive or a random search of compatible DNA-barcode combinations depending on the size of the DNA-barcode population and of the number of samples to be multiplexed.
#'
#' @usage 
#' experiment_design(file1, sample_number, mplex_level, chemistry = 4,
#'  file2 = NULL, export = NULL, metric = NULL, d = 3)
#'
#' @param file1 Path to the file containing the barcode names and the corresponding DNA sequences.
#' @param sample_number Number of libraries to be sequenced.
#' @param mplex_level The number at which the barcodes will be multiplexed.
#' @param file2 Path to the file containing the barcode names and the corresponding DNA sequences; used for dual-indexing.
#' @param chemistry An integer representing the number of channels (1, 2, 4) of the desired Illumina plateform.
#' @param export If not NULL, results are saved in a csv file at the specified path.
#' @param metric The type of distance (hamming or seqlev).
#' @param d The minimum value of the distance.
#' 
#' @details 
#' By specifying the total number of libraries and the number of libraries to be multiplexed, 
#' this function returns an optimal combination of DNA barcodes to be used for sequencing.
#' 
#' 
#' The inputs of the algorithm are a list of n distinct barcodes, the number N of required libraries, and the multiplex level k; N = ak, where a is the number of lanes of the flow cells to be used for the experiment.
#' 
#' 
#' * Step 1:
#'  
#'  This step consists of identifying a set of compatible barcode combinations. Given the number of barcodes and the multiplex level,
#'  the total number of barcode combinations (compatible or not) reads:  \deqn{{n}\choose{k}}
#'  If this number is not too large, the algorithm will perform an exhaustive search and output all compatible combinations of k barcodes. 
#'  Otherwise, it will proceed by picking up combinations at random, in order to identify a large enough set of compatible barcode combinations. 
#'  
#' * Step 2:
#'  
#'  We then select the N/k compatible combinations to be used in the experiment using a Shannon entropy maximization approach.
#'  It can be shown that the maximum value of the entropy that can be attained for a selection of N barcodes among n, with possible repetitions, reads: 
#'  \deqn{S_{max}=-(n-r)\frac{\lfloor N/n\rfloor}{N}\log(\frac{\lfloor N/n\rfloor}{N})-r\frac{\lceil N/n\rceil}{N}\log(\frac{\lceil N/n\rceil}{N})}
#'  
#'  where r denotes the rest of the division of N by n, and
#'  \deqn{\lfloor N/n\rfloor} and \deqn{\lceil N/n\rceil} denote
#'  the lower and upper integer parts of N/n, respectively.
#'
#'     + Case 1: number of lanes < number of compatible DNA-barcode combinations
#'     
#' This function seeks for compatible DNA-barcode combinations of highest entropy.
#' In brief this function uses a greedy descent algorithm to find an optimized selection. 
#' Note that the resulting optimized selection isn't necessary a globally optimal solution.
#'
#'     + Case 2: number of lanes >= number of compatible DNA-barcode combinations
#'     
#' In such a case, there are not enough compatible DNA-barcode combinations and redundancy is inevitable.
#'
#' @md
#'
#' @return 
#' Returns a dataframe containing compatible DNA-barcode combinations organized by lanes of the flow cell.
#'
#' @examples
#' write.table(DNABarcodeCompatibility::IlluminaIndexes,
#'  txtfile <- tempfile(), row.names = FALSE, col.names = FALSE, quote=FALSE)
#' experiment_design(file1=txtfile, sample_number=18, mplex_level=3, chemistry=4)
#' 
#'
#' 
#' 
#' @export
#' 
#' 

experiment_design = function (file1,
                              sample_number,
                              mplex_level,
                              chemistry = 4,
                              file2 = NULL,
                              export = NULL, 
                              metric = NULL, 
                              d = 3){
  
  if (is.null(file2)){
    file1  = file_loading_and_checking(file1)
    if (!is.null(file1)) {
      if (sample_and_multiplexing_level_check(sample_number, mplex_level) == TRUE){
        print("mlx and sample ok")
        result1 = final_result(file1,sample_number, mplex_level,chemistry, metric, d)
        if(!is.null(export)) {write.csv2(result1, file = export)}
        return(result1)
      }else{
        stop("the multiplexing level and / or the sample number are wrong")}
    }else{
      stop("An error occured on the first file")}
  }else{
    file1 = file_loading_and_checking(file1)
    if (!is.null(file1)){
      file2 = file_loading_and_checking(file2)
      if (!is.null(file2)){
        if(sample_and_multiplexing_level_check(sample_number, mplex_level)){
          result1 = get_result(file1, sample_number, mplex_level,metric, d)
          result2 = get_result(file2, sample_number, mplex_level,metric, d)
          result2 = check_for_duplicate(result1, result2)
          
          result1 = left_join(result1, select(file1, Id, sequence)) 
          print(result1)
          result2 = left_join(result2, select(file2, Id, sequence)) 
          print(result2)
          result = data.frame(sample = 1: sample_number %>% as.character(),
                              Lane = result1$Lane %>% as.character(),
                              Id1 = result1$Id %>% as.character(),
                              sequence1 = result1$sequence %>% as.character(),
                              Id2 = result2$Id %>% as.character(),
                              sequence2 = result2$sequence%>% as.character()) %>% arrange(Lane)
          
          if(!is.null(export)) {write.csv2(result1, file = export)}
          return(result)
        }else{
          stop("the multiplexing level and / or the sample number are wrong")}
      }else{
        stop("An error occured on the second file")
      }
    }else{
      stop("An error occured on the first file")
    }
  }
}
