#final function
#' @title 
#' Get compatible combinations of DNA barcodes
#'
#' @description 
#' Performs either an exhaustive or a random search of compatible DNA-barcode combinations depending on the size of the DNA-barcode population and of the number of samples to be multiplexed.
#'
#' @usage 
#' experiment_design(file1, sample_number, multiplexing_level, file2 = NULL, export = NULL,  filter = NULL, metric = 3)
#'
#' @param file1 path to the file containing the barcode names and the corresponding DNA sequences.
#' @param sample_number number of libraries to be sequenced
#' @param multiplexing_level the number of samples to be multiplexed for sequencing
#' @param file2 path to the file containing the barcode names and the corresponding DNA sequences; used for dual-indexing.
#' @param export TO BE DOCUMENTED
#' @param filter TO BE DOCUMENTED
#' @param metric minimal distance allowed between barcodes TO BE RENAMED dist
#' 
#' 
#' @details 
#' Illumina 4-channel sequencers use a green laser to read G/T nucleotides and a red laser to
#' read A/C nucleotides. With each sequencing cycle at least one of the two
#' nucleotides for each color channel must be read to ensure proper registration.
#' Since the DNA-barcode sequence of pooled, multiplexed libraries will be read
#' simultaneously, it is important to maintain color balance for each base of the
#' DNA-barcode sequences in the pooled library. Otherwise, DNA-barcode sequencing will fail due
#' to registration failure.
#'
#' By specifying the total number of libraries and the number of libraries to be multiplexed, 
#' this function returns an optimal combination of DNA barcodes to be used for sequencing.
#' This function first converts DNA barcodes into binary words: A and C will be converted in "0", G and T in "1", see \code{\link{sequence_binary_conversion}}. 
#' 
#' Second, it selects a series of combinations of binarized barcodes that are compatible for 2-channel sequencing, see \code{\link{is_a_good_combination}}, \code{\link{get_combinations}}.
#' 
#' Third, it reduces redundancy between DNA barcodes so that consumables can be used as much homogeneously as possible in order to extend the life time of expensive library-preparation kits. To this end, the maximum entropy of the distribution of barcodes is calculated given the available number of barcodes and libraries, see \code{\link{entropia_result}}, \code{\link{entropia_max}}. The maximum-entropy value sets the upper limit to be reached or approached when calculating the entropy of concurrent sets of compatible DNA barcodes.  
#'
#' Case 1: number of lanes < number of compatible DNA-barcode combinations
#' This function seeks for compatible DNA-barcode combinations of highest entropy. The search stops either if the maximum-entropy value is reached for a given set of DNA-barcode combinations or after 1000 iterations but returning the set of DNA-barcode combinations of the highest entropy.
#' The criterium of 1000 iterations has been choosen based on simulations.
#'
#' Case 2: number of lanes > number of compatible DNA-barcode combinations
#' In such a case, there are not enough compatible DNA-barcode combinations and redundancy is inevitable. This may dramatically decrease the life time of the library-preparation kits leading to the necessity to buy new expensive kits just for providing enough of some DNA barcodes.
#'
#' TO BE DOCUMENTED: 2-channel sequencers
#'
#' @return 
#' Returns a dataframe containing compatible DNA-barcode combinations organized by lanes of the flow cell.
#'
#' @examples
#' write.table(DNABarcodeCompatibility::illumina, file <- tempfile(), row.names = FALSE, col.names = FALSE, quote=FALSE)
#' results <- experiment_design(file, 6, 3)
#' 
#'
#' @seealso 
#' \code{\link{get_result}}
#' 
#' @export
#' 
experiment_design = function (index1,
                              sample_number,
                              multiplexing_level,
                              index2 = NULL,
                              export = NULL, 
                              filter = NULL, 
                              metric = 3){
  if (is.null(index2)){
    index1  = file_loading_and_checking(index1)
    if (!is.null(index1)) {
      if (sample_and_multiplexing_level_check(sample_number, multiplexing_level) == TRUE){
        print("mlx and sample ok")
        result1 = get_result(index1,sample_number, multiplexing_level, filter,metric)
        result1 = data.frame(sample = 1: sample_number %>% as.character(),
                             Lane = result1$Lane %>% as.character(),
                             Id = result1$Id %>% as.character(),
                             stringsAsFactors = FALSE)
        result1 = left_join(result1, select(index1, Id, sequence),by="Id") 
        if(!is.null(export)) {write.csv2(result1, file = export)}
        return(result1)
      }else{
        stop("the multiplexing level and / or the sample number are wrong")}
    }else{
      stop("An error occured on the first file")}
  }else{
    index1 = file_loading_and_checking(index1)
    if (!is.null(index1)){
      index2 = file_loading_and_checking(index2)
      if (!is.null(index2)){
        if(sample_and_multiplexing_level_check(sample_number, multiplexing_level)){
          result1 = get_result(index1, sample_number, multiplexing_level,filter,metric)
          result2 = get_result(index2, sample_number, multiplexing_level,filter,metric)
          result2 = check_for_duplicate(result1, result2)
          
          result1 = left_join(result1, select(index1, Id, sequence)) 
          print(result1)
          result2 = left_join(result2, select(index2, Id, sequence)) 
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


check_for_duplicate = function(result1, result2){
  check = data.frame(Id1 = result1$Id, Id2 = result2$Id)
  print(check)
  if (anyDuplicated(check) != 0){
    d = anyDuplicated(check)
    for(i in 1 : length(d)){
      id = result2[d[i],]$Id
      lane_to_change = result2 %>% filter(Lane == result2[d[i],]$Lane)
      lanes_to_keep  = result2 %>% filter(Lane != result2[d[i],]$Lane)
      j = which(lane_to_change$Id == id) %>% as.numeric()
      if (j<nrow(lane_to_change)){
        temp = lane_to_change[j,] 
        lane_to_change[j,] = lane_to_change[j+1,] 
        lane_to_change[j+1,] = temp
      }else {
        temp = lane_to_change[j,] 
        print(temp)
        lane_to_change[j,] = lane_to_change[1,] 
        lane_to_change[1,] = temp
      }
      result = bind_rows(lane_to_change,lanes_to_keep) %>% arrange(Lane)
    }
    print(result2)
    return (result)
  }else{
    return (result2)
  }
}

