#' @title 
#' Get a large set of compatible combinations.
#'
#' @description 
#' Finds a randomly generated set of at most 1000 combinations
#' of compatible barcodes.
#'
#' @usage 
#' get_random_combinations(index_df, mplex_level, platform)
#'
#' @param index_df A dataframe containing barcodes identifiers, corresponding
#' DNA sequences along with GC content and presence of homopolymers.
#' @param mplex_level The number at which the barcodes will be multiplexed.
#' Illumina recommends to not multiplex more than 96 libraries.
#' @param platform An integer representing the number of channels (1, 2, 4)
#' of the desired Illumina platform: 1 for iSeq; 2 for NextSeq, NovaSeq,
#' MiniSeq; 4 for HiSeq and MiSeq. 0 represents any other platform than
#' Illumina.
#'
#' @details 
#' This function is suited if the total number of possible combinations is too
#' high for an exhaustive search to be possible in a reasonable amount of time.
#'
#' @return 
#' A matrix containing the identifiers of compatible barcode combinations.
#'
#' @examples
#' get_random_combinations(DNABarcodeCompatibility::IlluminaIndexes, 3, 4)
#' 
#'
#' @seealso 
#' \code{\link{get_all_combinations}}
#' 
#' @export
#' 

get_random_combinations = function(index_df, mplex_level, platform){
    check_input_dataframe(  index_df, 
                            c("Id", "sequence", "GC_content", "homopolymer"))
    if (is.numeric(mplex_level)){
        if(mplex_level <= 96){
            if (mplex_level <= nrow(index_df)){
                combinations_m = NULL
                if (platform == 4 ){
                    combinations_m = 
                        get_random_combinations_4_channel(  index_df,
                                                            mplex_level)
                } else if (platform == 2){
                    combinations_m = 
                        get_random_combinations_2_channel(  index_df,
                                                            mplex_level)
                } else if (platform == 1){
                    combinations_m = 
                        get_random_combinations_1_channel(  index_df,
                                                            mplex_level)
                } else if (platform == 0){
                    combinations_m = 
                        get_random_combinations_0_channel(  index_df,
                                                            mplex_level)
                } else {
                    display_message(paste(  "Please choose a correct platform",
                                            "for your experiment"))
                    return(NULL)
                }
                if(nrow(as.matrix(combinations_m)) == 0){
                    display_message(paste(
                        "The programm didn't find any compatible", 
                        "combination, please check your index list",
                        "and your search parameters"))
                    return (NULL)
                }else{
                    return (combinations_m)
                }
            } else {display_message(paste(  "The value of mplex_level",
                                            " should not be higher than the",
                                            "number of available barcodes"))
            }
                
        }else{display_message(
            "The value of mplex_level should not be higher than 96")
        }
    }else{
        display_message("Please enter a number as mplex_level")
        return (NULL)
    }
}
