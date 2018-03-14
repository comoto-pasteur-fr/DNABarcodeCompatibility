# Inputs ------------------------------------------------------------------

#' Group of utility functions
#' 
#' Group of functions for analyzing and checking DNA barcodes
#' 
#' These functions allow the user to check the DNA content of the barcodes, and in particular the GC content, the homopolymer length and the distance between barcodes.
#'
#' 
#' @param index a dataframe containing a representation of the DNA-barcode dataset.
#' @param sequence a DNA sequence as string
#' @param index_distances a dataframe containing the calculted distances between pairs de DNA barcodes
#' @param metrict threshold distance between barcodes
#' 
#' @return The return value of each function is briefly documented thereafter for each function
#' @name Utilities_DNAbarcodeAnalysis
NULL


#' @rdname Utilities_DNAbarcodeAnalysis 
#' @md
#' @section Functions:
#' 
#' * unicity_check(): check the unicity of barcodes sequences and identifiers; Returns a boolean
#' 
#' 
#' @export

unicity_check = function(index){
  index$sequence <<- toupper(index$sequence)
  if(index$Id %>% anyDuplicated() != 0){#checks if the index Ids are unique
    v = paste("two duplicated indexes IDs,check row n°", anyDuplicated(index$Id))
    display_message(v)
    return(FALSE)
  }
  else if(index$sequence %>% anyDuplicated() != 0){#checks if the index sequences are unique
    v = paste("At least one of your sequences is not unique, check row n°", anyDuplicated(index$sequence))
    display_message(v) 
    return(FALSE) 
  }
  else{return (TRUE)}
}

#check the character for one sequence
#' @rdname Utilities_DNAbarcodeAnalysis 
#' @md
#' @section Functions:
#' 
#' * character_check(): check the character composition of one sequence; Returns a boolean
#' 
#' 
#' @export

character_check = function(sequence){
  wrong_letters = LETTERS[!(LETTERS) %in% c("A", "G", "C", "T")]
  check = all(!str_detect(sequence, wrong_letters))
  return(check)
}

#checks the characters of each sequence
#' @rdname Utilities_DNAbarcodeAnalysis 
#' @md
#' @section Functions:
#' 
#' * sequences_character_check(): check the character composition of each sequence in a vector of sequences; 
#' 
#' 
#' @export
sequences_character_check = function(sequences){
  return(map_lgl(sequences, character_check))
}

#checks the length of each sequence
#' @rdname Utilities_DNAbarcodeAnalysis 
#' @md
#' @section Functions:
#' 
#' * length_check(): check the length of each sequence; returns a boolean
#' 
#' 
#' @export
length_check = function(sequences){
  return(str_length(sequences) == str_length(sequences[1]))
}

#checks for character type and sequence  length
#' @rdname Utilities_DNAbarcodeAnalysis 
#' @md
#' @section Functions:
#' 
#' * character_and_length_check(): check for character type and sequence length; returns a boolean
#' 
#' 
#' @export
character_and_length_check = function(index){
  c_check =  sequences_character_check(index$sequence)
  if (length (which(c_check == FALSE)) > 0 ){
    display_message(paste("Your sequence contains a wrong charater, row n° : ", which(c_check == FALSE)))
    return(FALSE)
  }
  else {
    l_check = length_check(index$sequence)
    if (length(which(l_check == FALSE)) > 0){
      wrong_length = as.numeric(
        names(
          table(
            nchar
            (index$sequence)))[as.numeric(which(table(nchar(index$sequence))== min(table(nchar(index$sequence)))))])
      c = which(nchar(index$sequence) == wrong_length)
      display_message(paste("the indexes have not the same size, check out row n°", c))
      
      return(FALSE)
    }else{
      return(TRUE)
    }
  }
}



## Check the index for possible issues ####
#' @rdname Utilities_DNAbarcodeAnalysis 
#' @md
#' @section Functions:
#' 
#' * index_check(): check the index for possible issues; returns a boolean
#' 
#' 
#' @export
index_check = function(index){
  return (all(unicity_check(index), 
              character_and_length_check(index)))
}

## Calculate GC content of one sequence ####
#' @rdname Utilities_DNAbarcodeAnalysis 
#' @md
#' @section Functions:
#' 
#' * get_sequence_GC_content(): calculate GC content of one sequence; returns GC content as numerics
#' 
#' 
#' @export
get_sequence_GC_content =  function (sequence){
  GC_content = str_count(sequence,pattern = ("G|C"))/nchar(sequence) * 100
  return(round(GC_content, digits = 2))
}

## calculate GC content of each sequence and assign the sequence identifier ####
#' @rdname Utilities_DNAbarcodeAnalysis 
#' @md
#' @section Functions:
#' 
#' * get_index_GC_content(): calculate GC content of each sequence and assign the sequence identifier; returns TO BE COMPLETED
#' 
#' @export
get_index_GC_content = function(index){
  index_GC_content = map_dbl(index$sequence, get_sequence_GC_content)
  return (index_GC_content)
}

## calculate homopolymer length of one sequence ####
#' @rdname Utilities_DNAbarcodeAnalysis 
#' @md
#' @section Functions:
#' 
#' * get_sequence_homopolymer(): calculate homopolymer length of one sequence; returns a boolean
#' 
#' @export
get_sequence_homopolymer = function(sequence){
  if(length(str_which(sequence, "A{3,}|G{3,}|C{3,}|T{3,}")) > 0 ){
    return (TRUE)
  }else{return (FALSE)}
}

#' @rdname Utilities_DNAbarcodeAnalysis 
#' @md
#' @section Functions:
#' 
#' * get_index_homopolymer(): calculate homopolymer length of each sequence and assign the sequence identifier, returns 
#' 
#' @export
get_index_homopolymer = function(index){
  homopolymer  = map_lgl(index$sequence, get_sequence_homopolymer)
}


# generates all possible couples and caculates DNAbarcodes's hamming and seqlev distances
#' @rdname Utilities_DNAbarcodeAnalysis 
#' @md
#' @section Functions:
#' 
#' * index_distance(): generates all possible barcode pairs and caculates hamming and seqlev distances; returns a dataframe
#' 
#' @export
index_distance = function (index){
  index_couple = combn(index$sequence,2)  %>% t()%>%  as.data.frame(., stringsAsFactors = FALSE)
  index_couple = index_couple %>% mutate(n = 1 : nrow(index_couple))
  index_couple = index_couple %>% group_by(n) %>% mutate (hamming = DNABarcodes::distance(V1, V2, metric = "hamming"),
                                                          seqlev = DNABarcodes::distance (V1, V2, metric = "seqlev"))
  return(index_couple)
}

# couples with hamming distance under threshold
#' @rdname Utilities_DNAbarcodeAnalysis 
#' @md
#' @section Functions:
#' 
#' * low_hamming_distance():  TO BE DOCUMENTED
#' 
#' 
#' @export
low_hamming_distance = function(index_distances, metric){
  
  i_d = index_distances %>% filter(hamming < metric)
  i_d = i_d %>% group_by(n) %>% mutate(Id1 = index$Id[which(V1 == index$sequence)],
                                       Id2 = index$Id[which(V2 == index$sequence)])
  lhd = i_d %>% ungroup() %>% select(Id1, Id2) 
  return(lhd)
}

#couples with seqlev distance under threshold
#' @rdname Utilities_DNAbarcodeAnalysis 
#' @md
#' @section Functions:
#' 
#' * low_hamming_distance():  TO BE DOCUMENTED
#' 
#' 
#' @export
low_seqlev_distance = function(index_distances, metric){
  i_d = index_distances %>% filter(seqlev < metric)
  i_d = i_d %>% group_by(n) %>% mutate(Id1 = index$Id[which(V1 == index$sequence)],
                                       Id2 = index$Id[which(V2 == index$sequence)]) 
  lsd = i_d %>% ungroup() %>% select(Id1, Id2)
  return(lsd)
}
