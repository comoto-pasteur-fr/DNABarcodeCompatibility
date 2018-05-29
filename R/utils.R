
library("dplyr")
library("tidyr")
library("numbers")
library("purrr")
library ("stringr")
library("DNABarcodes")


error_messsage = ""


# Inputs ------------------------------------------------------------------


## reading of the file and creation of a data.frame containing the index Ids and corresponding sequence
read_index = function(file) {
  if(!file.exists(file)){
    display_message("Your file doesn't exist, please check the path", call. = FALSE )
    index <<- NULL
  }else{
    index <<- try(as.data.frame(read.table(file, 
                                           header = FALSE, 
                                           sep = "", 
                                           col.names = c("Id", "sequence"),
                                           colClasses = c("character", "character"), 
                                           stringsAsFactors = FALSE)), silent = TRUE)
    if (exists("index")){
      if(class(index) == "try-error"){
        index <<- NULL
        display_message("An error occurred, please check the content of your file")
      }
    }
  }
  return(index)
}


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
character_check = function(sequence){
  wrong_letters = LETTERS[!(LETTERS) %in% c("A", "G", "C", "T")]
  check = all(!str_detect(sequence, wrong_letters))
  return(check)
}

#checks the characters of each sequence
sequences_character_check = function(sequences){
  return(map_lgl(sequences, character_check))
}

#checks the length of each sequence
length_check = function(sequences){
  return(str_length(sequences) == str_length(sequences[1]))
}

#checks for character type and sequence  length
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



## Check the index for possible issues
index_check = function(index){
  return (all(unicity_check(index), 
              character_and_length_check(index)))
}

get_sequence_GC_content =  function (sequence){
  GC_content = str_count(sequence,pattern = ("G|C"))/nchar(sequence) * 100
  return(round(GC_content, digits = 2))
}

get_index_GC_content = function(index){
  index_GC_content = map_dbl(index$sequence, get_sequence_GC_content)
  return (index_GC_content)
}

get_sequence_homopolymer = function(sequence){
  if(length(str_which(sequence, "A{3,}|G{3,}|C{3,}|T{3,}")) > 0 ){
    return (TRUE)
  }else{return (FALSE)}
}

get_index_homopolymer = function(index){
  homopolymer  = map_lgl(index$sequence, get_sequence_homopolymer)
}










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

multiplexing_level_set = function (sample_number){
  v  = 2 : (sample_number-1)
  multiplexing_level_choices = v[sample_number %% v == 0]
  return (multiplexing_level_choices [multiplexing_level_choices< 96])
}

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


file_loading_and_checking = function(file){
  index = read_index(file) 
  if (!is.null(index) && index_check(index)){#  if no issue
    index_number <<- nrow(index)
    index  = index %>% mutate (GC_content = get_index_GC_content(index), 
                               homopolymer = get_index_homopolymer(index))
    #index_distances <<- index_distance(index)
    return (index)
  }else{
    return(NULL)
  }
}


# Binary conversion -------------------------------------------------------

# function used for the traduction of the nucleotide sequences into binary sequences according to the chemistry
# index_binary_conversion = function(index_df, chemistry){
#   if (chemistry == 4 ){
#     index_df = index_binary_conversion_4_channel(index_df)
#   } else if (chemistry == 2){
#     index_df = index_binary_conversion_2_channel(index_df)
#   } else if (chemistry == 1){
#     index_df = index_binary_conversion_1_channel(index_df)
#   } else {
#     display_message("Please choose a correct chemistry for your experiment ")
#   }
#   return (index_df)
# }


# sequence traduction for a 4_channel chemistry :
sequence_binary_conversion_4_channel = function(sequence){
  binary_sequence = 
    toupper(sequence) %>% 
    gsub("A|C", "0", .) %>% 
    gsub("G|T", "1", .)
  return(binary_sequence)
}
# index traduction for a 4_channel chemistry :
index_binary_conversion_4_channel = function(index){
  index = index %>% mutate (binary_4 = sequence_binary_conversion_4_channel(sequence))
  return (index)
}

# sequence traduction for a 2_channel chemistry for image 1 :
sequence_binary_conversion_2_channel_1 = function(sequence){
  binary_sequence = 
    toupper(sequence) %>% 
    gsub("G|C", "0", .) %>% 
    gsub("A|T", "1", .)
  return(binary_sequence)
}

# sequence traduction for a 2_channel chemistry for image 2 :
sequence_binary_conversion_2_channel_2 = function(sequence){
  binary_sequence = 
    toupper(sequence) %>% 
    gsub("G|T", "0", .) %>% 
    gsub("A|C", "1", .)
  return(binary_sequence)
}

# index traduction for a 2_channel chemistry :
index_binary_conversion_2_channel = function(index){
    index = index %>% mutate (binary_2_image_1 = sequence_binary_conversion_2_channel_1(sequence), 
                              binary_2_image_2 = sequence_binary_conversion_2_channel_2(sequence))
    return (index)
}


# sequence traduction for a 1_channel chemistry for image 1 :
sequence_binary_conversion_1_channel_1 = function(sequence){
  binary_sequence = 
    toupper(sequence) %>% 
    gsub("G|C", "0", .) %>% 
    gsub("A|T", "1", .)
  return(binary_sequence)
}

# sequence traduction for a 1_channel chemistry for image 2 :
sequence_binary_conversion_1_channel_2 = function(sequence){
  binary_sequence = 
    toupper(sequence) %>% 
    gsub("A|G", "0", .) %>% 
    gsub("C|T", "1", .)
  return(binary_sequence)
}
index_binary_conversion_1_channel = function(index){
  index = index %>% mutate (binary_1_image_1 = sequence_binary_conversion_1_channel_1(sequence), 
                            binary_1_image_2 = sequence_binary_conversion_1_channel_2(sequence))
  return (index)
}


# Compatibility -----------------------------------------------------------



#conversion of the string sequence into a vector of numeric
binary_word_into_numeric = function (binary_word){
  as.numeric(unlist(strsplit(as.character(binary_word),"")))
}

#conversion of each index sequence combination into a matrix
vectors_into_matrix = function (binary_word){
  m =  mapply(binary_word_into_numeric,binary_word)
  return(m)
}


#test if a column/line of a index combination is correct
any_different = function(binary_sequence){
  if (length(unique(binary_sequence)) > 1){
    return (TRUE)
  }else{
    return (FALSE)
  }
}

has_signal_in_both_channels = function(colored_sequence){
  if(any(as.logical(colored_sequence))){
    return (TRUE)
  }else{
    return (FALSE)
  }
}


# check if a combination of indexes is correct
is_a_good_combination = function (combination_matrix){
  all_combinations = vectors_into_matrix(combination_matrix)
  results = prod(apply(all_combinations,1,any_different))
  return(results)
}

is_a_good_combination_2 = function (combination_matrix){
  all_combinations = vectors_into_matrix(combination_matrix)
  results = prod(apply(all_combinations,1,has_signal_in_both_channels))
  return(results)
}

# keeps only the good ones :
list_of_good_combinations = function (matrix_id){
  list = apply(matrix_id,2,is_a_good_combination)
  return(list)
}

list_of_good_combinations_2 = function (matrix_id){
  list = apply(matrix_id,2,is_a_good_combination_2)
  return(list)
}

##super fast and furious
#matches an id to its binary_sequence
id_into_4_channel_binary_sequence = function (index_id_combination, index_df){
  index_rows = subset(index_df, Id == index_id_combination)
  index_binary_sequence_combination = as.character(index_rows$binary_4)
  return (index_binary_sequence_combination)
}

# matches an id to its binary_sequence for 2_channel image 1
id_into_2_channel_image_1_binary_sequence = function (index_id_combination, index_df){
  index_rows = subset(index_df, Id == index_id_combination)
  index_binary_sequence_combination = as.character(index_rows$binary_2_image_1)
  return (index_binary_sequence_combination)
}

# matches an id to its binary_sequence for 2_channel image 2
id_into_2_channel_image_2_binary_sequence = function (index_id_combination, index_df){
  index_rows = subset(index_df, Id == index_id_combination)
  index_binary_sequence_combination = as.character(index_rows$binary_2_image_2)
  return (index_binary_sequence_combination)
}

# matches an id to its binary_sequence for 1_channel image 1
id_into_1_channel_image_1_binary_sequence = function (index_id_combination, index_df){
  index_rows = subset(index_df, Id == index_id_combination)
  index_binary_sequence_combination = as.character(index_rows$binary_1_image_1)
  return (index_binary_sequence_combination)
}

# matches an id to its binary_sequence for 1_channel image 2
id_into_1_channel_image_2_binary_sequence = function (index_id_combination, index_df){
  index_rows = subset(index_df, Id == index_id_combination)
  index_binary_sequence_combination = as.character(index_rows$binary_1_image_2)
  return (index_binary_sequence_combination)
}



#matches an id to its binary_sequence
# id_into_binary_sequence_slow = function (an_index_id_combination, index){
#   result = as.matrix(index %>% filter (Id %in% an_index_id_combination) %>% select(binary))
#   return (result)
# }

# matches the matrix_id to the matrix_binary_sequence
matrix_id_to_binary_sequence = function(matrix_id, index_df){
  m = matrix_id
  n = nrow(m)
  m[1:n,] = sapply(X = m[1:n,], FUN = id_into_4_channel_binary_sequence, index_df)
  return (m)
}

# matches the matrix_id to the matrix_2_channel_immage_1_binary_sequence
matrix_id_to_2_channel_image_1_binary_sequence = function(matrix_id, index_df){
  m = matrix_id
  n = nrow(m)
  m[1:n,] = sapply(X = m[1:n,], FUN = id_into_2_channel_image_1_binary_sequence, index_df)
  return (m)
}

# matches the matrix_id to the matrix_2_channel_immage_2_binary_sequence
matrix_id_to_2_channel_image_2_binary_sequence = function(matrix_id, index_df){
  m = matrix_id
  n = nrow(m)
  m[1:n,] = sapply(X = m[1:n,], FUN = id_into_2_channel_image_2_binary_sequence, index_df)
  return (m)
}

# matches the matrix_id to the matrix_1_channel_immage_1_binary_sequence
matrix_id_to_1_channel_image_1_binary_sequence = function(matrix_id, index_df){
  m = matrix_id
  n = nrow(m)
  m[1:n,] = sapply(X = m[1:n,], FUN = id_into_1_channel_image_1_binary_sequence, index_df)
  return (m)
}

# matches the matrix_id to the matrix_1_channel_immage_2_binary_sequence
matrix_id_to_1_channel_image_2_binary_sequence = function(matrix_id, index_df){
  m = matrix_id
  n = nrow(m)
  m[1:n,] = sapply(X = m[1:n,], FUN = id_into_1_channel_image_2_binary_sequence, index_df)
  return (m)
}

# get all compatible combinations of an index for 4_channel chemistry
get_all_combinations_4_channel = function(index_df, multiplexing_level){
  index_df = index_binary_conversion_4_channel(index_df)
  index_df = index_df %>% arrange(Id)
  matrix_id = combn(x = index_df$Id, m = multiplexing_level)
  matrix_binary_sequence = matrix_id_to_binary_sequence(matrix_id = matrix_id, index_df = index_df)#matches Ids to binary sequences
  logical_rigth_combination = as.logical(x = list_of_good_combinations(m = matrix_binary_sequence))
  list_of_all_combinations = matrix_id[, logical_rigth_combination] %>% t()
  return(list_of_all_combinations)
}

# get all compatible combinations of an index for 2_channel chemistry
get_all_combinations_2_channel = function(index_df, multiplexing_level){
  index_df = index_binary_conversion_2_channel(index_df)
  index_df = index_df %>% arrange(Id)
  matrix_id = combn(x = index_df$Id, m = multiplexing_level)
  image_1_matrix_binary_sequence = matrix_id_to_2_channel_image_1_binary_sequence(matrix_id = matrix_id, index_df = index_df)#matches Ids to binary sequences
  image_1_logical_rigth_combination = as.logical(x = list_of_good_combinations_2(m = image_1_matrix_binary_sequence))
  image_2_matrix_binary_sequence = matrix_id_to_2_channel_image_2_binary_sequence(matrix_id = matrix_id, index_df = index_df)#matches Ids to binary sequences
  image_2_logical_rigth_combination = as.logical(x = list_of_good_combinations_2(m = image_2_matrix_binary_sequence))
  logical_rigth_combination = as.logical(image_1_logical_rigth_combination*image_2_logical_rigth_combination)
  list_of_all_combinations = matrix_id[, (logical_rigth_combination)] %>% t()
  return(list_of_all_combinations)
}

# get all compatible combinations of an index for 1_channel chemistry
get_all_combinations_1_channel = function(index_df, multiplexing_level){
  index_df = index_binary_conversion_1_channel(index_df)
  index_df = index_df %>% arrange(Id)
  matrix_id = combn(x = index_df$Id, m = multiplexing_level)
  image_1_matrix_binary_sequence = matrix_id_to_1_channel_image_1_binary_sequence(matrix_id = matrix_id, index_df = index_df)#matches Ids to binary sequences
  image_1_logical_rigth_combination = as.logical(x = list_of_good_combinations_2(m = image_1_matrix_binary_sequence))
  image_2_matrix_binary_sequence = matrix_id_to_1_channel_image_2_binary_sequence(matrix_id = matrix_id, index_df = index_df)#matches Ids to binary sequences
  image_2_logical_rigth_combination = as.logical(x = list_of_good_combinations_2(m = image_2_matrix_binary_sequence))
  logical_rigth_combination = as.logical(image_1_logical_rigth_combination*image_2_logical_rigth_combination)
  list_of_all_combinations = matrix_id[, (logical_rigth_combination)] %>% t()
  return(list_of_all_combinations)
}

get_all_combinations = function(index_df, multiplexing_level, chemistry){
  if (chemistry == 4 ){
    combinations_m = get_all_combinations_4_channel(index_df, multiplexing_level)
  } else if (chemistry == 2){
    combinations_m = get_all_combinations_2_channel(index_df, multiplexing_level)
  } else if (chemistry == 1){
    combinations_m = get_all_combinations_1_channel(index_df, multiplexing_level)
  } else {
    display_message("Please choose a correct chemistry for your experiment ")
  }
  return (combinations_m)
}

# get_all_combinations_2 = function(index, multiplexing_level){
#   index = index %>% arrange(Id)
#   matrix_id = combn(x = index$Id, m = multiplexing_level)
#   green_matrix_binary_sequence = green_matrix_id_to_binary_sequence(matrix_id = matrix_id, index = index)#matches Ids to binary sequences
#   green_logical_rigth_combination = as.logical(x = list_of_good_combinations_2(m = green_matrix_binary_sequence))
#   red_matrix_binary_sequence = red_matrix_id_to_binary_sequence(matrix_id = matrix_id, index = index)#matches Ids to binary sequences
#   red_logical_rigth_combination = as.logical(x = list_of_good_combinations_2(m = red_matrix_binary_sequence))
#   logical_rigth_combination = as.logical(green_logical_rigth_combination*red_logical_rigth_combination)
#   list_of_all_combinations = matrix_id[, (logical_rigth_combination)] %>% t()
#   return(list_of_all_combinations)
# }


# For a random search
get_random_combinations_4_channel = function (index_df, multiplexing_level){
  index_df = index_binary_conversion_4_channel(index_df)
  list_of_good_combs = matrix(nrow = 1000, ncol = multiplexing_level)
  j = 0
  for (i in 1:1000){
    combination =index_df[sample(x = 1:nrow(index_df), size = (multiplexing_level), replace = FALSE),]
    if(is_a_good_combination(combination$binary_4)){
      j = j+1
      combination = arrange(combination, Id) #facilitates the distance filter
      list_of_good_combs [j,] = combination$Id
    }
  }
  M = list_of_good_combs %>% unique() %>% na.omit()
  M = M[order(M[,1]),]#facilitates the distance filter
  return (M)
}

get_random_combinations_2_channel = function (index_df, multiplexing_level){
  index_df = index_binary_conversion_2_channel(index_df)
  list_of_good_combs = matrix(nrow = 1000, ncol = multiplexing_level)
  j = 0
  for (i in 1:1000){
    combination =index_df[sample(x = 1:nrow(index_df), size = (multiplexing_level), replace = FALSE),]
    if(is_a_good_combination_2(combination$binary_2_image_1)&&is_a_good_combination_2(combination$binary_2_image_2)){
      j = j+1
      combination = arrange(combination, Id) #facilitates the distance filter
      list_of_good_combs [j,] = combination$Id
    }
  }
  M = list_of_good_combs %>% unique() %>% na.omit()
  M = M[order(M[,1]),]#facilitates the distance filter
  return (M)
}

get_random_combinations_1_channel = function (index_df, multiplexing_level){
  index_df = index_binary_conversion_1_channel(index_df)
  list_of_good_combs = matrix(nrow = 1000, ncol = multiplexing_level)
  j = 0
  for (i in 1:1000){
    combination =index_df[sample(x = 1:nrow(index_df), size = (multiplexing_level), replace = FALSE),]
    if(is_a_good_combination_2(combination$binary_1_image_1)&&is_a_good_combination_2(combination$binary_1_image_2)){
      j = j+1
      combination = arrange(combination, Id) #facilitates the distance filter
      list_of_good_combs [j,] = combination$Id
    }
  }
  M = list_of_good_combs %>% unique() %>% na.omit()
  M = M[order(M[,1]),]#facilitates the distance filter
  return (M)
}

get_random_combinations = function(index_df, multiplexing_level, chemistry){
  if (chemistry == 4 ){
    combinations_m = get_random_combinations_4_channel(index_df, multiplexing_level)
  } else if (chemistry == 2){
    combinations_m = get_random_combinations_2_channel(index_df, multiplexing_level)
  } else if (chemistry == 1){
    combinations_m = get_random_combinations_1_channel(index_df, multiplexing_level)
  } else {
    display_message("Please choose a correct chemistry for your experiment ")
  }
  return (combinations_m)
}


# gets the rights combinations according to the number of possible combinations
get_combinations = function (index, multiplexing_level, chemistry){

  if (choose(nrow(index),multiplexing_level) <= 2024){
    return(get_all_combinations(index, multiplexing_level, chemistry))
  }else {
    return(get_random_combinations(index, multiplexing_level, chemistry))
  }
}





# Filtering  --------------------------------------------------------------

# generates all possible couples and caculates DNAbarcodes's hamming and seqlev distances
index_distance = function (index_df){
  index_distance_df = combn(index_df$sequence,2)  %>% t()%>%  as.data.frame(., stringsAsFactors = FALSE)
  index_distance_df = index_distance_df %>% mutate(n = 1 : nrow(index_distance_df))
  index_distance_df = index_distance_df %>% group_by(n) %>% mutate (hamming = distance(V1, V2, metric = "hamming"),
                                                          seqlev = distance (V1, V2, metric = "seqlev"))
  return(index_distance_df)
}


# couples with hamming distance under threshold
low_hamming_distance = function(index_df, index_distance_df, d){
  i_d = index_distance_df %>% filter(hamming < d)
  i_d = i_d %>% group_by(n) %>% mutate(Id1 = index_df$Id[which(V1 == index_df$sequence)],
                                       Id2 = index_df$Id[which(V2 == index_df$sequence)])# matches the sequence to the id
  low_distance_tab = i_d %>% ungroup() %>% select(Id1, Id2) 
  return(low_distance_tab)
  
}


#couples with seqlev distance under threshold
low_seqlev_distance = function(index_df,index_distance_df, d){
  i_d = index_distance_df %>% filter(seqlev > d) # suoerieur ou superieur ou egal
  i_d = i_d %>% group_by(n) %>% mutate(Id1 = index_df$Id[which(V1 == index_df$sequence)],
                                       Id2 = index_df$Id[which(V2 == index_df$sequence)]) # matches the sequence to the id
  low_distance_tab = i_d %>% ungroup() %>% select(Id1, Id2)
  return(low_distance_tab)
}

# low distance tab = hamming rejection table or seq lev rejection table
filter_combinations = function(combinations_m, low_distance_tab){
  #browser()
  id1 = combinations_m %>% as.data.frame() %>% select_all() == unlist(low_distance_tab$Id1)
  id1 = apply(id1, 1, any)
  id2 = combinations_m %>% as.data.frame() %>% select_all() == unlist(low_distance_tab$Id2)
  id2 = apply(id2, 1, any)
  to_remove = id1 * id2
  return(combinations_m [!as.logical(to_remove),])
  
}

distance_constraints_filter = function(index_df, combinations_m, metric, d){
  index_distance_df =  index_distance(index_df)
  if (metric == "hamming"){
    hamming_rejection_table = low_hamming_distance(index_df, index_distance_df, d)
    filtered_combinations_m = filter_combinations(combinations_m, hamming_rejection_table)
  }else if (metric == "seqlev") {
    seqlev_rejection_table = low_seqlev_distance(index_df, index_distance_df, d) 
    filtered_combinations_m = filter_combinations(combinations_m, seqlev_rejection_table)
  }
  return (filtered_combinations_m)
}





# Result ------------------------------------------------------------------


# Shannon's entropy
shannon_entropy = function(frequence){
  return(-1 * sum(frequence*log(frequence)))
}

#for a matrix of combination
entropy_result = function (index_combination){
  d = index_combination  %>% table()
  d = d/sum(d)
  entropy = shannon_entropy(d)
  return (entropy)
}

# Celine's entropy for given parameters
celine_entropy = function (index_number,sample_number){
  k = index_number
  n = sample_number
  entropy =
    - (k - (n %% k)) * (floor(n/k)/n) * log(floor(n/k)/n) -
    (n %% k) * (ceiling(n/k)/n) * log(ceiling(n/k)/n)
  return(entropy)
}


entropy_max = function (index_number,sample_number){
  if(sample_number > index_number){
    return (celine_entropy(index_number, sample_number))
  }
  else {
    return(log(sample_number))
  }
}

comb_result = function (index_combination, sample_number, index_number, nb_lane){
  # maximal entropy determination for the experiment
  max = entropy_max(index_number,sample_number)
  print(max)
  if(nb_lane < nrow(index_combination)){
    # best combination determination
    one_combination = index_combination [sample(x = 1:nrow(index_combination),nb_lane),]
    an_entropy = entropy_result(one_combination)
    if (an_entropy == max){
      return(one_combination)# we stop if we reach the maximal value of entropy
    }
    else{
      for (i in 1 : 5000){# to debate
        another_combination = index_combination [sample(x = 1:nrow(index_combination),nb_lane),]
        another_entropy = entropy_result(another_combination)
        if(another_entropy == max){
          one_combination = another_combination
          i = 5000 # we stop if we reach the maximal value of entropy
        }
        else if(another_entropy > an_entropy){
          one_combination = another_combination # we keep the best one
          an_entropy = another_entropy # we keep the higher one
        }
      }
    }
  }
  else if (nb_lane >=  nrow(index_combination)){#calcul entopie ????
    {
      one_combination = index_combination[sample(x = 1:nrow(index_combination), nrow(index_combination)),]
      for(i in 2:(floor(nb_lane/nrow(index_combination)))){
        one_combination = rbind(one_combination,
                                index_combination[sample(x = 1:nrow(index_combination), nrow(index_combination)),])
      }
      one_combination = rbind(one_combination,
                              index_combination[sample(x = 1:nrow(index_combination), nb_lane %% nrow(index_combination)),] )
    }
  }
  print(entropy_result(one_combination))
  return (one_combination)
}


# gets the result
get_result = function (index_df, sample_number, multiplexing_level, chemistry, metric = NULL, d = 3){
  combinations_m = get_combinations(index_df, multiplexing_level, chemistry)
  if(!is.null(metric)){
    combinations_m = distance_constraints_filter (index_df, combinations_m, metric, d)
  }
  nb_lane = sample_number %>% as.numeric() / multiplexing_level %>% as.numeric()
  cb = comb_result(combinations_m, sample_number, index_number, nb_lane) %>% as.data.frame()
  result = data.frame(Id = as.vector(cb%>% t() %>% as.vector),
                      Lane = (rep(1:nb_lane, length.out = sample_number, each = multiplexing_level)))
  result$Id = as.character(result$Id)
  return(result)
}


# Experiment Design (single or dual) -----------------------------------------------------------

# change the position of index if there is any duplicate

# Illumina 770-2017-004-C|page 3
# Using unique dual index combinations is a
# best practice to make sure that reads with incorrect indexes do not
# impact variant calling or assignment of gene expression counts.
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

#final function
experiment_design = function (index1,
                              sample_number,
                              multiplexing_level,
                              chemistry = 4,
                              index2 = NULL,
                              export = NULL, 
                              metric = NULL, 
                              d = 3){
  if (is.null(index2)){
    index1  = file_loading_and_checking(index1)
    if (!is.null(index1)) {
      if (sample_and_multiplexing_level_check(sample_number, multiplexing_level) == TRUE){
        print("mlx and sample ok")
        result1 = final_result(index1,sample_number, multiplexing_level,system, metric, d)
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
          result1 = get_result(index1, sample_number, multiplexing_level,metric, d)
          result2 = get_result(index2, sample_number, multiplexing_level,metric, d)
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

display_message = function (a_message){
  error_messsage <<- a_message
  print(a_message)
}



# For java ----------------------------------------------------------------


is_a_prime_number = function (sample_number){
  result = isPrime(sample_number) %>% as.numeric()
  return(result)
}

final_result = function(index1, sample_number, multiplexing_level,system,filter,metric){
  result1 = get_result(index1,sample_number, multiplexing_level,system, filter,metric)
  result1 = data.frame(sample = 1: sample_number %>% as.character(),
                       Lane = result1$Lane %>% as.character(),
                       Id = result1$Id %>% as.character(),
                       stringsAsFactors = FALSE)
  result1 = left_join(result1, select(index1, Id, sequence),by="Id") 
  return (result1)
}


final_result_dual = function(index1,index2, sample_number, multiplexing_level, system,filter,metric){
  result1 = get_result(index1, sample_number, multiplexing_level,system,filter,metric)
  result2 = get_result(index2, sample_number, multiplexing_level,system,filter,metric)
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
  return(result)

}

recursive_entropy = function(combination_m, nb_lane){
 #browser()
  #k = 0
  while(nrow(combination_m)> nb_lane){
    ind = combn(nrow(combination_m), nrow(combination_m)-1)
    a = vector(length = nrow(combination_m))
    #  a=entropy_result(index_combination[ind[,],])
    for(i in 1:nrow(combination_m)){
      a[i]=entropy_result(combination_m[ind[,i],])
      #print(i)
      #b = index_combination[ind[,i],]
      #a = 2
      #k = k + 1
    }
    x = which.max(a)
    z = which(a == a[x])
   # print(z)
    #print(" ")
    #sample(1:length(z),1)
    combination_m = combination_m[ind[,z[sample(1:length(z),1)]],] 
  }
  print(entropy_result(combination_m))
 return (combination_m)
}

#g %>% filter(n == s[,1]) %>% select(-n) %>% unlist() %>% entropy_result()

exchange_entropy = function (combination_m, nb_lane){
  comb_m = sample
}

optimized_combination = function (combination_m, nb_lane, index_number){
  max = entropy_max(index_number, length(combination_m) * nb_lane)
  if(nb_lane < nrow(combination_m)){
    if (nrow(combination_m >80)){
      a_combination = combination_m[sample(1:nrow(combination_m),80)]
      a_combination = recursive_entropy(a_combination,nb_lane)
      
      while (i < 10 || entropy_result(a_combination) < entropy_max)
        temp_combination = recursive_entropy(combination_m[sample(1:nrow(combination_m),80)], nb_lane)
      if (entropy_result(temp_combination) > entropy_result(a_combination) ){
        a_combination = temp_combination
      }
      
    }else { 
      a_combination = recursive_entropy(combination_m,nb_lane)
      while (i < 10 || entropy_result(combination_m) < entropy_max)
        temp_combination = recursive_entropy(combination_m,nb_lane)
      if (entropy_result(temp_combination) > entropy_result(a_combination) ){
        a_combination = temp_combination
      }
    }
    
  }
}

    if(combination_m > 100){
      possible_combination = 
      for (i in 1:10){
      sample_combination = combination_m[sample(1:nrow(combination_m),100)]
      one_combination = recursive_entropy()
      }
    }
    
  }
}
one_combination = index_combination [sample(x = 1:nrow(index_combination),nb_lane),]
an_entropy = entropy_result(one_combination)
if (an_entropy == max){
  return(one_combination)# we stop if we reach the maximal value of entropy
}
else{
  for (i in 1 : 10000){# to debate
    another_combination = index_combination [sample(x = 1:nrow(index_combination),nb_lane),]
    another_entropy = entropy_result(another_combination)
    if(another_entropy == max){
      one_combination = another_combination
      i = 10000 # we stop if we reach the maximal value of entropy
    }
    else if(another_entropy > an_entropy){
      one_combination = another_combination # we keep the best one
      an_entropy = another_entropy # we keep the higher one
    }




q = t[sample(1:nrow(t),100),]
entropy_result(q)
n = 10
recursive_entropy(q,n)
recursive_entropy(q,n)
recursive_entropy(q,n)
recursive_entropy(q,n)
recursive_entropy(q,n)