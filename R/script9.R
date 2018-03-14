
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




file_loading_and_checking = function(file){
  index = read_index(file) 
  if (!is.null(index) && index_check(index)){#  if no issue
    index_number <<- nrow(index)
    index  = index %>% mutate (GC_content = get_index_GC_content(index), 
                               homopolymer = get_index_homopolymer(index))
    index_distances <<- index_distance(index)
    return (index)
  }else{
    return(NULL)
  }
}


# generates all possible couples and caculates DNAbarcodes's hamming and seqlev distances
index_distance = function (index){
  index_couple = combn(index$sequence,2)  %>% t()%>%  as.data.frame(., stringsAsFactors = FALSE)
  index_couple = index_couple %>% mutate(n = 1 : nrow(index_couple))
  index_couple = index_couple %>% group_by(n) %>% mutate (hamming = DNABarcodes::distance(V1, V2, metric = "hamming"),
                                                          seqlev = DNABarcodes::distance (V1, V2, metric = "seqlev"))
  return(index_couple)
}

# couples with hamming distance under threshold
low_hamming_distance = function(index_distances, metric){
  
  i_d = index_distances %>% filter(hamming < metric)
  i_d = i_d %>% group_by(n) %>% mutate(Id1 = index$Id[which(V1 == index$sequence)],
                                       Id2 = index$Id[which(V2 == index$sequence)])
  lhd = i_d %>% ungroup() %>% select(Id1, Id2) 
  return(lhd)
}


#couples with seqlev distance under threshold
low_seqlev_distance = function(index_distances, metric){
  i_d = index_distances %>% filter(seqlev < metric)
  i_d = i_d %>% group_by(n) %>% mutate(Id1 = index$Id[which(V1 == index$sequence)],
                                       Id2 = index$Id[which(V2 == index$sequence)]) 
  lsd = i_d %>% ungroup() %>% select(Id1, Id2)
  return(lsd)
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
  return (multiplexing_level_choices)
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

# Binary conversion -------------------------------------------------------


#A function is used for the traduction of a sequence :
sequence_binary_conversion = function(sequence){
  binary_sequence = 
    toupper(sequence) %>% 
    gsub("A|C", "0", .) %>% 
    gsub("G|T", "1", .)
  return(binary_sequence)
}


#creates and add to the index the corresponding binary sequence to each row
index_binary_conversion = function(index){
  index = index %>% mutate (binary = sequence_binary_conversion(sequence))
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
  }
  else{
    return (FALSE)
  }
}

#check if a combination of indexes is correct
is_a_good_combination = function (combination_matrix){
  all_combinations = vectors_into_matrix(combination_matrix)
  results = prod(apply(all_combinations,1,any_different))
  return(results)
}

#keeps only the good ones :
list_of_good_combinations = function (matrix_id){
  list = apply(matrix_id,2,is_a_good_combination)
  return(list)
}

##super fast and furious
#matches an id to its binary_sequence
id_into_binary_sequence = function (index_id_combination, index){
  index_rows = subset(index, Id == index_id_combination)
  index_binary_sequence_combination = as.character(index_rows$binary)
  return (index_binary_sequence_combination)
}

#matches an id to its binary_sequence
id_into_binary_sequence_slow = function (an_index_id_combination, index){
  result = as.matrix(index %>% filter (Id %in% an_index_id_combination) %>% select(binary))
  return (result)
}

#matches the matrix_id to the matrix_binary_sequence
matrix_id_to_binary_sequence = function(matrix_id, index){
  m = matrix_id
  n = nrow(m)
  m[1:n,] = sapply(X = m[1:n,], FUN = id_into_binary_sequence, index)
  return (m)
}

# get all compatible combinations of an index
get_all_combinations = function(index, multiplexing_level){
  index = index %>% arrange(Id)
  matrix_id = combn(x = index$Id, m = multiplexing_level)
  matrix_binary_sequence = matrix_id_to_binary_sequence(matrix_id = matrix_id, index = index)#matches Ids to binary sequences
  logical_rigth_combination = as.logical(x = list_of_good_combinations(m = matrix_binary_sequence))
  list_of_all_combinations = matrix_id[, logical_rigth_combination] %>% t()
  return(list_of_all_combinations)
}


# For a random search
get_random_combinations = function (index, multiplexing_level){
  list_of_good_combs = matrix(nrow = 700, ncol = multiplexing_level)
  j = 0
  for (i in 1:700){
    combination =index[sample(x = 1:nrow(index), size = (multiplexing_level), replace = FALSE),]
    if(is_a_good_combination(combination$binary)){
      j = j+1
      combination = arrange(combination, Id) #facilitates the distance filter
      list_of_good_combs [j,] = combination$Id
    }
  }
  M = list_of_good_combs %>% unique() %>% na.omit()
  M = M[order(M[,1]),]#facilitates the distance filter
  return (M)
}




# gets the rights combinations according to the number of possible combinations
get_combinations = function (index, multiplexing_level){
  if (choose(nrow(index),multiplexing_level) <= 2024){
    return(get_all_combinations(index, multiplexing_level))
  }else {
    return(get_random_combinations(index, multiplexing_level))
  }
}
# low distance tab = hamming rejection table or seq lev rejection table
filter_combinations = function(combinations, low_distance_tab){
  
  id1 = combinations %>% as.data.frame() %>% select_all() == unlist(low_distance_tab$Id1)
  id1 = apply(id1, 1, any)
  id2 = combinations%>% as.data.frame() %>% select_all() == unlist(low_distance_tab$Id2)
  id2 = apply(id2, 1, any)
  
  to_remove = id1 * id2
  return(combinations [!as.logical(to_remove),])
  
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


# Max entropy
entropy_max = function (index_number,sample_number){
  if(sample_number > index_number){
    k = index_number
    n = sample_number
    entropy =
      - (k - (n %% k)) * (floor(n/k)/n) * log(floor(n/k)/n) -
      (n %% k) * (ceiling(n/k)/n) * log(ceiling(n/k)/n)
    return (entropy)
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
      for (i in 1 : 1000){# to debate
        another_combination = index_combination [sample(x = 1:nrow(index_combination),nb_lane),]
        another_entropy = entropy_result(another_combination)
        if(another_entropy == max){
          one_combination = another_combination
          i = 1000 # we stop if we reach the maximal value of entropy
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
get_result = function (index, sample_number, multiplexing_level, filter = NULL, metric = 3){
  index = index %>% index_binary_conversion()
  index_combinations = get_combinations(index, multiplexing_level)
  if(!is.null(filter)){
    if(filter == "hamming"){
      hamming_rejection_table = low_hamming_distance(index_distances, metric) 
      index_combinations <<- filter_combinations(index_combinations, hamming_rejection_table)
    }else {
      seqlev_rejection_table = low_seqlev_distance(index_distances, metric) 
      index_combinations <<- filter_combinations(index_combinations, seqlev_rejection_table)
    }
  }
  nb_lane = sample_number %>% as.numeric() / multiplexing_level %>% as.numeric()
  cb = comb_result(index_combinations, sample_number, index_number, nb_lane) %>% as.data.frame()
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

display_message = function (a_message){
  error_messsage <<- a_message
  print(a_message)
}

# For java ----------------------------------------------------------------


is_a_prime_number = function (sample_number){
  result = isPrime(sample_number) %>% as.numeric()
  return(result)
}





