

library("dplyr")
library("tidyr")
library("numbers")
library("purrr")
library ("stringr")
library("DNABarcodes")


# initialize variable for CRAN/Bioc repository compatibility
globalVariables(
    c(
        "Id",
        "Id1",
        "Id2",
        "Lane",
        "V1",
        "V2",
        "barcodes",
        "bcID",
        "bc_gp",
        "combID",
        "hamming",
        "index",
        "isDropping",
        "matchOcc",
        "pairID",
        "seqlev",
        "phaseshift",
        "."
    )
)


# Inputs ------------------------------------------------------------------
index_env <- new.env() # ?bindenv() for help
read_index = function(file) {
    if (!file.exists(file)) {
        display_message("Your file doesn't exist, please check the path")
        return(NULL)
        assign("index", NULL, envir = index_env) #index <<- NULL
    } else{
        input <- try(as.data.frame(
            read.table(
                file,
                header = FALSE,
                sep = "",
                col.names = c("Id", "sequence"),
                colClasses = c("character", "character"),
                stringsAsFactors = FALSE
            )
        ), silent = TRUE)
        if (exists("input")) {
            if (! "try-error" %in% is(input)) {
                assign("index", input, envir = index_env)  
                return(input)
            } else {
                assign("index", NULL, envir = index_env) #index <<- NULL
                display_message(
"An error occurred, please check the content of your file")
                return(NULL)

            }
        }
    }
}


unicity_check = function(index) {
    index$sequence = toupper(index$sequence)
    assign("index_df", index, envir = index_env)
    if (index$Id %>% anyDuplicated() != 0) {
        #checks if the index Ids are unique
        v = paste("two duplicated indexes IDs,check row number",
                anyDuplicated(index$Id))
        display_message(v)
        return(FALSE)
    }
    else if (index$sequence %>% anyDuplicated() != 0) {
        #checks if the index sequences are unique
        v = paste(
            "At least one of your sequences is not unique, check row number",
            anyDuplicated(index$sequence)
        )
        display_message(v)
        return(FALSE)
    }
    else{
        return (TRUE)
    }
}


#check the character for one sequence
character_check = function(sequence) {
    wrong_letters = LETTERS[!(LETTERS) %in% c("A", "G", "C", "T")]
    check = all(!str_detect(sequence, wrong_letters))
    return(check)
}


#checks the characters of each sequence
sequences_character_check = function(sequences) {
    return(map_lgl(sequences, character_check))
}


#checks the length of each sequence
length_check = function(sequences) {
    return(str_length(sequences) == str_length(sequences[1]))
}


#checks for character type and sequence  length
character_and_length_check = function(index) {
    c_check =  sequences_character_check(index$sequence)
    if (length (which(c_check == FALSE)) > 0) {
        display_message(paste(
            "Your sequence contains a wrong character, row number : ",
            which(c_check == FALSE)
        ))
        return(FALSE)
    }
    else {
        l_check = length_check(index$sequence)
        if (length(which(l_check == FALSE)) > 0) {
            wrong_length = as.numeric(
                names(table(nchar(index$sequence)))[
                    as.numeric(
                    which(table(nchar(index$sequence)) == 
                            min(table(nchar(index$sequence)))))
                    ])
            c = which(nchar(index$sequence) == wrong_length)
            display_message(paste("the indexes have not the same size,
                            check out row number", c))
            
            return(FALSE)
        } else{
            return(TRUE)
        }
    }
}



## Check the index for possible issues
index_check = function(index) {
    return (all(unicity_check(index),
                character_and_length_check(index)))
}


get_sequence_GC_content =  function (sequence) {
    GC_content = str_count( sequence,
                            pattern = ("G|C")) / nchar(sequence) * 100
    return(round(GC_content, digits = 2))
}


get_index_GC_content = function(index) {
    index_GC_content = map_dbl(index$sequence, get_sequence_GC_content)
    return (index_GC_content)
}


get_sequence_homopolymer = function(sequence) {
    if (length(str_which(sequence, "A{3,}|G{3,}|C{3,}|T{3,}")) > 0) {
        return (TRUE)
    } else{
        return (FALSE)
    }
}


get_index_homopolymer = function(index) {
    homopolymer  = map_lgl(index$sequence, get_sequence_homopolymer)
}



sample_number_check = function (sample_number) {
    if (is.na(sample_number)) {
        display_message("you have to enter an integer value")
        return(FALSE)
    }
    else if (sample_number != floor(sample_number)) {
        display_message("you have to enter an integer value")
        return(FALSE)
    }
    else if (sample_number < 2) {
        display_message(paste("You need at least 2 samples in order to",
                        "use multiplexing"))
        return(FALSE)
    } else if (sample_number > 1000) {
        display_message(paste("The sample number is too high,",
                        "please enter a value under 1000"))
        return(FALSE)
    }
    else {
        return(TRUE)
    }
}



multiplexing_level_set = function (sample_number) {
    v  = 2:(sample_number)
    multiplexing_level_choices = v[sample_number %% v == 0]
    return (multiplexing_level_choices [multiplexing_level_choices < 96])
}


sample_and_multiplexing_level_check = function(sample_number, mplex_level) {
    if (sample_number_check(sample_number)) {
        possible_multiplexing_level = multiplexing_level_set(sample_number)
        if (!(mplex_level %in% possible_multiplexing_level)) {
            display_message(paste(
                "The sample number isn't a multiple of the multiplexing level.",
        "Here are the possible multiplexing levels :"
            ))
            display_message(possible_multiplexing_level)
            return(FALSE)
        } else {
            return (TRUE)
        }
    } else {
        return(FALSE)
    }
}




# sequence traduction for a 4_channel platform :
sequence_binary_conversion_4_channel = function(sequence) {
    binary_sequence =
        toupper(sequence) %>%
        gsub("A|C", "0", .) %>%
        gsub("G|T", "1", .)
    return(binary_sequence)
}



# index traduction for a 4_channel platform :
index_binary_conversion_4_channel = function(index) {
    index = index %>% 
        mutate (binary_4 = sequence_binary_conversion_4_channel(sequence))
    return (index)
}



# sequence traduction for a 2_channel platform for image 1 :
sequence_binary_conversion_2_channel_1 = function(sequence) {
    binary_sequence =
        toupper(sequence) %>%
        gsub("G|C", "0", .) %>%
        gsub("A|T", "1", .)
    return(binary_sequence)
}



# sequence traduction for a 2_channel platform for image 2 :
sequence_binary_conversion_2_channel_2 = function(sequence) {
    binary_sequence =
        toupper(sequence) %>%
        gsub("G|T", "0", .) %>%
        gsub("A|C", "1", .)
    return(binary_sequence)
}



# index traduction for a 2_channel platform :
index_binary_conversion_2_channel = function(index) {
    index = index %>% 
        mutate (
        binary_2_image_1 = sequence_binary_conversion_2_channel_1(sequence),
        binary_2_image_2 = sequence_binary_conversion_2_channel_2(sequence)
        )
    return (index)
}



# sequence traduction for a 1_channel platform for image 1 :
sequence_binary_conversion_1_channel_1 = function(sequence) {
    binary_sequence =
        toupper(sequence) %>%
        gsub("G|C", "0", .) %>%
        gsub("A|T", "1", .)
    return(binary_sequence)
}



# sequence traduction for a 1_channel platform for image 2 :
sequence_binary_conversion_1_channel_2 = function(sequence) {
    binary_sequence =
        toupper(sequence) %>%
        gsub("A|G", "0", .) %>%
        gsub("C|T", "1", .)
    return(binary_sequence)
}



index_binary_conversion_1_channel = function(index) {
    index = index %>%
        mutate (
        binary_1_image_1 = sequence_binary_conversion_1_channel_1(sequence),
        binary_1_image_2 = sequence_binary_conversion_1_channel_2(sequence)
        )
    return (index)
}


# Compatibility -----------------------------------------------------------




#conversion of the string sequence into a vector of numeric
binary_word_into_numeric = function (binary_word) {
    as.numeric(unlist(strsplit(as.character(binary_word), "")))
}



#conversion of each index sequence combination into a matrix
vectors_into_matrix = function (binary_word) {
    m =  mapply(binary_word_into_numeric, binary_word)
    return(m)
}



#test if a column/line of a index combination is correct
any_different = function(binary_sequence) {
    if (length(unique(binary_sequence)) > 1) {
        return (TRUE)
    } else{
        return (FALSE)
    }
}



has_signal_in_both_channels = function(colored_sequence) {
    if (any(as.logical(colored_sequence))) {
        return (TRUE)
    } else{
        return (FALSE)
    }
}



# check if a combination of indexes is correct
is_a_good_combination = function (combination_matrix) {
    all_combinations = vectors_into_matrix(combination_matrix)
    results = prod(apply(all_combinations, 1, any_different))
    return(results)
}



is_a_good_combination_2 = function (combination_matrix) {
    all_combinations = vectors_into_matrix(combination_matrix)
    results = prod(apply(all_combinations, 1, has_signal_in_both_channels))
    return(results)
}



# keeps only the good ones :
list_of_good_combinations = function (matrix_id) {
    list = apply(matrix_id, 2, is_a_good_combination)
    return(list)
}



list_of_good_combinations_2 = function (matrix_id) {
    list = apply(matrix_id, 2, is_a_good_combination_2)
    return(list)
}



##super fast and furious
#matches an id to its binary_sequence
id_into_4_channel_binary_sequence = function (index_id_combination, index_df) {
    index_rows = subset(index_df, Id == index_id_combination)
    index_binary_sequence_combination = as.character(index_rows$binary_4)
    return (index_binary_sequence_combination)
}



# matches an id to its binary_sequence for 2_channel image 1
id_into_2_channel_image_1_binary_sequence = function (index_id_combination,
                                                    index_df) {
    index_rows = subset(index_df, Id == index_id_combination)
    index_binary_sequence_combination = 
        as.character(index_rows$binary_2_image_1)
    return (index_binary_sequence_combination)
}



# matches an id to its binary_sequence for 2_channel image 2
id_into_2_channel_image_2_binary_sequence = function (index_id_combination,
                                                    index_df) {
    index_rows = subset(index_df, Id == index_id_combination)
    index_binary_sequence_combination = 
        as.character(index_rows$binary_2_image_2)
    return (index_binary_sequence_combination)
}



# matches an id to its binary_sequence for 1_channel image 1
id_into_1_channel_image_1_binary_sequence = function (index_id_combination,
                                                    index_df) {
    index_rows = subset(index_df, Id == index_id_combination)
    index_binary_sequence_combination = 
        as.character(index_rows$binary_1_image_1)
    return (index_binary_sequence_combination)
}



# matches an id to its binary_sequence for 1_channel image 2
id_into_1_channel_image_2_binary_sequence = function (index_id_combination,
                                                    index_df) {
    index_rows = subset(index_df, Id == index_id_combination)
    index_binary_sequence_combination = 
        as.character(index_rows$binary_1_image_2)
    return (index_binary_sequence_combination)
}





# matches the matrix_id to the matrix_binary_sequence
matrix_id_to_binary_sequence = function(matrix_id, index_df) {
    m = matrix_id
    n = nrow(m)
    # m[1:n, ] = sapply(X = m[1:n, ], FUN = id_into_4_channel_binary_sequence,
    #                 index_df)
    m[seq(1,n), ] = vapply(X = m[seq(1,n), ],
                        FUN = id_into_4_channel_binary_sequence,
                        "",
                        index_df)
    return (m)
}



# matches the matrix_id to the matrix_2_channel_immage_1_binary_sequence
matrix_id_to_2_channel_image_1_binary_sequence = function(matrix_id, index_df){
    m = matrix_id
    n = nrow(m)
    m[seq(1,n), ] = vapply(X = m[seq(1,n), ],
                    FUN = id_into_2_channel_image_1_binary_sequence,
                    "",
                    index_df)
    return (m)
}



# matches the matrix_id to the matrix_2_channel_immage_2_binary_sequence
matrix_id_to_2_channel_image_2_binary_sequence = function(matrix_id, index_df){
    m = matrix_id
    n = nrow(m)
    m[seq(1,n), ] = vapply(X = m[seq(1,n), ],
                    FUN = id_into_2_channel_image_2_binary_sequence,
                    "",
                    index_df)
    return (m)
}



# matches the matrix_id to the matrix_1_channel_immage_1_binary_sequence
matrix_id_to_1_channel_image_1_binary_sequence = function(matrix_id, index_df){
    m = matrix_id
    n = nrow(m)
    m[seq(1,n), ] = vapply(X = m[seq(1,n), ],
                    FUN = id_into_1_channel_image_1_binary_sequence,
                    "",
                    index_df)
    return (m)
}



# matches the matrix_id to the matrix_1_channel_immage_2_binary_sequence
matrix_id_to_1_channel_image_2_binary_sequence = function(matrix_id, index_df){
    m = matrix_id
    n = nrow(m)
    m[seq(1,n), ] = vapply(X = m[seq(1,n), ],
                    FUN = id_into_1_channel_image_2_binary_sequence,
                    "",
                    index_df)
    return (m)
}



# get all compatible combinations of an index for 4_channel platform
get_all_combinations_4_channel = function(index_df, mplex_level) {
    index_df = index_binary_conversion_4_channel(index_df)
    index_df = index_df %>% arrange(Id)
    matrix_id = combn(x = index_df$Id, m = mplex_level)
    matrix_binary_sequence = matrix_id_to_binary_sequence(matrix_id, index_df)
    logical_rigth_combination = as.logical (
        x = list_of_good_combinations(matrix_id = matrix_binary_sequence)
    )
    list_of_all_combinations = matrix_id[, logical_rigth_combination] %>% t()
    return(list_of_all_combinations)
}



# get all compatible combinations of an index for 2_channel platform
get_all_combinations_2_channel = function(index_df, mplex_level) {
    index_df = index_binary_conversion_2_channel(index_df)
    index_df = index_df %>% arrange(Id)
    matrix_id = combn(x = index_df$Id, m = mplex_level)
    #matches Ids to binary sequences
    image_1_matrix_binary_sequence = 
        matrix_id_to_2_channel_image_1_binary_sequence(matrix_id = matrix_id,
                                                        index_df = index_df)
    image_1_logical_rigth_combination = as.logical (
        x = list_of_good_combinations_2(image_1_matrix_binary_sequence)
    )
    #matches Ids to binary sequences
    image_2_matrix_binary_sequence = 
        matrix_id_to_2_channel_image_2_binary_sequence(matrix_id = matrix_id,
                                                        index_df = index_df)
    image_2_logical_rigth_combination = as.logical (
        x = list_of_good_combinations_2(image_2_matrix_binary_sequence)
    )
    
    logical_rigth_combination = as.logical(
        image_1_logical_rigth_combination * image_2_logical_rigth_combination)
    
    list_of_all_combinations = matrix_id[, (logical_rigth_combination)] %>%
        t()
    
    return(list_of_all_combinations)
}



# get all compatible combinations of an index for 1_channel platform
get_all_combinations_1_channel = function(index_df, mplex_level) {
    index_df = index_binary_conversion_1_channel(index_df)
    index_df = index_df %>% arrange(Id)
    matrix_id = combn(x = index_df$Id, m = mplex_level)
    #matches Ids to binary sequences
    image_1_matrix_binary_sequence = 
        matrix_id_to_1_channel_image_1_binary_sequence(matrix_id = matrix_id,
                                                        index_df = index_df)
    image_1_logical_rigth_combination = as.logical(
        x = list_of_good_combinations_2(image_1_matrix_binary_sequence)
    )
    #matches Ids to binary sequences
    image_2_matrix_binary_sequence =
        matrix_id_to_1_channel_image_2_binary_sequence(matrix_id = matrix_id,
                                                        index_df = index_df)
    image_2_logical_rigth_combination = as.logical(
        x = list_of_good_combinations_2(image_2_matrix_binary_sequence)
    )
    
    logical_rigth_combination = as.logical(
        image_1_logical_rigth_combination * image_2_logical_rigth_combination)
    
    list_of_all_combinations = matrix_id[, (logical_rigth_combination)] %>% 
        t()
    
    return(list_of_all_combinations)
}
#get all combination for other platform than Illumina (no constraint)
get_all_combinations_0_channel = function(index_df, mplex_level) {
    list_of_all_combinations  = combn(x = index_df$Id, m = mplex_level) %>% 
        t()
    return(list_of_all_combinations)
}

# For a random search
get_random_combinations_4_channel = function (index_df, mplex_level) {
    index_df = index_binary_conversion_4_channel(index_df)
    list_of_good_combs = matrix(nrow = 1000, ncol = mplex_level)
    j = 0
    for (i in seq(1,1000)) {
        combination = index_df[sample(
            x = seq(1,nrow(index_df)),
            size = (mplex_level),
            replace = FALSE
        ), ]
        if (is_a_good_combination(combination$binary_4)) {
            j = j + 1
            #facilitates the distance filter
            combination = arrange(combination, Id) 
            list_of_good_combs [j, ] = combination$Id
        }
    }
    M = list_of_good_combs %>% unique() %>% na.omit()
    #facilitates the distance filter
    M = M[order(M[, 1]), ]
    return (M)
}



get_random_combinations_2_channel = function (index_df, mplex_level) {
    index_df = index_binary_conversion_2_channel(index_df)
    list_of_good_combs = matrix(nrow = 1000, ncol = mplex_level)
    j = 0
    for (i in seq(1,1000)) {
        combination = index_df[sample(
            x = seq(1,nrow(index_df)),
            size = (mplex_level),
            replace = FALSE
        ), ]
        if (is_a_good_combination_2(combination$binary_2_image_1) &&
            is_a_good_combination_2(combination$binary_2_image_2)) {
            j = j + 1
            #facilitates the distance filter
            combination = arrange(combination, Id)
            list_of_good_combs [j, ] = combination$Id
        }
    }
    M = list_of_good_combs %>% unique() %>% na.omit()
    #facilitates the distance filter
    M = M[order(M[, 1]), ]
    return (M)
}



get_random_combinations_1_channel = function (index_df, mplex_level) {
    index_df = index_binary_conversion_1_channel(index_df)
    list_of_good_combs = matrix(nrow = 1000, ncol = mplex_level)
    j = 0
    for (i in seq(1,1000)) {
        combination = index_df[sample(
            x = seq(1,nrow(index_df)),
            size = (mplex_level),
            replace = FALSE
        ), ]
        if (is_a_good_combination_2(combination$binary_1_image_1) &&
            is_a_good_combination_2(combination$binary_1_image_2)) {
            j = j + 1
            #facilitates the distance filter
            combination = arrange(combination, Id) 
            list_of_good_combs [j, ] = combination$Id
        }
    }
    M = list_of_good_combs %>% unique() %>% na.omit()
    #facilitates the distance filter
    M = M[order(M[, 1]), ]
    return (M)
}

get_random_combinations_0_channel = function (index_df, mplex_level) {
    list_of_good_combs = matrix(nrow = 1000, ncol = mplex_level)
        for (i in seq(1,1000)) {
            combination = index_df[sample(
                x = seq(1,nrow(index_df)),
                size = (mplex_level),
                replace = FALSE
            ), ]
            
            combination = arrange(combination, Id) 
            list_of_good_combs [i, ] = combination$Id
        }
    M = list_of_good_combs %>% unique()
    #facilitates the distance filter
    
    M = M[order(M[, 1]), ]
    return (M)
}



# gets the rights combinations according to the number of possible combinations
get_combinations = function (index_df, mplex_level, platform){
    if (choose(nrow(index_df),mplex_level) <= 2024 || mplex_level == 2){
        return(get_all_combinations(index_df, mplex_level, platform))
    } else {
        return(get_random_combinations(index_df, mplex_level, platform))
    }
}





# Filtering  --------------------------------------------------------------



# generates all possible couples 
# and calculates DNAbarcodes's hamming and seqlev distances
index_distance = function (index_df) {
    index_distance_df = combn(index_df$sequence, 2) %>% t() %>% 
        as.data.frame(stringsAsFactors = FALSE)
    index_distance_df = index_distance_df %>%
        mutate(n = seq(1,nrow(index_distance_df)))
    index_distance_df = index_distance_df %>% group_by(n) %>% mutate (
        hamming = distance(V1, V2, metric = "hamming"),
        seqlev = distance (V1, V2, metric = "seqlev"),
        phaseshift = distance(V1, V2, metric = "phaseshift")
    )
    return(index_distance_df)
}



# couples with hamming distance under threshold
low_hamming_distance = function(index_df, index_distance_df, d) {
    i_d = index_distance_df %>% filter(hamming < d)
    i_d = i_d %>% group_by(n) %>%
        # matches the sequence to the id
        mutate(Id1 = index_df$Id[which(V1 == index_df$sequence)],
                Id2 = index_df$Id[which(V2 == index_df$sequence)])
    low_distance_tab = i_d %>% ungroup() %>% select(Id1, Id2)
    return(low_distance_tab)
    
}



#couples with seqlev distance under threshold
low_seqlev_distance = function(index_df, index_distance_df, d) {
    i_d = index_distance_df %>% filter(seqlev < d)
    i_d = i_d %>% group_by(n) %>% 
        # matches the sequence to the id
        mutate(Id1 = index_df$Id[which(V1 == index_df$sequence)],
                Id2 = index_df$Id[which(V2 == index_df$sequence)]) 
    low_distance_tab = i_d %>% ungroup() %>% select(Id1, Id2)
    return(low_distance_tab)
}

#couples with phaseshift distance under threshold
low_phaseshift_distance = function(index_df, index_distance_df, d) {
    i_d = index_distance_df %>% filter(phaseshift < d)
    i_d = i_d %>% group_by(n) %>% 
        # matches the sequence to the id
        mutate( Id1 = index_df$Id[which(V1 == index_df$sequence)],
                Id2 = index_df$Id[which(V2 == index_df$sequence)]) 
    low_distance_tab = i_d %>% ungroup() %>% select(Id1, Id2)
    return(low_distance_tab)
}



# low distance tab = hamming rejection table or seqlev rejection table
filter_combinations = function(combinations_m, low_distance_tab) {
    # browser()
    combinations_df <-
        combinations_m %>% as.data.frame(., stringsAsFactors = FALSE) %>%
        mutate(combID =
                    0:((nrow(.) > 0) * (nrow(.) - 1)))
    # head(combinations_df)
    
    combinations_df_long <- gather(combinations_df, 
                                    bc_gp, 
                                    barcodes,
                                    -combID) %>%
        select(-bc_gp) %>% arrange(combID)
    # head(combinations_df_long)
    
    dist_tab_long <- gather(
        low_distance_tab %>% mutate(pairID = 0:((nrow(low_distance_tab) > 0)
                                            *(nrow(low_distance_tab) - 1))),
        bcID,
        barcodes,
        -pairID) %>%
        select(-bcID) %>% arrange(pairID)
    
    droppedComb <- combinations_df_long %>%  group_by(combID) %>%
        do ({
            single_comb <- .
            dist_tab_long %>% group_by(pairID) %>%
                summarise(
                    matchOcc = sum(barcodes %in% single_comb$barcodes)
                ) %>%
                ungroup() %>%
                mutate(isDropping = as.logical(matchOcc == 2)) %>%
                summarise(isDropping = any(isDropping))
        }) %>% ungroup()
    
    
    keptComb <-
        inner_join(combinations_df, droppedComb, by = "combID") %>%
        filter(!isDropping) %>% select(-c(combID, isDropping)) %>% as.matrix()
    
    # if (nrow(keptComb)<1) stop("No combination left after filtering,
    # please reduce the threshold.")
    
    return(keptComb)
    
}


# Result ------------------------------------------------------------------



# Shannon's entropy
shannon_entropy = function(frequency) {
    return(-1 * sum(frequency * log(frequency)))
}



#for a matrix of combination
entropy_result = function (index_combination) {
    d = index_combination  %>% table()
    d = d / sum(d)
    entropy = shannon_entropy(d)
    return (entropy)
}



# Celine's entropy for given parameters
entropy_n_k = function (index_number, sample_number) {
    k = index_number
    n = sample_number
    entropy =
        -(k - (n %% k)) * (floor(n / k) / n) * log(floor(n / k) / n) -
        (n %% k) * (ceiling(n / k) / n) * log(ceiling(n / k) / n)
    return(entropy)
}



entropy_max = function (index_number, sample_number) {
    if (sample_number > index_number) {
        return (entropy_n_k(index_number, sample_number))
    }
    else {
        return(log(sample_number))
    }
}



recursive_entropy = function(combination_m,
                            nb_lane,
                            method = "greedy_exchange") {
    
    if (!method %in% c("greedy_exchange", "greedy_descent")) {
        display_message(
            paste0(
                "recursive_entropy(): the selected method '",
                method,
                "' does not exist. Please select 'greedy_exchange' 
                or 'greedy_descent'"
            )
        )
        return(combination_m)
    }
    
    if (method == "greedy_descent") {
        # Celine's implementation
        while (nrow(combination_m) > nb_lane) {
            ind = combn(nrow(combination_m), nrow(combination_m) - 1)
            a = vector(length = nrow(combination_m))
            for (i in seq(1,nrow(combination_m))) {
                a[i] = entropy_result(combination_m[ind[, i], ])
            }
            x = which.max(a)
            z = which(a == a[x])
        combination_m = combination_m[ind[, z[sample(seq(1,length(z)), 1)]], ]
        }
        return (combination_m)
    }
    
    if (method == "greedy_exchange") {
        # Jacques' implementation from his Matlab code
        clist = combination_m
        ncomb = nb_lane
        c0 = combination_m[sample(seq(1,nrow(combination_m)), nb_lane), ]
        N = nrow(combination_m)
        c = c0
        ctest = c
        Sc = entropy_result(c)
        test = 1
        
        while (test) {
            test = 0
            Stest = Sc
            
            n = 1
            i = 1
            
            while ((n <= ncomb) && (i <= N) && (test == 0)) {
                ctp = c
                
                ctp[n, ] = clist[i, ]
                
                Stp = entropy_result(ctp)
                
                if (Stp > Sc) {
                    test = 1
                    Stest = Stp
                    ctest = ctp
                }
                if (i < N) {
                    i = i + 1
                } else {
                    i = 1
                    n = n + 1
                }
            }
            if (Stest > Sc) {
                c = ctest
                Sc = Stest
            }
        }
        return (ctest)
    }
}


# gets the result
get_result = function ( index_df,
                        sample_number,
                        mplex_level,
                        platform,
                        metric = NULL,
                        d = 3,
                        thrs_size_comb = 120,
                        max_iteration = 50,
                        method = "greedy_exchange") {
    
    combinations_m = get_combinations(index_df, mplex_level, platform)
    if (!is.null(metric)) {
        combinations_m = distance_filter (index_df, combinations_m, metric, d)
    }
    if(!is.null(combinations_m)){
        if (sample_number ==  mplex_level){
            result = combinations_m[sample(seq(1,nrow(combinations_m)), 1),]
            result = data.frame(Id = result,
                                Lane = (rep(
                                    seq(1,1),
                                    length.out = sample_number,
                                    each = mplex_level
                                )))
            result$Id = as.character(result$Id)
            return (result)
        }else{
            nb_lane = sample_number %>% 
                as.numeric() / mplex_level %>% as.numeric()
            index_number = nrow(index_df)
            cb = optimize_combinations(combinations_m,
                                        nb_lane,
                                        index_number,
                                        thrs_size_comb,
                                        max_iteration,
                                        method) %>% as.data.frame()
            result = data.frame(Id = as.vector(cb %>% t() %>% as.vector),
                                Lane = (rep(
                                    seq(1,nb_lane),
                                    length.out = sample_number,
                                    each = mplex_level
                                )))
            result$Id = as.character(result$Id)
            return(result)}
    } else {
        return(NULL)
    }
}



# Experiment Design (single or dual) -----------------------------------------

# change the position of index if there is any duplicate


# Illumina 770-2017-004-C|page 3
# Using unique dual index combinations is a
# best practice to make sure that reads with incorrect indexes do not
# impact variant calling or assignment of gene expression counts.
check_for_duplicate = function(result1, result2) {
    check = data.frame(Id1 = result1$Id, Id2 = result2$Id)
    if (anyDuplicated(check) != 0) {
        d = anyDuplicated(check)
        for (i in seq(1,length(d))) {
            id = result2[d[i], ]$Id
            lane_to_change = result2 %>% filter(Lane == result2[d[i], ]$Lane)
            lanes_to_keep  = result2 %>% filter(Lane != result2[d[i], ]$Lane)
            j = which(lane_to_change$Id == id) %>% as.numeric()
            if (j < nrow(lane_to_change)) {
                temp = lane_to_change[j, ]
                lane_to_change[j, ] = lane_to_change[j + 1, ]
                lane_to_change[j + 1, ] = temp
            } else {
                temp = lane_to_change[j, ]
                print(temp)
                lane_to_change[j, ] = lane_to_change[1, ]
                lane_to_change[1, ] = temp
            }
            result = bind_rows(lane_to_change, lanes_to_keep) %>% arrange(Lane)
        }
        return (result)
    } else{
        return (result2)
    }
}




is_a_prime_number = function (sample_number) {
    result = isPrime(sample_number) %>% as.numeric()
    return(result)
}



final_result = function(index_df,
                        sample_number,
                        mplex_level,
                        platform,
                        metric,
                        d,
                        thrs_size_comb = 120,
                        max_iteration = 50,
                        method = "greedy_exchange") {
    if (mplex_level<= nrow(index_df)){
    result1 = get_result(
        index_df,
        sample_number,
        mplex_level,
        platform,
        metric,
        d,
        thrs_size_comb,
        max_iteration,
        method
    )}else {
        result1 = NULL
        display_message(
        "the number of barcodes is not sufficient for the multiplexing level")
    }

    if(!is.null(result1)){
    result1 = data.frame(
        sample = seq(1,sample_number) %>% as.character(),
        Lane = result1$Lane %>% as.character(),
        Id = result1$Id %>% as.character(),
        stringsAsFactors = FALSE
    )

    result1 = left_join(result1, select(index_df, Id, sequence), by = "Id")
    return (result1)
    } else {
        return(NULL)
    }
}



final_result_dual = function(   index_df_1,
                                index_df_2,
                                sample_number,
                                mplex_level,
                                platform = 4,
                                metric = NULL,
                                d = 3,
                                thrs_size_comb = 120,
                                max_iteration = 50,
                                method = "greedy_exchange") {
    result1 = get_result(
        index_df_1,
        sample_number,
        mplex_level,
        platform,
        metric,
        d,
        thrs_size_comb,
        max_iteration,
        method
    )
    print(result1)
    result2 = get_result(
        index_df_2,
        sample_number,
        mplex_level,
        platform,
        metric,
        d,
        thrs_size_comb,
        max_iteration,
        method
    )
    print(result2)
    if(!is.null(result1) && !is.null(result2)){
        result2 = check_for_duplicate(result1, result2)
        result1 = left_join(result1,
                            select(index_df_1, Id, sequence),
                            by = "Id")
        print(result1)
        result2 = left_join(result2,
                            select(index_df_2, Id, sequence),
                            by = "Id")
        print(result2)
        result = data.frame(
            sample = seq(1,sample_number) %>% as.numeric(),
            Lane = result1$Lane %>% as.character(),
            Id1 = result1$Id %>% as.character(),
            sequence1 = result1$sequence %>% as.character(),
            Id2 = result2$Id %>% as.character(),
            sequence2 = result2$sequence %>% as.character(),
            stringsAsFactors = FALSE
        ) %>% arrange(sample)
        result$sample = as.character(result$sample)
        return(result)
    } else {
        return(NULL)
    }
}


display_message <- function (a_message) {
    assign("error_message", a_message, envir = .GlobalEnv)
    print(a_message)
}

check_input_dataframe <- function(df, column_names){
    
    if ( !("data.frame" %in% is(df))) {
        display_message(paste(  "wrong input format: a dataframe is",
                                "expected"))
    } else if (!all(column_names %in% names(df))) {
            missing_column_names <- setdiff(column_names, names(df))
            display_message(paste(  "wrong input format,",
                                    "missing column", missing_column_names))
    }
}

