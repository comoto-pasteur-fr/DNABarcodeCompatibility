#!/usr/bin/Rscript --no-environ
options(width=120)
## Setup  I/O ####
argv = commandArgs(TRUE)

if((length(argv) < 2) | (length(argv) > 13)){
  stop("Usage: time_simulation_s.R <rep_number> <thrs_size_comb> <max_iteration> <inputfile> <nb_lane> <barcode_number> <chemistry> <outputfile1> <outputfile2> <outputfile3>")
}else{
  rep_number=as.integer(argv[1])
  thrs_size_comb=as.integer(argv[2])
  max_iteration=as.integer(argv[3])
  inputfile=normalizePath(argv[4])
  nb_lane=as.integer(argv[5])
  barcode_number=as.integer(argv[6])
  chemistry=as.integer(argv[7])
  outputfile1=normalizePath(argv[8])
  outputfile2=normalizePath(argv[9])
  outputfile3=normalizePath(argv[10])
  outputfile4=normalizePath(argv[11])
  mplex_level=as.integer(argv[12])
  algorithm=as.character(argv[12])
}


# max_iteration=0
# rep_number=1
# chemistry=4
# thrs_size_comb=130
# combination_m=DNABarcodeCompatibility::get_all_combinations(index_df = DNABarcodeCompatibility::IlluminaIndexes[1:18,], mplex_level = 4, chemistry = 4)
# nb_lane=4
# barcode_number=18
# mplex_level = 4
# algorithm="greedy_exchange"

## Parametrize the optimize_combination function for the sake of the simulation
optimize_combinations = function (combination_m, nb_lane, index_number, thrs_size_comb=80, max_iteration=10, algorithm){
  # browser()
  if (nrow(as.matrix(combination_m)) == 0){
    DNABarcodeCompatibility:::display_message("No combinations have been found")
  }else {
    if(is.numeric(index_number)){
      if (is.numeric(nb_lane)){
        
        max = DNABarcodeCompatibility:::entropy_max(index_number, ncol(combination_m) * nb_lane)
        print(paste("Theoretical max entropy:",round(max, 3)))
        if(nb_lane < nrow(combination_m)){
          if (nrow(combination_m) > thrs_size_comb){ 
            init_combination_m = combination_m[sample(1:nrow(combination_m),thrs_size_comb),]
            random_pick_combination_init_m = random_pick_combination_m = init_combination_m[sample(1:nrow(init_combination_m), size = nb_lane),]
            a_combination = DNABarcodeCompatibility:::recursive_entropy(init_combination_m,nb_lane, method=algorithm)
            print(paste("Entropy of the entire subset:",round(DNABarcodeCompatibility:::entropy_result(init_combination_m), 3)))
            print(paste("Entropy of a random pick:",round(DNABarcodeCompatibility:::entropy_result(random_pick_combination_m), 3)))
            i = 0
            while ((i < max_iteration) && (DNABarcodeCompatibility:::entropy_result(a_combination) < max) ){
              temp_combination = DNABarcodeCompatibility:::recursive_entropy(combination_m[sample(1:nrow(combination_m),thrs_size_comb),], nb_lane, method=algorithm)
              if (DNABarcodeCompatibility:::entropy_result(temp_combination) > DNABarcodeCompatibility:::entropy_result(a_combination) ){
                a_combination = temp_combination
              }
              temp_random_combination=init_combination_m[sample(1:nrow(init_combination_m), size = nb_lane),]
              if (DNABarcodeCompatibility:::entropy_result(temp_random_combination) > DNABarcodeCompatibility:::entropy_result(random_pick_combination_m) ){
                random_pick_combination_m = temp_random_combination
              }
              i = i+1
            }
          } else {
            init_combination_m = combination_m
            random_pick_combination_init_m = random_pick_combination_m = init_combination_m[sample(1:nrow(init_combination_m), size = nb_lane),]
            print(paste("Entropy of the entire subset:",round(DNABarcodeCompatibility:::entropy_result(init_combination_m), 3)))
            print(paste("Entropy of a random pick:",round(DNABarcodeCompatibility:::entropy_result(random_pick_combination_m), 3)))
            a_combination = DNABarcodeCompatibility:::recursive_entropy(combination_m,nb_lane, method=algorithm)
            i = 0
            while ((i < max_iteration) && (DNABarcodeCompatibility:::entropy_result(a_combination) < max) ){
              temp_combination = DNABarcodeCompatibility:::recursive_entropy(combination_m, nb_lane, method=algorithm)
              if (DNABarcodeCompatibility:::entropy_result(temp_combination) > DNABarcodeCompatibility:::entropy_result(a_combination) ){
                a_combination = temp_combination
              }
              temp_random_combination=init_combination_m[sample(1:nrow(init_combination_m), size = nb_lane),]
              if (DNABarcodeCompatibility:::entropy_result(temp_random_combination) > DNABarcodeCompatibility:::entropy_result(random_pick_combination_m) ){
                random_pick_combination_m = temp_random_combination
              }
              i = i+1
            }
          }
        } else {
          # stop("nb_lane > nrow(combination_m): case not implemented for simulations")
          n = nb_lane - nrow(combination_m)
          combination_m = combination_m[sample(1:nrow(combination_m),nrow(combination_m)),,drop = F]
          part_combination = combination_m[sample(1:nrow(combination_m),n, replace = TRUE),]
          a_combination = rbind(combination_m, part_combination)
          init_combination_m = random_pick_combination_init_m = random_pick_combination_m = a_combination
          print(paste("Entropy of the entire subset:",round(DNABarcodeCompatibility:::entropy_result(init_combination_m), 3)))
          print(paste("Entropy of a random pick:",round(DNABarcodeCompatibility:::entropy_result(random_pick_combination_m), 3)))
          i = 0
          while ((i < 30) && (DNABarcodeCompatibility:::entropy_result(a_combination) < max) && n > 0){
            part_combination = combination_m[sample(1:nrow(combination_m),n, replace = TRUE),]
            temp_combination = rbind(combination_m, part_combination)
            if (DNABarcodeCompatibility:::entropy_result(temp_combination) > DNABarcodeCompatibility:::entropy_result(a_combination) ){
              a_combination = temp_combination
            }
            i = i + 1
          }
          a_combination = a_combination[sample(1:nrow(a_combination),nrow(a_combination)),]
        }
        print(paste("Entropy of the optimized set:", round(DNABarcodeCompatibility:::entropy_result(a_combination),3)))
        return(list("init_comb"=init_combination_m, "opt_comb"=a_combination, "random_comb"=random_pick_combination_init_m, "opt_random_comb"=random_pick_combination_m))
      } else {
        display_message("Please enter a number as nb_lane")
      }
    } else {
      display_message("Please enter a number as index_number")
    }
  }
}


## Main ####
combination_m <- local(get(load(inputfile)))

S_max = DNABarcodeCompatibility:::entropy_max(barcode_number, ncol(combination_m) * nb_lane)
elapsed <- as.numeric(system.time({out_comb <- optimize_combinations(combination_m, nb_lane, index_number=barcode_number, thrs_size_comb=thrs_size_comb, max_iteration=max_iteration, algorithm)}))[3]
S_init=DNABarcodeCompatibility:::entropy_result(unlist(out_comb$init_comb))
S_opt=DNABarcodeCompatibility:::entropy_result(unlist(out_comb$opt_comb))
S_random=DNABarcodeCompatibility:::entropy_result(unlist(out_comb$random_comb))
S_random_opt=DNABarcodeCompatibility:::entropy_result(unlist(out_comb$opt_random_comb))

# print(paste(rep_number, round(elapsed,0), nrow(combination_m), barcode_number, thrs_size_comb, chemistry, nb_lane), quote=F)
data.frame(rep_number=rep_number, time=round(elapsed,0), nb_comp_comb=nrow(combination_m), barcode_set_size=barcode_number, barcode_subset_size=nrow(unlist(out_comb$init_comb)), chemistry=chemistry, nb_lane=nb_lane, mplex_level=mplex_level, algorithm=algorithm)


## Save barcode occurrence from random pick, with input simulation parameters
df <- dplyr::mutate(as.data.frame(table(unlist(out_comb$random_comb))),
                    rep_number=rep_number,
                    thrs_size_comb=nrow(unlist(out_comb$init_comb)),
                    barcode_number=barcode_number,
                    chemistry=chemistry, 
                    nb_lane=nb_lane,
                    mplex_level=mplex_level, 
                    algorithm=algorithm)
write.table(df, file = outputfile1 , quote = F, row.names = F, col.names = F, append = T)


## Save barcode occurrence from optimized set, with input simulation parameters
df <- dplyr::mutate(as.data.frame(table(unlist(out_comb$opt_comb))),
                    rep_number=rep_number,
                    thrs_size_comb=nrow(unlist(out_comb$init_comb)),
                    barcode_number=barcode_number,
                    chemistry=chemistry, 
                    nb_lane=nb_lane,
                    mplex_level=mplex_level,
                    algorithm=algorithm)
write.table(df, file = outputfile2 , quote = F, row.names = F, col.names = F, append = T)


## Save simulation parameters and elasped time
df <- data.frame(time=round(elapsed,3), nb_comp_combinations=nrow(combination_m), rep_number=rep_number, barcode_number=barcode_number, 
                 thrs_size_comb=nrow(unlist(out_comb$init_comb)), S_max, S_init, S_random, S_opt, S_random_opt, chemistry=chemistry, nb_lane=nb_lane, mplex_level=mplex_level, algorithm=algorithm)
write.table(df, file = outputfile3 , quote = F, row.names = F, col.names = F, append = T)

## Save barcode occurrence from random pick + optimization, with input simulation parameters
df <- dplyr::mutate(as.data.frame(table(unlist(out_comb$opt_random_comb))),
                    rep_number=rep_number,
                    thrs_size_comb=nrow(unlist(out_comb$init_comb)),
                    barcode_number=barcode_number,
                    chemistry=chemistry, 
                    nb_lane=nb_lane,
                    mplex_level=mplex_level, 
                    algorithm=algorithm)
write.table(df, file = outputfile4 , quote = F, row.names = F, col.names = F, append = T)
