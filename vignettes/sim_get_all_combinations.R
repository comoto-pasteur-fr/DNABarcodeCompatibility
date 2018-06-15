#!/usr/bin/Rscript --no-environ

## Setup  I/O ####
argv = commandArgs(TRUE)

if((length(argv) < 2) | (length(argv) > 5)){
  stop("Usage: time_simulation_s.R <mplex_level> <barcode_number> <chemistry> <outputfile1> <outputfile2>")
}else{
  mplex_level=as.integer(argv[1])
  barcode_number=as.integer(argv[2])
  chemistry=as.integer(argv[3])
  outputfile1=normalizePath(argv[4])
  outputfile2=normalizePath(argv[5])
}

## Main ####
nb_total_combinations <- as.numeric(choose(barcode_number, mplex_level))
index <- dplyr::sample_n(DNABarcodeCompatibility::IlluminaIndexes, nrow(DNABarcodeCompatibility::IlluminaIndexes), replace = FALSE)
elapsed <- as.numeric(system.time(compatible_comb <- DNABarcodeCompatibility::get_all_combinations(index_df = index[1:barcode_number,], mplex_level = mplex_level, chemistry = chemistry))[3])
nb_comp_combinations <- nrow(compatible_comb)

## Save barcode occurrence from the set of compatible barcodes, with input simulation parameters
df <- dplyr::mutate(as.data.frame(table(compatible_comb)),
                    mplex_level=mplex_level,
                    barcode_number=barcode_number,
                    chemistry=chemistry)
write.table(df, file = outputfile1 , quote = F, row.names = F, col.names = F, append = T)



## Save simulation parameters and elasped time
print(paste(round(elapsed,3), nb_total_combinations, nb_comp_combinations, mplex_level, barcode_number, chemistry), quote=F)
df <- data.frame(time=round(elapsed,3), nb_total_combinations=nb_total_combinations, nb_comp_combinations=nb_comp_combinations, mplex_level=mplex_level, barcode_number=barcode_number, chemistry=chemistry)
write.table(df, file = outputfile2 , quote = F, row.names = F, col.names = F, append = T)

