
distance <- function(sequence1, sequence2, metric = c("hamming", "seqlev", "levenshtein", 
    "phaseshift"), cost_sub = 1, cost_indel = 1) {
    metric <- match.arg(metric)
    
    if ((metric == "hamming" || metric == "phaseshift") && (nchar(sequence1) != nchar(sequence2))) {
        msg <- sprintf("%s metric is not defined for strings of unequal length (|sequence1|==%i but |sequence2|==%i).", 
            metric, nchar(sequence1), nchar(sequence2))
        warning(msg)
    }
    
    return(.distance(sequence1, sequence2, metric, cost_sub, cost_indel))
}
 
