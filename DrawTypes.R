# first created: 28 Sep 2023
# last updated: 18 Oct 2023
# author: Andreas Uthemann

DrawTypes <- function(M, init, p_class, seed = 1){

    # draw submitter types (weak or strong) for M sample paths
    # returns a MM*S matrix of types
    # types[m,s] = 1 : submitter s for draw m is weak (0 if strong)
    # MM <= M as duplicate rows are dropped

    old <- .Random.seed
    on.exit( { .Random.seed <<- old } ) 

    set.seed(seed)

    types_0 <- init$weak   # types_0[i] == TRUE if submitter i is classified as weak using kmeans clustering on ts stddev of submissions
    S <- length(types_0)   # number of submitters

    types <- matrix(NA, nrow = M, ncol = S)
    for (s in 1:S){
        q <- ifelse(types_0[s], p_class, 1 - p_class)  # if s is classed as weak, P(type_s == 1) = 1 - p_class, if strong then P(type_s == 1) = p_class > 1/2
        types[,s] = rbinom(M, 1, q)
    }

    types <- unique(types)  # remove duplicate rows
    return(types)
}