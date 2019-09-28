##' global alignment 
##'
##' 
##' @title global_aln
##' @param X sequence 1
##' @param Y sequence 2
##' @param score score system
##' @return smuAlign object
##' @export
##' @author Guangchuang Yu
global_aln <- function(X, Y, score = list(match = 5, mismatch = -2, indel = -6)) {
    seq.x <- unlist(strsplit(X, ''))
    seq.y <- unlist(strsplit(Y, ''))
    seq.x <- c(0,seq.x)
    seq.y <- c(0,seq.y)
    match <- score$match
    mismatch <- score$mismatch
    indel <- score$indel
    ## initial the score matrix
    mat <- matrix(NA, length(seq.x), length(seq.y))
    mat[,1] <- sapply(1:length(seq.x)-1, function(x) x * indel)
    mat[1,] <- sapply(1:length(seq.y)-1, function(x) x * indel) 
    ## The dynamic programming, global alignment recursion
    for (i in 2:length(seq.x)) {
        for (j in 2:length(seq.y)){
            ## seq.x[i] , seq.y[j] are aligned
            if ( seq.x[i] == seq.y[j]) {
                mat[i,j] <- mat[i-1, j-1] + match
            } else {
                mat[i,j] <- mat[i-1, j-1] + mismatch
            }         ## seq.x[i] aligned to -
            sc <- mat[i-1,j] + indel
            if (sc > mat[i,j])
                mat[i,j] = sc 
            ## seq.y[j] aligned to -
            sc <- mat[i,j-1] + indel   
            if (sc > mat[i,j])
                mat[i,j] = sc
        }
    } 
    ## Traceback
    i <- length(seq.x)
    j <- length(seq.y)
    ax <- character()
    ay <- character() 
    while (i > 1 && j >1){     
        ## case 1: best was seq.x[i] aligned to seq.y[j]
        sc <- mat[i-1,j-1]     
        if (seq.x[i] == seq.y[j]) {
            sc <- sc + match
        } else {
            sc <- sc + mismatch
        } 
        if (sc == mat[i,j]) {
            ax <- c(seq.x[i], ax)
            ay <- c(seq.y[j], ay)
            i <- i -1
            j <- j-1
            next
        }     
        ## case 2: best was seq.x[i] aligned to -
        if ((mat[i-1,j] + indel) == mat[i,j]) {
            ax <- c(seq.x[i], ax)
            ay <- c("-", ay)
            i <- i-1
            next
        }     
        ## case 3: best was seq.y[j] aligned to -
        if ((mat[i,j-1] + indel) == mat[i,j]) {
            ax <- c("-", ax)
            ay <- c(seq.y[j], ay)
            j <- j-1
            next
        }
    }

    ax <- paste(ax, collapse='')
    ay <- paste(ay, collapse='')

    structure(
        list(seq = c(X, Y),
             aln = c(ax, ay),
             score = score,
             matrix = mat
             ),
        class = "smuAlign"
    )
}

##' @method print smuAlign
##' @export
print.smuAlign <- function(x, ...) {
    mat <- x$matrix
    cat("Sequence X: ", x$seq[1],"\n")
    cat("Sequence Y: ", x$seq[2],"\n")
    cat("Scoring system: ", x$score$match, " for match; ",
        x$score$mismatch, " for mismatch; ",
        x$score$indel, " for gap", "\n\n")

    cat("Dynamic programming matrix:\n")
    print(mat)
    cat("\nAlignment:\n")
    cat(x$aln[1], "\n")
    cat(paste0(rep("|", nchar(x$aln[1])), collapse = ""), "\n")
    cat(x$aln[2],"\n\n")
    cat("Optimum alignment mat: ", mat[length(mat)],"\n")
}
