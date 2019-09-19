##' extract accession number and sequence from genbank file
##'
##' 
##' @title read.genbank
##' @param file input genbank file
##' @return sequence object
##' @export
##' @author Guangchuang Yu
read.genbank <- function(file) {
    x <- readLines(file)
    acc <- sub("\\w+\\s+(\\w+)$", "\\1", x[grep("ACCESSION", x)])
    i <- grep("ORIGIN", x)
    ss <- x[(i+1):length(x)]
    ss <- ss[1:(grep("//", ss) -1)]
    ss <- gsub("\\s+\\d+", "", ss)
    ss <- gsub("\\s+", "", ss)
    ss <- paste0(ss, collapse = "")
    f <- tempfile(fileext = ".fasta")
    cat(">",  file = f,append = TRUE)
    cat(acc,  file = f,append = TRUE, sep = "\n")
    cat(ss,  file = f,append = TRUE, sep = "\n")
    read.fasta(f)
}

