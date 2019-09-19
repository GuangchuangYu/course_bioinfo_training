##' @importFrom treeio read.fasta
##' @export
treeio::read.fasta

##' @importFrom ape base.freq
##' @export
ape::base.freq

##' write DNA or protein sequences to FASTA file
##'
##' 
##' @title write.fasta
##' @param x DNA or protein sequence object
##' @param file output file
##' @return NULL
##' @export
write.fasta <- function(x, file) {
    ape::write.FASTA(x, file)
}
