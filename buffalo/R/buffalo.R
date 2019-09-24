#' Print message with buffalo
#' 
#' @title buffalo_say
#' @param who Who love Bioinformatics? Defaults to SMU
#' @importFrom emojifont emoji
#' @importFrom cowsay say
#' @export
buffalo_say <- function(who="SMU") {
    msg <- paste(who, emoji("heart"), 
            " Bioinformatics", emoji("muscle"))
    say(msg, 'buffalo', 
        by_color = rainbow(10))
}