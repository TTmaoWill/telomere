library(tidyverse)
# filename <- 'CB4856\\IIIL_CB4856.ape'
#' Read sequences from APE and FASTA files
#'
#' This function reads DNA sequences from files with extensions `.ape` and `.fasta`.
#' It returns a tibble containing the name, length, and sequence of each DNA sequence.
#'
#' @param filenames A character vector of file paths to the `.ape` and `.fasta` files.
#' @return A tibble with columns: `name` (the name of the sequence), `length` (the length of the sequence), and `sequence` (the DNA sequence).
#' @importFrom tibble tibble add_row
#' @importFrom Biostrings readDNAStringSet
#' @importFrom stringr str_split
#' @importFrom dplyr %>%
#' @examples
#' \dontrun{
#' filenames <- c("sequence1.ape", "sequence2.fasta")
#' sequences <- read_ape_fasta(filenames)
#' print(sequences)
#' }
#' @export

read_ape_fasta <- function(filenames){
  out <- tibble(name = character(), length = numeric(), sequence = character())
  for (filename in filenames){
    type <- strsplit(filename,'\\.')[[1]]
    type <- type[length(type)]
    if (type == 'ape') {
      ls <- readLines(con = file(filename, open='r'))
      
      seq <- ls[(which(ls == "ORIGIN") + 1):(length(ls)-1)]
      seq <- strsplit(paste(seq,collapse=''),split=' ')[[1]]
      seq <- seq[seq != '']
      seq <- seq[is.na(as.numeric(seq))]
      seq <- toupper(paste(seq,collapse=''))
      name <- strsplit(filename,'\\.ape')[[1]][1]
      length <- nchar(seq)
      
      stopifnot(nchar(seq) == length)
    } else {
      seq <- toupper(as.character(readDNAStringSet(filename)))
      
      name <- strsplit(filename,'\\.fasta')[[1]][1]
      length <- nchar(seq)
    }
    out <- out %>% add_row(name = name, length = length, sequence = seq)
  }
  out
}

