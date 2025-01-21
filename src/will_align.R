# Borrowed from Evan
# Will Tang

# Install required packages if not already installed
# required_packages <- c("seqinr", "Biostrings", "stringr", "gtsummary", "ggplot2", "ggthemes", "ggrepel", "beepr", "tidyverse", "stringi")

# new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
# if(length(new_packages)) install.packages(new_packages)

# Load required libraries
library(seqinr)
library(Biostrings)
library(stringr)
library(gtsummary)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(tidyverse)
library(stringi)
source(file.path('src','read_seqs.R'))

# Function to get edit operations between two strings
#' @param string1 The first string
#' @param string2 The second string
#' @return A list of edit operations
getOp <- function(string1, string2) {
  n <- nchar(string1)
  m <- nchar(string2)
  d <- matrix(0, n + 1, m + 1)
  
  for (i in 1:(n + 1))
    d[i, 1] <- i
  
  for (j in 1:(m + 1))
    d[1, j] <- j
  
  for (j in 2:(m + 1)) {
    for (i in 2:(n + 1)) {
      if (substr(string1, i - 1, i - 1) != substr(string2, j - 1, j - 1))
        d[i, j] <- min(d[i - 1, j] + 1, d[i, j - 1] + 1, d[i - 1, j - 1] + 1)
      else
        d[i, j] <- d[i - 1, j - 1]
    }
  }
  
  i <- n + 1
  j <- m + 1
  operations <- list()
  
  while (i > 1 && j > 1) {
    if (substr(string1, i - 1, i - 1) != substr(string2, j - 1, j - 1)) {
      if (d[i, j] == d[i - 1, j] + 1) {
        operations <- c(operations, paste("Delete", i - 1))
        i <- i - 1
      } else if (d[i, j] == d[i, j - 1] + 1) {
        operations <- c(operations, paste("Insert", j - 1))
        j <- j - 1
      } else {
        operations <- c(operations, paste("Substitute", i - 1))
        i <- i - 1
        j <- j - 1
      }
    } else {
      i <- i - 1
      j <- j - 1
    }
  }
  
  if (i > 1) {
    for (k in 1:(i - 1)) {
      operations <- c(operations, paste("Delete", k))
    }
  } else if (j > 1) {
    for (k in 1:(j - 1)) {
      operations <- c(operations, paste("Insert", k))
    }
  }
  operations <- rev(operations)
  return(operations)
}

# Function to get the type of edit operations
#' @param opList A list of edit operations
#' @return The number of unique operation types
getOpType <- function(opList) {
  a <- c()
  for (i in 1:length(opList)) {
    a[i] <- strsplit(opList[[i]], ' ')[[1]][1]
  }
  return(length(unique(a)))
}

# Function to check if there is a substitution in the edit operations
#' @param opList A list of edit operations
#' @return TRUE if there is a substitution, FALSE otherwise
hasSub <- function(opList) {
  a <- c()
  for (i in 1:length(opList)) {
    a[i] <- strsplit(opList[[i]], ' ')[[1]][1]
  }
  return('Substitute' %in% a)
}

# Function to get edit operations for a given string
#' @param string1 The string to check
#' @param tandem The tandem repeat sequence
#' @return A string describing the edit operations
getEditOperations <- function(string1,tandem) {
  ref1 <- tandem
  ref2 <- paste0(tandem, tandem)
  
  op1 <- getOp(ref1, string1)
  op2 <- getOp(ref2, string1)
  
  if (nchar(string1) <= 6) {
    operations <- paste(unlist(op1), collapse = ', ')
  } else {
    if (getOpType(op1) < getOpType(op2)) {
      operations <- paste(unlist(op1), collapse = ', ')
    } else if (getOpType(op1) > getOpType(op2)) {
      operations <- paste(unlist(op2), collapse = ', ')
    } else {
      operations <- ifelse(length(op1) < length(op2), paste(unlist(op1), collapse = ', '), paste(unlist(op2), collapse = ', '))
    }
  }
  return(operations)
}

# Function to check the type of mutation
#' @param str The string to check
#' @param tandem The tandem repeat sequence
#' @return A string describing the type of mutation
check_type <- function(str,tandem) {
  if (is.na(match(str, tandem))) {
    if (nchar(str) > 12) {
      'Long Insertion'
    } else {
      getEditOperations(str, tandem)
    }
  } else {
    'No Mutation'
  }
}

# Function to find mutations in the sequences
#' @param data_name The name of the data folder
#' @param tandem The tandem repeat sequence
find_muts <- function(data_name, tandem) {
  dir <- file.path('data','raw',data_name)
  files <- list.files(dir)
  setwd(dir)
  seqs <- read_ape_fasta(files)
  setwd(file.path('..','..','processed'))

  if (!dir.exists(paste0(data_name,'_result'))) {
    dir.create(paste0(data_name,'_result'))
  }
  
  for (i in 1:nrow(seqs)) {
    print(paste0('Processing ', i, ' of ', nrow(seqs)))
    file <- files[i]
    row <- seqs[i, ]
    dnaseq <- row$sequence
    length <- row$length
    name <- row$name
    
    print(name)

    t_index <- as.vector(gregexpr(tandem, dnaseq)[[1]])
    freq <- length(t_index[t_index > 0])
    print(freq)
    
    if (freq < 10) {
      original_sequence <- "ATCG"
      
      # Define the complement mapping
      complement_mapping <- "TAGC"
      dnaseq <- chartr(original_sequence, complement_mapping, dnaseq)
      dnaseq <- stri_reverse(dnaseq)
      t_index <- as.vector(gregexpr(tandem, dnaseq)[[1]])
    }
    print(t_index)

    if (t_index[1] != -1){
      checked_index <- c()
      
      if (length(t_index) != 1) {
        for (i in 1:(length(t_index) - 1)) {
          if (t_index[i + 1] - t_index[i] >= 6) {
            checked_index <- append(checked_index, t_index[i])
          }
        }
      }
      checked_index <- append(checked_index, t_index[length(t_index)])
      t_index <- checked_index
      
      if (t_index[1] != 1) {
        splited <- c(substr(dnaseq, 1, t_index[1] - 1))
      } else {
        splited <- c()
      }
      
      for (i in 1:(length(t_index) - 1)) {
        if (t_index[i + 1] - t_index[i] <= 9) {
          splited <- append(splited, substr(dnaseq, t_index[i], t_index[i + 1] - 1))
        } else {
          splited <- append(splited, substr(dnaseq, t_index[i], t_index[i] + 5))
          splited <- append(splited, substr(dnaseq, t_index[i] + 6, t_index[i + 1] - 1))
        }
      }
      if (nchar(dnaseq) - t_index[length(t_index)] <= 9) {
        splited <- append(splited, substr(dnaseq, t_index[length(t_index)], nchar(dnaseq)))
      } else {
        splited <- append(splited, substr(dnaseq, t_index[length(t_index)], t_index[length(t_index)] + 5))
        splited <- append(splited, substr(dnaseq, t_index[length(t_index)] + 6, nchar(dnaseq)))
      }
      temp <- paste(splited, collapse = '')
      
      out <- tibble(sequence = splited)
      out$type <- sapply(out$sequence, check_type,tandem)
      out$length <- sapply(out$sequence, nchar)
      out$pos <- lag(out$length)
      out[1, ]$pos <- 1
      out <- out %>% mutate(pos = cumsum(pos)) %>% select(-length)
      
      write.table(out, file.path(paste0(data_name,'_result'), paste0(name, '_processed.txt')), quote = F, row.names = F,sep="\t")
      
      # Generate new output table
      new_out <- tibble(start = integer(), end = integer(), type = character())
      new_out <- add_row(new_out, start = 1, end = nchar(dnaseq), type = 'Range')
      for (j in 1:nrow(out)) {
        if (out$type[j] != 'No Mutation') {
          ops <- strsplit(out$type[j], ', ')[[1]]
          for (op in ops) {
            print(op)
            parts <- strsplit(op, ' ')[[1]]
            pos <- as.integer(parts[2])
            if (parts[1] == 'Substitute') {
              new_out <- add_row(new_out, start = out$pos[j] + pos - 1, end = out$pos[j] + pos, type = 'Substitution')
            } else if (parts[1] == 'Insert') {
              new_out <- add_row(new_out, start = out$pos[j] + pos - 1, end = out$pos[j] + pos, type = 'Insertion')
            } else if (parts[1] == 'Delete') {
              new_out <- add_row(new_out, start = out$pos[j] + pos - 1, end = out$pos[j] + pos, type = 'Deletion')
            } else if (parts[1] == 'Long' && parts[2] == 'Insertion') {
                new_out <- add_row(new_out, start = out$pos[j], end = out$pos[j] + nchar(out$sequence[j]) - 1, type = 'Long Insertion')
            }
          }
        }
      }
      
      write.table(new_out, file.path(paste0(data_name,'_result'), paste0(name, '_mut.txt')), quote = F, row.names = F,sep="\t")
    }
    else {
      write.table(data.frame(), file.path(paste0(data_name,'_result'), paste0(name, '_processed.txt')), quote = F, row.names = F,sep="\t")
      write.table(data.frame(), file.path(paste0(data_name,'_result'), paste0(name, '_mut.txt')), quote = F, row.names = F,sep="\t")
    } 
  }
}


plot_mutations <- function(data_name) {
  # List all files ending with "_mut.txt" in the folder
  folder_path <- file.path(paste0(data_name,'_result'))
  files <- list.files(path = folder_path, pattern = "_mut\\.txt$", full.names = TRUE)
  
  # Check if there are matching files
  if (length(files) == 0) {
    message("No files ending with '_mut.txt' found in the folder.")
    return(NULL)
  }
  
  # Iterate over each file
  for (file in files) {
    # Read the data
    data <- read_tsv(file)
    
    # Check if the data is empty
    if (nrow(data) == 0) {
      # Generate a white plot
      plot <- ggplot() +
        theme_void() +
        labs(title = paste0("No mutations in ", basename(file)))
    } else {
      # Assign colors to mutation types
      data$color <- factor(data$type, 
                           levels = c("Range", "Long Insertion", "Substitution", "Deletion", "Insertion"),
                           labels = c("black", "purple", "blue", "red", "green"))
      
      # Generate the plot
      plot <- ggplot() +
        geom_segment(data = data, aes(x = start, xend = end, y = 0, yend = 0, color = color), size = ifelse(data$type == "Range", 1, 10)) +
        scale_color_identity() +
        labs(x = "Position", y = "", title = paste0("Mutations in ", basename(file))) +
        theme_minimal() +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.background = element_rect(fill = "white"), panel.grid = element_blank())
    }
    
    # Save the plot as PNG
    output_file <- file.path(folder_path, paste0(tools::file_path_sans_ext(basename(file)), "_plot.png"))
    ggsave(output_file, plot, width = 10, height = 5)
  }
  
  message("Plots have been created and saved.")
}

# Main script to find mutations for different data sets
# name <- 'hawaiian'
# name <- 'CB4856'
name <- 'N2'
# for (name in c('hawaiian', 'CB4856', 'N2')) {
#   find_muts(name, tandem)
# }
tandem <- 'TTAGGC'
find_muts(name, tandem)
plot_mutations(name)
