# Script for plotting deletions and insertions in telomere sequence. ----
# Written by Evan Lister-Shimauchi, 2022, in the lab of Shawn Ahmed
# Load required libraries (not sure if I used all of these in the final version)
library(seqinr)
library(Biostrings)
library(stringr)
library(gtsummary)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(beepr)
library(tidyverse)

# Manually paste in the actual sequence of your telomere you want analyzed. ----
# This current version requires the sequence be in the correct orientation (TTAGGC rather than GCCTAA). 
# A future improvement would be a few lines to allow either orientation to be pasted in.
# A further improvement would be to make R run this script on all FASTA files in a given folder.

telomereseq <- c("AGGGACCagtaatacaattttcaaatggcTATATATGtatcataaatttttctattacttgtaataaaatatttatagtaATCATTTctgatcgatttttttcctgaaaaattatgatataAGGGTTTAAATCATAGTTTGGAGCTCTTACCAtttcttaggcttaggcttttaggcttaggcttaggcttaggcttagggcttaggcttaggcttaggcttaggcttaggctttaggcttaggcttagttaggcttaggcttaggctttagggcttaggcttaggccttaggcttaggctttaggcttaggcttaggcttagggcttagggcttaggcttaggcttaggcttaggcttaggcttaggcttggctTAGGGCCTTAGGGCTTATTAGGCTTGCTTAGCTAGGGCTTAGGCGGCGgctggcttaggcttaggaggCTTAGgctggcttaggcttaggcttaggcttaggcttaggcttagggcttagggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggctttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggctttaggcttaggcttaggcttaggcttaggcttaggcttagggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttagggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggggcttaggcttaggcttaggcttaggcttaggcttaggcaggcttaggcttaggcttagcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttagggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggctttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttagggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttagcgCTTAGGCTTagcggcttaggcttaggcagcttaggcttaggcttaggcttggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggctttaggcttaggcttaggcttaggcttagggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttagggcttaggcttaggcttaggcttaggcttagggcttaggcttaggcttaggcttaggcttagggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttagggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttagggcttaggcttaggcttaggcttagggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttagggcttaggcttaggcttaggcttaggcttaggcttagggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggctttaggcttaggcttaggcttaggcttaggcttaggctttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggccttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcgcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttagcagcttaggcttaggcaggcttaggcttaggcttaggcttaggcttaggcttaggctttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggctttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggctttaggctttaggcttaggcttaggcttagggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttaggcttagcttaggcttaggcttaggcttaggcttaggcttaggcttagcttaggcttaggcttaaggcTTAG")
telomereseq <- toupper(telomereseq)

# Enter a name that will appear in the final graph and output text file.

telomerename = "IIL_DL238"

# Split the telomere sequence string into a listof single characters

telomeresplit <- unlist(strsplit(as.character(telomereseq), split = ""))

# Automatically make reference sequence of similar size with perfect telomere repeats, with a few extra tacked on. ----

refseq <- c()

for (i in 1:(2 + ceiling(length(telomeresplit) / 6))) {
  refseq <- append(refseq,"TTAGGC")
}

refseq <- paste(unlist(strsplit(as.character(refseq), split = "")),collapse = "")

# Align actual telomere sequence (telomereseq) to idealized reference telomere sequence (refseq). ----
# This relies on the pairwiseAlignment function from the Biostrings package.
# The current strategy is to run the alignment twice (with the telomere sequence as the reference or subject), but there may be a more elegant way to do this.

alignmentR <- pairwiseAlignment(refseq, telomereseq, scoreOnly = FALSE)
alignmentV <- pairwiseAlignment(telomereseq, refseq, scoreOnly = FALSE)

telomereseq2 <- unlist(strsplit(as.character(alignmentV), split = ""))
refseq2 <- unlist(strsplit(as.character(alignmentR), split = ""))

changeseqD <- c()
changeseqI <- c()
for (i in 1:length(telomereseq2)) {
  changeseqD <- append(changeseqD, if(telomereseq2[i] == "-") refseq2[i] else "." )
  changeseqI <- append(changeseqI, if(refseq2[i] == "-") telomereseq2[i] else "." )
}

# Group individual bases into connected mutations. ----

biglistD <- c()
as(biglistD,"list")
biglistDloc <- c()
biglistI <- c()
as(biglistI,"list")
biglistIloc <- c()

for (i in 2:length(changeseqD)) {
  if(changeseqD[i] == ".") 
    if(changeseqD[i-1] != ".") (biglistD <- append(biglistD,str_flatten(tempA)))
  if(changeseqD[i] == ".") 
    if(changeseqD[i-1] != ".") (biglistDloc <- append(biglistDloc,i))
  if(changeseqD[i] != ".")
    if(changeseqD[i-1] != ".") (tempA <- append(tempA,changeseqD[i])) else (tempA <- list(changeseqD[i]))
}

for (i in 2:length(changeseqI)) {
  if(changeseqI[i] == ".") 
    if(changeseqI[i-1] != ".") (biglistI <- append(biglistI,str_flatten(tempA)))
  if(changeseqI[i] == ".") 
    if(changeseqI[i-1] != ".") (biglistIloc <- append(biglistIloc,i))
  if(changeseqI[i] != ".")
    if(changeseqI[i-1] != ".") (tempA <- append(tempA,changeseqI[i])) else (tempA <- list(changeseqI[i]))
}

# You've generated lists of the 5' position and identity of each deletion (D) and insertion (I)

biglistD
biglistDloc
biglistI
biglistIloc

# A problem I have not yet worked out is that you will have a bunch of mismatches at the beginning or end of the sequence.
# This is due to the reference and actual sequences not being the exact same size and starting at the same base.
# This will show up in the final table as a bunch of perfect telomere repeats.
# You can check this by eye by pulling up the alignments. 

alignmentR
alignmentV

# Make this all into a big list! ----
# Combined table ordered by location.
# This is in "tidy" data format, which means that each row is one observation and each column is one variable.
# This allows us to more easily display the data using ggplot2.

bigtableD <- data.frame(biglistD)
bigtableD <- setNames(bigtableD,c("bases"))
bigtableD$condition <- "deletion"
bigtableD$location <- biglistDloc
bigtableI <- data.frame(biglistI)
bigtableI <- setNames(bigtableI,c("bases"))
bigtableI$condition <- "insertion"
bigtableI$location <- biglistIloc
overview <- rbind(bigtableD,bigtableI)
overview$background <- telomerename
overview[order(overview$location),]
overview

# You can create pretty tables for export, though they do not include location.

tbl_summary(bigtableD)
tbl_summary(bigtableI)

# Create nice looking plot showing locations, color coded insertions vs deletions, and base pair changes. ----
# You can find the tutorial used in the creation of this graph at https://rafalab.github.io/dsbook/ggplot2.html

p <- overview %>% ggplot(aes(x = location, y = background))
p +
  geom_point(aes(col = condition), size = 1) +
  geom_text_repel(aes(label = bases, col = condition), nudge_y = 0, max.overlaps = 20) +
  ylab("") +
  xlab("Distance from subtelomere (bp)") +
  ggtitle(telomerename) +
  scale_color_discrete(name = "")

# Create a text file that has the info for the current telomere.
# Ideally, we could generate a bunch of these, then use them to compare across telomeres or extract metrics

write_tsv(overview, paste("table", telomerename, ".txt", sep = ""))

beep()

stop()

# Once text files have been generated for multiple telomeres, they can be plotted together ----
# You'll need to customize the filepaths, change the table names depending on the number, etc.
# Could be more generalized

tableA <- read_delim("text_outputs/tableWS274_chrIL.txt", delim = "\t")
tableB <- read_delim("text_outputs/tableWS274_chrIIL.txt", delim = "\t")
tableC <- read_delim("text_outputs/tableWS274_chrIIIL.txt", delim = "\t")
tableE <- read_delim("text_outputs/tableWS274_chrVL.txt", delim = "\t")
tableF <- read_delim("text_outputs/tableWS274_chrXL.txt", delim = "\t")
tableH <- read_delim("text_outputs/tableWS274_chrIIR.txt", delim = "\t")
tableJ <- read_delim("text_outputs/tableWS274_chrIVR.txt", delim = "\t")
tableL <- read_delim("text_outputs/tableWS274_chrXR.txt", delim = "\t")

overview2 <- rbind(tableA, tableB, tableC, tableD, tableE, tableF, tableG, tableH, tableI, tableJ, tableK, tableL)

p <- overview2 %>% ggplot(aes(x = location, y = background))
p +
  geom_point(aes(col = condition), size = 2, alpha = .75, pch = 15) +
  ylab("") +
  xlab("Distance from subtelomere (bp)") +
  scale_color_discrete(name = "")

beep()