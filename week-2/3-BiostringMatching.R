library(Biostrings)
library(BSgenome)
library(BSgenome.Scerevisiae.UCSC.sacCer2)



# Biostring-Matching-Strings: 
# string to string
dnaseq <- DNAString("ACGTACGT")
dnaseq
seqnames(Scerevisiae)
matchPattern(dnaseq, Scerevisiae[["chrI"]])
countPattern(dnaseq, Scerevisiae[["chrI"]])

# set of string to set of strings
vmatchPattern(dnaseq, Scerevisiae)
dnaseq == reverseComplement(dnaseq)

# pairwiseALignment
dnaseq2 <- DNAString("ACGTTGGT")
pwalign::pairwiseAlignment(dnaseq, dnaseq2)

# Create a Position Weight Matrix (PWM) from multiple sequences.
motifs <- DNAStringSet(c("ACGT", "ACGA", "ACGG"))
pwm <- PWM(motifs)

matchPWM(pwm, Scerevisiae$chrI)

# string to set of strings
# set of string to string






