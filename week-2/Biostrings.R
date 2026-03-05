library(Biostrings)
browseVignettes("Biostrings")

# DNAstring
dna1 <- DNAString("ACGCGGTCT")
dna1

# subset DNAstring
dna1[2:4]

# DNAStingSet
dna2 <- DNAStringSet(c("AGT", "AGG", "AGC"))
dna2

# subset DNAStringSet
dna2[1:2]
dna2[[1]]

# names DNA set
names(dna2) = paste0("seq", 1:length(dna2))
dna2

# IUPAC_CODE_MAP
IUPAC_CODE_MAP

# function for DNAstring
rev(dna1)
reverse(dna1)
reverseComplement(dna1)
translate(dna1)
alphabetFrequency(dna1)
# GC contecnt
letterFrequency(dna1, letters = "GC")
# repeat neuclotide frequency
dinucleotideFrequency(dna1)




# function for DNAStringset
width(dna2)
sort(dna2)
rev(dna2)
reverse(dna2)
reverseComplement(dna2)
translate(dna2)
alphabetFrequency(dna2)
# GC contecnt
letterFrequency(dna2, letters = "GC")
# repeat neuclotide frequency
dinucleotideFrequency(dna2)
consensusMatrix(dna2)
