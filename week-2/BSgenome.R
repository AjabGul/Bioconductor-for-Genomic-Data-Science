# intallation
# BiocManager::install("BSgenome")
browseVignettes("BSgenome")
library(BSgenome)
library(Biostrings)

# available genomes in BSgenome data
available.genomes()
# load the geonome of
library(BSgenome.Scerevisiae.UCSC.sacCer2)
Scerevisiae


# list the names of all chromosome sequences in the yeast genome
seqnames(Scerevisiae)
# show the length (number of base pairs) of each chromosome
seqlengths(Scerevisiae)
# extract chromosome I sequence from the genome object
chrI <- Scerevisiae$chrI
chrI
# count GC-content of chrI
letterFrequency(chrI, "GC")
# as.prob = TRUE => proportion of GC bases relative to the chromosome length.
letterFrequency(chrI, "GC", as.prob = TRUE)



# create a BSParams object specifying the genome and 
# function to apply to each chromosome
param <- new("BSParams",
             X = Scerevisiae,
             FUN = letterFrequency)

# apply letterFrequency to every chromosome to count GC content
bsapply(param, "GC")
# convert the list of GC counts from all chromosomes into a single numeric vector
unlist(bsapply(param, "GC"))
# compute overall GC proportion of the whole genome 
# using total GC / total genome length
sum(unlist(bsapply(param, "GC"))) / sum(seqlengths(Scerevisiae))
# compute GC proportion separately for each chromosome
unlist(bsapply(param, "GC",  as.prob =TRUE))
