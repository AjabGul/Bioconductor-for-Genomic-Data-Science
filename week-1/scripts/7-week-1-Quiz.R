library(AnnotationHub)
library(GenomicRanges)

ah <- AnnotationHub()

# Q1- ----
# Use the AnnotationHub package to obtain data on “CpG Islands” in the human genome.

# search specie homo sapiens
query(ah, "Homo sapiens")

# search CpG avalaible data in homo sapiens
query(ah, c("CpG", "Homo sapiens"))

# extract the exact CpG island data in homo sapiens
query_results <- query(ah, c("CpG Islands", "Homo sapiens"))
query_results

# mcols (Metadata Columns) extracts extra information (annotations) attached to an object.
#       display only three column title, genome, rdataDateAdded
mcols(query_results, c("title", "genome", "rdatadateadded"))

# check either cpg hg38 data is available
query(ah, c("CpG", "hg38"))

# Load CpG island data from the above 4 set of data
cpg_islands <- ah[["AH5086"]]
cpg_islands

seqnames(cpg_islands)
ranges(cpg_islands)
strand(cpg_islands)

# Inspect first 10 CpG islands
head(cpg_islands, 10)

# Check total number
length(cpg_islands)

# Check chromosome distribution
table(seqnames(cpg_islands))


# Autosomes in human genome
# paste0() concatenates strings and numbers without a separator.
autosomes <- paste0("chr", 1:22)
autosomes

# Keep only CpG islands on chr1–chr22
autosomal_cpg <- cpg_islands[seqnames(cpg_islands) %in% autosomes]
autosomal_cpg
# Total number of CpG islands on autosomes
length(autosomal_cpg)


# Q2: CpG islands on chromosome 4 ----
chr4_cpg <- cpg_islands[seqnames(cpg_islands) == "chr4"]
chr4_cpg
# Total number of CpG islands on chr4
length(chr4_cpg)


# Q3: Obtain the data for the H3K4me3 histone modification for the H1 cell line ----
#     from Epigenomics Roadmap, using AnnotationHub. Subset these regions
#     to only keep regions mapped to the autosomes (chromosomes 1 to 22).
#     Question: How many bases does these regions cover?

# step: 0 extract the H1 H3K4me3 cell line
query_H1 <- query(ah, c("H3K4me3", "E003", "Homo sapiens"))
query_H1

query_H1[grep("H3K4me3", query_H1$title)]

h1_cellline <- ah[["AH29884"]]
h1_cellline  

# Step 1: Define autosomes
autosomes <- paste0("chr", 1:22)

# Step 2: Subset H3K4me3 peaks to autosomes
auto_h1 <- h1_cellline[ seqnames(h1_cellline) %in% autosomes ]
auto_h1

# Step 3: Merge overlapping regions
# reduce() merges overlapping genomic regions into continuous regions.
auto_h1_reduced <- reduce(auto_h1)
auto_h1

# Step 4: Calculate total bases 
# Width = number of base pairs covered by that genomic region.
total_bases <- sum(width(auto_h1_reduced))
total_bases


# Q4: Obtain the data for the H3K27me3 histone modification for the H1 ----
# cell line from Epigenomics Roadmap, using the AnnotationHub package. Subset these 
#    regions to only keep regions mapped to the autosomes. In the return data, 
#    each region has an associated "signalValue". Question: What is the mean 
#    signalValue across all regions on the standard chromosomes?

# query H3K27me3 peaks for H1 (E003)
query_H3K27 <- query(ah, c("H3K27me3", "E003", "Homo sapiens"))
query_H3K27

H3K27me3_data <- ah[["AH29892"]]
H3K27me3_data

autosomal_h3k27 <- paste0("chr", 1:22)
autosomal_h3k27

h1_h3k27_auto <- H3K27me3_data[seqnames(H3K27me3_data) %in% autosomal_h3k27]
h1_h3k27_auto

# Note: mcols only extract metaData Columns not Grnages of genomic positions
mcols(h1_h3k27_auto)
mcols(h1_h3k27_auto)$signalValue
mean(mcols(h1_h3k27_auto)$signalValue)

mean_signal <- mean(mcols(h1_h3k27_auto)$signalValue)
mean_signal



# Question 5 ----
# Bivalent regions are bound by both H3K4me3 and H3K27me3.
# Question: Using the regions we have obtained above, how many bases on 
# the standard chromosomes are bivalently marked?

# A bivalent region is a DNA region that has two different histone modifications at the same time:
# 🟢 H3K4me3 → activating mark
# 🔴 H3K27me3 → repressive mark
# When both marks are present at the same genomic region, it is called Bivalent chromatin

# install.packages("magick")
library(magick)
img <- image_read("/home/ajab/AAjabGul/CS-Data/1-DataScience/2-BI_and_GDS/GDS-Coursera/5-Bioconductor_R/week-1/scripts/bivalent_mark.png")
print(img)
  
# Q6: We will examine the extent to which bivalent regions overlap CpG Islands. ----

# Question 7 ----
# Question: How big a fraction (expressed as a number between 0 and 1) of 
# the bases which are part of CpG Islands, are also bivalent marked.


# Q8: ----
# Question: How many bases are bivalently marked within 10kb of CpG Islands?
# Tip: consider using the "resize()"" function.


# Question 9 ----
# Question: How big a fraction (expressed as a number between 0 and 1) of the human genome is contained in a CpG Island?
# Tip 1: the object returned by AnnotationHub contains "seqlengths".
# Tip 2: you may encounter an integer overflow. As described in the session on R Basic Types, you can address this by converting integers to numeric before summing them, "as.numeric()".


# Question 10 ----
# Question: Compute an odds-ratio for the overlap of bivalent marks with CpG islands.

