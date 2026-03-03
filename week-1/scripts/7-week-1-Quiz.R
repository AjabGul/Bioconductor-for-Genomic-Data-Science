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
auto_h1_reduced

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

# A bivalent region is a DNA region that has two different histone 
# modifications at the same time:
# 🟢 H3K4me3 → activating mark
# 🔴 H3K27me3 → repressive mark
# When both marks are present at the same genomic region, 
# it is called Bivalent chromatin

# we already have autocomal reduced of h3k4me3
auto_h1_reduced
# We reduce to merge overlapping regions (important before intersection).
h3k27me3_reduced <- reduce(h1_h3k27_auto)
h3k27me3_reduced

# bivalent marked regions
# intersect() returns the genomic coordinates that are 
# shared (overlapping) between two sets of genomic regions.
bivalent_region <- intersect(auto_h1_reduced, h3k27me3_reduced)
bivalent_region

# calculate the bases on the standard chromosomes are bivalent marked?
bivalent_total_bases <- sum(width(bivalent_region))
bivalent_total_bases

# Q6: We will examine the extent to which bivalent regions overlap CpG Islands. ----
# Question: how big a fraction (expressed as a number between 0 and 1) of 
# the bivalent regions, overlap one or more CpG Islands?

# bivalent_region: GRanges object containing genomic regions
# that are marked by BOTH H3K4me3 and H3K27me3 (bivalent chromatin)
bivalent_region

# autosomal_cpg: GRanges object containing CpG Islands
# located only on chromosomes 1–22 (autosomes)
autosomal_cpg


# findOverlaps() identifies genomic regions in 'bivalent_region'
# that overlap with any region in 'autosomal_cpg'
# It returns a Hits object showing which regions overlap
hits <- findOverlaps(bivalent_region, autosomal_cpg)
# Display the overlap mapping between query (bivalent) 
# and subject (CpG islands)
hits


# queryHits(hits) extracts indices of bivalent regions
# that have at least one overlap with a CpG island
# unique() ensures each overlapping bivalent region
# is counted only once (even if overlapping multiple CpGs)
bivalent_with_cpg <- length(unique(queryHits(hits)))
# Print number of bivalent regions overlapping CpG islands
bivalent_with_cpg


# length() gives total number of bivalent regions
# regardless of whether they overlap CpG islands or not
total_bivalent <- length(bivalent_region)
# Print total number of bivalent regions
total_bivalent


# Compute fraction of bivalent regions that overlap
# one or more CpG Islands
# Result is a number between 0 and 1
fraction_overlap <- bivalent_with_cpg / total_bivalent
# Print final fraction
fraction_overlap

# Question 7 ----
# Question: How big a fraction (expressed as a number between 0 and 1) of 
# the bases which are part of CpG Islands, are also bivalent marked.

# Merge overlapping CpG Islands on autosomes into continuous regions
# This avoids double counting bases in overlapping CpG Islands
cpg_reduced <- reduce(autosomal_cpg)
# Display the reduced CpG Island regions
cpg_reduced


# Merge overlapping bivalent regions into continuous genomic ranges
# Ensures each base is counted only once
bivalent_reduced <- reduce(bivalent_region)
# Display the reduced bivalent regions
bivalent_reduced


# Find genomic regions that are both CpG Islands and bivalently marked
# intersect() returns only the overlapping portions of the two GRanges objects
cpg_bivalent_overlap <- intersect(cpg_reduced, bivalent_reduced)
# Display the overlapping regions (CpG + bivalent)
cpg_bivalent_overlap


# Count the total number of bases in the overlapping regions
# width() gives the length of each region; sum() adds them together
overlap_bases <- sum(width(cpg_bivalent_overlap))
# Print number of bases in CpG Islands that are also bivalent
overlap_bases


# Count total number of bases in all CpG Islands
total_cpg_bases <- sum(width(cpg_reduced))
# Print total CpG bases on autosomes
total_cpg_bases


# Compute fraction of CpG bases that are also bivalently marked
# This is a number between 0 and 1
fraction_cpg_bivalent <- overlap_bases / total_cpg_bases
# Print final fraction
fraction_cpg_bivalent


# Q8: ----
# Question: How many bases are bivalently marked within 10kb of CpG Islands?
# Tip: consider using the "resize()"" function.

# Resize each CpG Island to 10 kb upstream and downstream
# width = 10,000 * 2 + original width? Actually we use "flank"
# But easiest: extend by 10kb on both sides
cpg_flank_10kb <- resize(cpg_reduced, width = width(cpg_reduced) + 20000, fix = "center")
cpg_flank_10kb

# Find bases in bivalent regions that fall within the ±10 kb of CpG Islands
bivalent_near_cpg <- intersect(bivalent_reduced, cpg_flank_10kb)
# Display overlapping regions
bivalent_near_cpg


# Answer: Total bases in bivalent regions within 10 kb of CpG Islands
bivalent_near_cpg_bases <- sum(width(bivalent_near_cpg))
bivalent_near_cpg_bases


# Optional Step: Fraction of total bivalent bases that are within 10 kb of CpG Islands
fraction_near_cpg <- bivalent_near_cpg_bases / sum(width(bivalent_reduced))
fraction_near_cpg

# Question 9 ----
# Question: How big a fraction (expressed as a number between 0 and 1) of 
# the human genome is contained in a CpG Island?
# Tip 1: the object returned by AnnotationHub contains "seqlengths".
# Tip 2: you may encounter an integer overflow. As described in the 
# session on R Basic Types, you can address this by converting integers 
# to numeric before summing them, "as.numeric()".


# Fraction of human genome in CpG Islands
# (Using AUTOSOMES ONLY: chr1–22)

# we already have autosomal cpg
autosomal_cpg

# Reduce overlapping CpG Islands
autosomal_cpg_reduced <- reduce(autosomal_cpg)

# Calculate total CpG bases (convert to numeric)
total_cpg_bases <- sum(as.numeric(width(autosomal_cpg_reduced)))

# Load genome object
genome <- ah[["AH5018"]]

genome_size <- sum(as.numeric(seqlengths(genome)[autosomes]))

# Compute fraction
fraction_genome_cpg <- total_cpg_bases / genome_size

# Print result
fraction_genome_cpg

# Question 10 ----
# Question: Compute an odds-ratio for the overlap of bivalent marks with CpG islands.

# Odds Ratio for overlap of
# bivalent marks with CpG islands
# (Base-level enrichment, autosomes only)


# Make sure regions are reduced
bivalent_reduced <- reduce(bivalent_region)
cpg_reduced <- reduce(autosomal_cpg)

# Compute overlap bases (a)
overlap <- intersect(bivalent_reduced, cpg_reduced)
a <- sum(as.numeric(width(overlap)))

# Compute totals
total_cpg <- sum(as.numeric(width(cpg_reduced)))
total_bivalent <- sum(as.numeric(width(bivalent_reduced)))

autosomes <- paste0("chr", 1:22)
genome <- ah[["AH5018"]]
genome_size <- sum(as.numeric(seqlengths(genome)[autosomes]))

# Remaining cells of 2x2 table
b <- total_bivalent - a         # Bivalent but NOT CpG
c <- total_cpg - a              # CpG but NOT bivalent
d <- genome_size - (a + b + c)  # Neither CpG nor bivalent

#  Compute Odds Ratio
odds_ratio <- (a * d) / (b * c)

odds_ratio
