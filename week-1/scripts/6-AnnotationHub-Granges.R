# Open the documentation (vignettes) for the AnnotationHub package
browseVignettes("AnnotationHub")


# Load the AnnotationHub package to access genomic datasets
library(AnnotationHub)
# Load GenomicRanges to work with genomic interval objects (GRanges)
library(GenomicRanges)


# Connect to the online AnnotationHub resource database
ahub <- AnnotationHub()
# Display the AnnotationHub object summary
ahub
# Count total number of resources available in AnnotationHub
length(ahub)


# List all unique data providers in AnnotationHub (e.g., UCSC, ENCODE, etc.)
unique(ahub$dataprovider)
# List all species available in AnnotationHub
unique(ahub$species)
# List all different R object classes stored (e.g., GRanges, ChainFile, etc.)
unique(ahub$rdataclass)


# Search AnnotationHub for human ChainFiles from UCSC (note: typo in 'Homo sapien')
dm <- query(ahub, c("ChainFile", "UCSC", "Homo sapien"))

# Display the search results for chain files
dm

# Extract metadata columns from the search results
df <- mcols(dm)

# Display metadata information about the returned datasets
df

# Restrict AnnotationHub to only Homo sapiens datasets
ahub <- subset(ahub, species == "Homo sapiens")

# Search for H3K4me3 ChIP-seq datasets from the GM12878 cell line
qhs <- query(ahub, c("H3K4me3", "Gm12878"))

# Display the query results for H3K4me3 datasets
qhs

# Load the second dataset from the search results into gr1 (as a GRanges object)
gr1 <- qhs[[2]]

# Load rtracklayer package for importing/exporting genomic data formats
require("rtracklayer")

# Load the fourth dataset from results into gr2
gr2 <- qhs[[4]]

# Display the genomic ranges object stored in gr1
gr1

# Summarize the distribution of peak widths in gr1
summary(width(gr1))

# Display the genomic ranges object stored in gr2
gr2

# Summarize the distribution of peak widths in gr2
summary(width(gr2))

# Create a frequency table of peak widths in gr2
table(width(gr2))

# Assign gr2 to a new object called peaks (for downstream overlap analysis)
peaks = gr2

# Display the H3K4me3 query results again
qhs

# Display only the 4th dataset entry from query results
qhs[4]



# ref seq ----

# Search AnnotationHub for datasets related to RefSeq annotations
qhs <- query(ahub, "RefSeq")

# Display the search results returned from AnnotationHub
qhs

# Check which genome build (e.g., hg19, hg38) each dataset belongs to
qhs$genome

# Load the first RefSeq dataset from the search results into an object called 'genes'
genes <- qhs[[1]]

# Display the gene annotation object (usually a GRanges object)
genes

# Count how many times each gene name appears in the dataset
table(genes$name)

# Count how many genes have 1 transcript, 2 transcripts, etc. (frequency of frequencies)
table(table(genes$name))

# Display the input arguments and default parameter values of the promoters() function
args(promoters)

# Create promoter regions (e.g., 2000 bp upstream, 200 bp downstream)
prom <- promoters(genes, upstream=2000, downstream=200)
prom
# Check the distribution of promoter widths (lengths of promoter regions)
table(width(prom))

ov <- findOverlaps(prom, peaks)
# Find genomic overlaps between promoter regions (prom) and peaks

ov
# Display the overlap object (shows which promoters overlap which peaks)

length(unique(queryHits(ov)))
# Count number of unique promoters that overlap at least one peak

length(unique(subjectHits(ov)))
# Count number of unique peaks that overlap at least one promoter

length(subsetByOverlaps(peaks, prom, ignore.strand = TRUE )) / length(peaks)
# Proportion of peaks that overlap promoters (strand ignored)

length(subsetByOverlaps(prom, peaks, ignore.strand = TRUE )) / length(prom)
# Proportion of promoters that overlap peaks (strand ignored)

sum(width(reduce(peaks, ignore.strand = TRUE))) 
# Total genomic base pairs covered by peaks (after merging overlapping peaks)

sum(width(reduce(peaks, ignore.strand = TRUE))) / 10^6
# Convert total peak coverage into megabases (Mb)

sum(width(reduce(prom, ignore.strand = TRUE))) / 10^6
# Total promoter coverage in megabases (Mb)

sum(width(intersect(peaks, prom, ignore.strand = TRUE))) / 106
# (Likely typo) Calculates overlapping base pairs but divides by 106 instead of 10^6

sum(width(intersect(peaks, prom, ignore.strand = TRUE))) / 10^6
# Total overlapping base pairs between peaks and promoters in Mb

inOut = matrix(0, ncol = 2, nrow = 2)
# Create 2x2 matrix initialized with zeros (for contingency table)

colnames(inOut) = c("in", "out")
# Name columns as inside vs outside promoters

rownames(inOut) = c("in", "out")
# Name rows as peaks inside vs outside promoters

inOut
# Display empty contingency table

inOut[1,1] = sum(width(intersect(peaks, prom, ignore.strand = TRUE)))
# Fill cell: base pairs of peaks overlapping promoters

inOut[1,2] = sum(width(setdiff(peaks, prom, ignore.strand = TRUE)))
# Fill cell: base pairs of peaks NOT overlapping promoters

inOut[2,1] = sum(width(setdiff(prom, peaks, ignore.strand = TRUE)))
# Fill cell: base pairs of promoters NOT overlapping peaks

inOut
# Display partially filled contingency table

colSums(inOut)
# Show column totals (total in vs out)

rowSums(inOut)
# Show row totals (total peaks vs promoters)

inOut[2,2] = 3*10^9 - sum(inOut)
# Fill remaining genome space assuming human genome size ≈ 3 billion bp

inOut
# Display full contingency table

fisher.test(inOut)$statistic
# Perform Fisher's exact test to measure enrichment association

oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2])
# Manually calculate odds ratio from contingency table

oddsRatio
# Print odds ratio

# half the human genome 
inOut[2,2] = 0
# Reset background cell

inOut[2,2] = 1.5 * 10^9 - sum(inOut)
# Now assume genome size is 1.5 billion bp (half genome scenario)

inOut
# Display updated contingency table

oddsRatio = inOut[1,1] * inOut[2,2] / (inOut[2,1] * inOut[1,2])
# Recalculate odds ratio under half-genome assumption

oddsRatio
# Print updated odds ratio
