library(GenomicRanges)

# basics of GRange ----
gr <- GRanges(
  seqnames = c("chr1"), 
  ranges = IRanges(start = c(1,3,5), width = 3),
  strand = c("+", "-", "+")
  )
gr


# flank():
# + strand → flank is before your region (upstream).
# - strand → flank is after your region (downstream).
flank(gr, 5)


# promoters():
# should be applied to GRanges objects where each range represents 
# a gene (or transcript), and the range start (for + strand) 
# or end (for − strand) corresponds to the transcription start site (TSS).
promoters(gr)

# seqinfor tell about name of chr, length, circular or linear
seqinfo(gr)
seqlengths(gr)
seqlevels(gr)
seq_along(gr)
seq_len(gr)


gaps(gr)


gr1 <- GRanges(
  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
  strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score = 1:10,
  GC = seq(1, 0, length=10))
gr1

gaps(gr1)
seqlengths(gr1)
sort(gr1)
seqlevels(gr1)
seqnames(gr1)
gaps(gr1)
genome(gr1) = "hg19"
genome(gr1)
seqinfo(gr1)
ranges(gr1)
granges(gr1)


# findOverlape 
gr2 = gr1
genome(gr2) = "hg18"
gr2
findOverlaps(gr2,gr1)

# Dataframe and Granges ----

# DataFrame:
# All columns must have same length
# ✔ Column names are required
# ✔ Rows can be named (optional)
# ✔ Columns can be numeric, character, logical, factor, etc.


# creating DataFrame using IRanges
ir <- IRanges(start = 1:3, width = 2)
ir

df <- DataFrame(ir = ir, score = rnorm(3))
df

df[1:1]

df$ir

df$score

df2 <- data.frame(ir = ir)
df2


# GRanges and DataFrame

gr3 <- GRanges(
  seqnames = "chr1",
  strand = c("+", "-", "+"),
  ranges = IRanges(start = c(1,3,5), width = 3)
)

gr3

# the entire sets are called values or column metaData also called Data Frame
values(gr3) = DataFrame(score = rnorm(3))
gr3

values(gr3)

mcols(gr3)


# findOverlapes and subsetbyvoerlapes in GRanges ----
library(GenomicRanges)

gr4 <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(start = c(1, 10, 20), width = 5),
  strand = "+"
)

gr4

gr5 <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(start = c(3, 12), width = 5),
  strand = "+"
)
gr5

findOverlaps(gr4, gr5)

subsetByOverlaps(gr4, gr5)

# Making GRanges from dataframe ----

df3 <- DataFrame(seqnames = "chr1", start = 1:3, end = 3:5, score = rnorm(3))
df3

makeGRangesFromDataFrame(df3)
makeGRangesFromDataFrame(df3, keep.extra.columns = TRUE)
