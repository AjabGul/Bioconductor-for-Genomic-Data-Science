library(GenomicRanges)

# Run length Encoding
rl <- Rle(c(1,1,1,1,1,1,2,2,2,2,2,4,4,2))
rl

# display lengths and values of the vector
runLength(rl)
runValue(rl)

# convert the Rle back to numeric vector
as.numeric(rl)


ir <- IRanges(start = c(2,8), width = 4)
ir

# aggregate
aggregate(rl, ir, FUN = mean)

# covrage
ir2 <- IRanges(start = 1:5, width = 3)
ir2

coverage(ir2)

# slice
slice(rl, 2)
slice(rl, 3)

# view
vi <- Views(rl, IRanges(2,8))
vi
vi <- Views(rl, IRanges(c(2,8),  width = 2))
vi
mean(vi)

# GRnageView
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:10, width = 3))
gr

rl <- coverage(gr)
rl


vi <- Views(rl, as(GRanges("chr1", ranges = IRanges(3,7)), "RangesList"))
vi

vi$chr1
