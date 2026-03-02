# G Range definition  ----

# it is an alternative to bedtools of command line in R

# A special data structure for storage of genomic intervels from 
# the GenomicRanges package (Bioconductor).

# It represents genomic intervals/ranges 
# (like chromosome regions, start–end positions, and strand).
# Think of it as a smart data frame for genomic coordinates. 


# IRange  ----
# is a vector that contains integers intervals 

library(IRanges)

# we need to provide 2 out of 3 start, end, width arguments

# width = end - start + 1
ir1 <- IRanges(start = c(1,2,3), end = c(3,5,7)) 
ir1

# end = start + width -1
ir2 <- IRanges(start = c(1,2,3), width = 3)
ir2

# start = end - width + 1
ir3 <- IRanges(end = c(3,7,8), width = 3)
ir3

# get the rang from IRange matrix
start(ir1)
end(ir1)
width(ir1)
dim(ir1)
length(ir1)

# concatinate two irange
c(ir1, ir2)

# function for plot the i range ----
plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col = "black", sep = 0.5, ...) {
  height <- 1
  if (is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins) * (height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x) - 0.5, ybottom,
       end(x) + 0.5, ybottom + height,
       col = col, ...)
  title(main)
  axis(1)
}

par(mfrow = c(2,1))   # single plot
# par(mar = c(4,4,2,2)) # set margins: bottom, left, top, right



ir <- IRanges(start = c(1,3,7,9), end = c(4,4,8,10))
ir

# plotRanges() function:
# Plots the original intervals
# Each range is drawn separately
# Overlaps are shown as overlapping bars
plotRanges(ir)

# reduce() function:
# merges overlapping or adjacent ranges into the smallest 
# possible non-overlapping ranges

# application in genomics:
# 1-Total exon length
# 2-Union of peaks (ChIP-seq)
# 3-Covered genome regions
# 4-CNV merging
plotRanges(reduce(ir))

# disjoin() function: 
# Break everything into non-overlapping basic units.
# splits ranges into the smallest non-overlapping pieces
# Instead of merging, it cuts overlaps into atomic intervals.

# application in genomics:
# 1-Counting reads per exon part
# 2-Avoiding double-counting overlaps
# 3-Differential exon usage
plotRanges(disjoin(ir))

# -------------Manipulation of IRange------------------------------- ----

# -------------------resize()

resize(ir, width = 1, fix = "start")

resize(ir, width = 1, fix= "center")

# union and intersection of the ranges

ir5 <- IRanges(start = c(1,3,5), width = 1)
ir6 <- IRanges(start = c(4,5,6), width = 1)
ir5
ir6

# ------------UNION:
# Merges all ranges, combining overlaps
# Produces non-overlapping ranges that cover all positions in either set
union(ir5, ir6)
# union followed by reduce
reduce(c(ir5,ir6))

# -----------INTERSECTION:
# Which regions are common to both sets
intersect(ir5, ir6)

# -----------findOverlape() function
# it will show the postion of the over lap bases in query and subject
ov <- findOverlaps(ir5,ir6)
ov
# hit of query
queryHits(ov)
unique(queryHits(ov))
# hit of subject
subjectHits(ov)
unique(subjectHits(ov))

#------- countOverlap()
countOverlaps(ir5,ir6)

# how to show the arguments/parameters of a function
args(findOverlaps)

# ------nearest(snps, genes)
nearest(ir5,ir6)

