library(GenomicRanges)

# tow GRanges object from chr1 and chr2
gr1 <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(1,5,10), width = 4))
gr1

gr2 <- GRanges(seqnames = "chr2", ranges = IRanges(start = c(2,8,15), width = 3))
gr2

# GRangeList: combine the the two GRnages in to GRangesList
grl <- GRangesList(gr1 = gr1, gr2 = gr2)
print(grl)


# first element of the GrangesList
grl[1]
grl[[1]]
# acced the element by name
grl["gr1"]


# get the start position and end position of each GRanges in the list
start(grl)
end(grl)
# get the seq name of the GRanges in the list
seqnames(grl)

# elementLegnths to get number element in each GRanges element
elementNROWS(grl)
length(grl)



# shift all ranges in a GRangesList by 10 bases
shifted_grl <- shift(grl, 10)
shifted_grl


# findOverlapes
# GRangeList and single RGrange Object
gr_single <- GRanges(seqnames = "chr1", ranges = IRanges(start = 3, width = 5))
gr_single
findOverlaps(grl, gr_single)

# findOverlapes between two GRnagesLists
findOverlaps(grl, grl)
