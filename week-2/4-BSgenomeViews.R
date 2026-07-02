library(BSgenome)
library(BSgenome.Scerevisiae.UCSC.sacCer2)
library(AnnotationHub)


dnaseq <- DNAString("ACGTACGT")
dnaseq

seqnames(Scerevisiae)

# biostring matchPattern single 
vi <- matchPattern(dnaseq, Scerevisiae[["chrI"]])
vi

ranges(vi)
Scerevisiae$chrI[57932:57939]
alphabetFrequency(vi)
shift(vi, 10)

# Biostring vmatchPattern dnaseq within entire genome
gr <- vmatchPattern(dnaseq, Scerevisiae)
gr

## view the basis too
vi2 <- Views(Scerevisiae, gr)
vi2

ahub <- AnnotationHub()
ahub

qh <- query(ahub, c("sacCer2", "genes"))
qh

genes <- ahub[["AH7048"]]
genes

prom <- promoters(genes)
prom

prom <- trim(prom)
prom

seqinfo(prom) <- seqinfo(Scerevisiae)[seqlevels(prom)]
promViews <- Views(Scerevisiae, prom)
promViews

gcProm <- letterFrequency(promViews, "GC", as.prob = TRUE)
gcProm

plot(density(gcProm))
abline(v = 0.38)
