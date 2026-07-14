library(GenomicFeatures)
#BiocManager::install("GenomicFeatures")
library(GenomicRanges)


# load transcript database for human hg19
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
# rename for convenience
txdb <- tolower(txdb)
txdb

#extracting genomic features fron the txdb ----

# get the genes ranges from txdb
genes <- genes(txdb)
genes

# extract promoters
prom <- promoters(txdb)
prom

# extract transcript ranges from txdb
transcripts <- transcripts(txdb)
transcripts

# extract exones ranges
exons <- exons(txdb)
exons
# exons grouped by transcripts
exonByText <- exonsBy(txdb, by="tx")
exonByText


# extract cds
cds <- cds(txdb)
cds
# coding seq cds grouped by transcripts
cdsByTx <- cdsBy(txdb, by="tx")
cdsByTx

# extract intronsByTranscript()
intronByTranscript <- intronsByTranscript(txdb)
intronByTranscript

# subset genes, transcripts and exon by overlaps:  ----
# finds all genes that overlaps on chr1:100000-200000 region
subset_gene <- subsetByOverlaps(genes, GRanges("chr1", IRanges(start =100000, end = 200000)))
subset_gene
# see how many genes overlaps
length(subset_gene)

# find all transcript that overlaps on chr1:100000:200000
subset_transcripts <- subsetByOverlaps(exonByText, GRanges("chr1", IRanges(start =100000, end = 200000)))
subset_transcripts
length(subset_transcripts)

# find all exon overlaps on chr1:100k-200k
subset_exons <- subsetByOverlaps(exons, GRanges("chr1", IRanges(start =100000, end = 200000)))
subset_exons
length(subset_exons)


# get the transcripts lengths including cds length
tx_length <- transcriptLengths(txdb, with.cds_len = TRUE)
tx_length

sum_width <- sum(width(cdsByTx[[2]]))
sum_width

# to make your own transcriptional data base using ----
#  UCSC, biomark, UCSC, ensemble, or GFF data 
# we use makeTxDB()

# BiocManager::install("txdbmaker") # installation
library(txdbmaker)

# makeTxDbFromBiomart()
# makeTxDbFromEnsembl()
# makeTxDbFromUCSC()
# makeTxDbFromGFF()