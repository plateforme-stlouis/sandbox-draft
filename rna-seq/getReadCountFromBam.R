library(Rsamtools)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)

cwd = "star_out3/"
setwd(cwd)

files = list.files(path=".", pattern="^[A-Za-z1-90_-]+.sortedByCoord.out.bam")

#files  <-  "DS4-SANG-CD4MEMOIRE_S8Aligned.sortedByCoord.out.bam"
## files  <-  "DS4-SANG-CD4NAIVE_S7Aligned.sortedByCoord.out.bam"
## files  <-  "KESAMOCD4M_S4Aligned.sortedByCoord.out.bam"
## files  <-  "KESAMOCD4_S1Aligned.sortedByCoord.out.bam"
## files  <-  "KESAMOCD8_S2Aligned.sortedByCoord.out.bam"
## files  <-  "KESASangCD8N_S5Aligned.sortedByCoord.out.bam"
## files  <-  "KESASangCD4N_S3Aligned.sortedByCoord.out.bam"

files  <-  "KESASangCD8M_S6Aligned.sortedByCoord.out.bam"


print(files)
print("Files ok.")

bamfiles <- BamFileList(files)
print("Bam Files ok.")

ebg <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by="gene")
print("Exons By ok.")

## register(MulticoreParam(16))
## print("Register ok.")

print("Start...")
se <- summarizeOverlaps(features=ebg,
                        reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
print("Done.")

expr = assay(se)
save(expr, file="read.counts.RData")

setwd("..")
print(getwd())
print("All done.")
