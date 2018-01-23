
source("http://bioconductor.org/biocLite.R")


pkgs <- c(
    "Rsamtools"
  , "TxDb.Hsapiens.UCSC.hg19.knownGene"
  , "GenomicFeatures"
  , "GenomicAlignments"
  , "BiocParallel"
  , "edgeR"
  , "limma"
  , "cluster"
  , "foreign"
  , "nnet"
  , "rpart"
  , "Glimma" # Manual Depend: foreign nnet rpart
  , "gplots"
  , "RColorBrewer"
  , "org.Hs.eg.db"
)

for (pkg in pkgs) {
    print(pkg)
    biocLite(pkg)
}
