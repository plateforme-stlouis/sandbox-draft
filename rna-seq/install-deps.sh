#!/bin/bash

pkgs=()
pkgs+=( "Rsamtools" )
pkgs+=( "TxDb.Hsapiens.UCSC.hg19.knownGene" )
pkgs+=( "GenomicFeatures" )
pkgs+=( "GenomicAlignments" )
pkgs+=( "BiocParallel" )
pkgs+=( "edgeR" )
pkgs+=( "limma" )
pkgs+=( "foreign" )
pkgs+=( "nnet" )
pkgs+=( "rpart" )
pkgs+=( "Glimma" ) # manual Depend: foreign nnet rpart
pkgs+=( "gplots" )
pkgs+=( "RColorBrewer" )
pkgs+=( "org.Hs.eg.db" )

curl -O http://bioconductor.org/biocLite.R
for pkg in ${pkgs[@]}
do
    echo $pkg
    Rscript --no-init-file -e "source(\"biocLite.R\") ; biocLite(\"$pkg\")"
done

[ -f biocLite.R ] && rm biocLite.R
