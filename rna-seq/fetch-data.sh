#!/bin/bash

files=()
files+=( '2cells_1.fastq' )
files+=( '2cells_2.fastq' )
files+=( '6h_1.fastq' )
files+=( '6h_2.fastq' )
files+=( 'zebrafish-rna-seq.pdf' )
from=http://ftp.ebi.ac.uk/pub/training/Train_online/RNA-seq_exercise

DIR_DAT=data-ftp.ebi.ac.uk


mkdir -p $DIR_DAT
for file in ${files[@]}
do
    echo "Downloading $file..."
    curl -C - -o $DIR_DAT/$file $from/$file
done
echo "Done."
