#!/bin/bash

NPROC=40


NAMES=()
NAMES+=( 'KESAMOCD4_S1' )
NAMES+=( 'KESAMOCD8_S2' )
NAMES+=( 'KESASangCD4N_S3' )
NAMES+=( 'KESAMOCD4M_S4' )
NAMES+=( 'KESASangCD8N_S5')
NAMES+=( 'KESASangCD8M_S6' )
NAMES+=( 'DS4-SANG-CD4NAIVE_S7' )
NAMES+=( 'DS4-SANG-CD4MEMOIRE_S8' )

DIR_PROJ=/data2/tmp/P17019/fastq/

DIR_OUT=star_out3


STAR_INDEX=/data2/fdb/star/Homo_sapiens/UCSC/hg19/res
FILE_GTF=/data2/fdb/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf


module load star/2.5.2b

mkdir -p $DIR_OUT

for NAME in ${NAMES[@]}
do
    echo $NAME

    R1=${NAME}_R1_001.fastq.gz
    R2=${NAME}_R2_001.fastq.gz

    echo $R1
    echo $R2

STAR                                        \
    --genomeDir $STAR_INDEX                 \
    --runThreadN $NPROC                     \
    --readFilesIn                           \
    ${DIR_PROJ}/${R1}                       \
    ${DIR_PROJ}/${R2}                       \
    --readFilesCommand zcat                 \
    --runMode alignReads                    \
    --outSAMtype BAM SortedByCoordinate     \
    --outSAMunmapped Within                 \
    --quantMode GeneCounts TranscriptomeSAM \
    --sjdbGTFfile ${FILE_GTF}               \
    --outFileNamePrefix "${DIR_OUT}/${NAME}"

done

#   --outBAMsortingThreadN $NPROC           \

#--genomeLoad LoadAndKeep                \

# --twopassMode Basic                     \

# --outFilterScoreMinOverLread 0.3
# --outFilterMatchNminOverLread 0.3

# --outFilterScoreMinOverLread 0
# --outFilterMatchNminOverLread 0
# --outFilterMatchNmin 0
# # --outFilterMismatchNmax 2
