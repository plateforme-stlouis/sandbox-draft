#!/usr/bin/bash


mkdir -p subdir1
touch subdir1/G630_RH_C0006FZ_7_1_HTTLKBBXX.IND12.fastq
touch subdir1/G630_RH_C0006FZ_7_2_HTTLKBBXX.IND12.fastq
touch subdir1/G630_RH_C0006G5_7_1_HTTLKBBXX.IND19.fastq
touch subdir1/G630_RH_C0006G5_7_2_HTTLKBBXX.IND19.fastq

mkdir -p subdir2
touch subdir2/G630_RH_C0006G6_8_1_HTTLKBBXX.IND5.fastq
touch subdir2/G630_RH_C0006G6_8_2_HTTLKBBXX.IND5.fastq
touch subdir2/G630_RH_C0006G9_7_1_HTTLKBBXX.IND2.fastq
touch subdir2/G630_RH_C0006G9_7_2_HTTLKBBXX.IND2.fastq
touch subdir2/G630_RH_C0006GL_8_1_HTTLKBBXX.IND2.fastq
touch subdir2/G630_RH_C0006GL_8_2_HTTLKBBXX.IND2.fastq


python --version
echo "python star-all.py --in-dir=. --out-dir=align # Cmd-line"
echo "python star-all.py                            # GUI"
