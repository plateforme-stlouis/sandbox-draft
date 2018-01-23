
Attempt to glue a RNA-seq workflow.


1. Fetch the data samples,
```bash
    bash fetch-data.sh
```
It will create the directory `data-ftp.ebi.ac.uk` with 4 FASTQ files.

2. Install Miniconda and add Bioconda
```bash
    curl -C - -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    conda config --add channels r
    conda config --add channels conda-forge
    conda config --add channels bioconda
```
Edit your `~/.condarc` to obtain (because the order matters),
```
channels:
  - bioconda
  - conda-forge
  - r
  - defaults
```

3. Create one R environment and go inside
```bash
    conda create -n r-env r-base r-xml -y
    # should take some time, be patient...
    source activate r-env
```
Then, you can return whenever to this environment with the command
`source activate r-env`.

4. Install locally the R dependencies
```bash
    Rscript install-deps.R
    Rscript check-deps.R
```


Requirements
------------

STAR

R:
```bash
    grep library *.R | cut -d'(' -f2 | cut -d')' -f1 | sort | uniq
```
 - Rsamtools
 - TxDb.Hsapiens.UCSC.hg19.knownGene
 - GenomicFeatures
 - GenomicAlignments
 - BiocParallel
 - edgeR
 - limma
 - Glimma
 - gplots
 - RColorBrewer
 - org.Hs.eg.db
 - xlsx
