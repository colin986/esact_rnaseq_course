[![DOI](https://zenodo.org/badge/225317346.svg)](https://zenodo.org/badge/latestdoi/225317346)

# ESACT Metabolic and Bioprocess Modelling for Animal Cells Course 2023

## Introduction
---
This respository enables the reproduction of the differential expression analysis of the RNA-seq data described in: 

Tzani *et al.* **Sub physiological temperature induces pervasive alternative splicing in Chinese hamster ovary cells**
*Biotechnology and Bioengineering 2020* [https://doi.org/10.1002/bit.27365]( https://doi.org/10.1002/bit.27365)

### Data availability:
[NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA593052/)  
[European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/view/PRJNA593052)

### Dependancies
All the programmes must be added to the PATH to run the workflow
- [Python 2.7.12](https://www.python.org/download/releases/2.7/)
- [trimmomatic 0.36](http://www.usadellab.org/cms/?page=trimmomatic) 
- [cutadpat 1.18](https://cutadapt.readthedocs.io/en/stable/)
- [STAR-2.7.11a](https://github.com/alexdobin/STAR)

## RNASeq data preprocesssing
---
### Download the data from ENA
This is a simple way to dowload from ENA, for higher speed download use the Aspera client
Total data download size: **~95G** 
```bash
mkdir -p data/ena
wget -q "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/057/SRR10572657/*" -P data/ena 
wget -q "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/058/SRR10572658/*" -P data/ena
wget -q "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/059/SRR10572659/*" -P data/ena 
wget -q "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/060/SRR10572660/*" -P data/ena 
wget -q "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/061/SRR10572661/*" -P data/ena 
wget -q "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/062/SRR10572662/*" -P data/ena 
wget -q "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/063/SRR10572663/*" -P data/ena 
wget -q "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/064/SRR10572664/*" -P data/ena
```
```bash
conda activate esact_rnaseq
```

```bash
fastqc -t 70 data/ena/* -o data/quality_test/before/
multiqc data/quality_test/ --filename raw_data_qc
```

```bash
mkdir data/trimmed
IN_DIR=data/ena
OUT_DIR=data/trimmed

cat sample_info.txt | cut -f 2 | tail -n 8 | while read SAMPLE_ID; do
    cutadapt  \
    -A AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
    -o $OUT_DIR/"$SAMPLE_ID"_1.fastq.gz \
    -p $OUT_DIR/"$SAMPLE_ID"_2.fastq.gz \
    $IN_DIR/"$SAMPLE_ID"_1.fastq.gz $IN_DIR/"$SAMPLE_ID"_2.fastq.gz
done
```

```bash
mkdir -p data/preprocessed/paired data/preprocessed/unpaired
IN_DIR=data/trimmed
OUT_DIR=data/preprocessed

cat sample_info.txt | cut -f 2 | tail -n 8 | while read SAMPLE_ID; do
    trimmomatic PE \
    -threads 70 \
    $IN_DIR/"$SAMPLE_ID"_1.fastq.gz $IN_DIR/"$SAMPLE_ID"_2.fastq.gz \
    $OUT_DIR/paired/"$SAMPLE_ID"_1.fastq.gz $OUT_DIR/unpaired/"$SAMPLE_ID"_1.fastq.gz \
    $OUT_DIR/paired/"$SAMPLE_ID"_2.fastq.gz $OUT_DIR/unpaired/"$SAMPLE_ID"_2.fastq.gz \
    SLIDINGWINDOW:4:20 \
    MINLEN:36 \
    -trimlog $OUT_DIR/"$SAMPLE_ID".trimmomatic.log
done
```

```bash
fastqc -t 70 data/preprocessed/* -o data/quality_test/after/
multiqc data/quality_test/after --filename preprocessed_qc
```


## Download CHO cell PICR reference genome

```bash
mkdir reference_genome

genome_ftp=https://ftp.ensembl.org/pub/release-110/fasta/cricetulus_griseus_picr/
wget $genome_ftp/dna/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel.fa.gz -P reference_genome

gtf_ftp=https://ftp.ensembl.org/pub/release-110/gtf/cricetulus_griseus_picr
wget $gtf_ftp/Cricetulus_griseus_picr.CriGri-PICRH-1.0.110.gtf.gz -P reference_genome

gunzip reference_genome/*
```


## Create STAR genome index 

```bash
mkdir reference_genome/star_index

FASTA=reference_genome/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel.fa
GTF=reference_genome/Cricetulus_griseus_picr.CriGri-PICRH-1.0.110.gtf

STAR --runThreadN 70 \
     --runMode genomeGenerate \
     --sjdbOverhang 124\
     --genomeChrBinNbits 16 \
     --genomeDir reference_genome/star_index \
     --genomeFastaFiles $FASTA \
     --sjdbGTFfile $GTF
```

## Align reads to index
```bash
mkdir -p data/mapped
OUT_DIR=data/mapped
IN_DIR=data/preprocessed/paired

cat sample_info.txt | cut -f 2 | tail -n 8 | while read SAMPLE_ID; do

    STAR \
    --runThreadN 16 \
    --readFilesIn $IN_DIR/"$SAMPLE_ID"_1.fastq.gz $IN_DIR/"$SAMPLE_ID"_2.fastq.gz \
    --genomeDir reference_genome/star_index \
    --readFilesCommand gunzip -c \
    --outFileNamePrefix $OUT_DIR/"$SAMPLE_ID" \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic
done
```
```bash
mkdir -p data/counts
OUT_DIR=data/counts
IN_DIR=data/mapped

GTF=reference_genome/Cricetulus_griseus_picr.CriGri-PICRH-1.0.110.gtf

cat sample_info.txt | cut -f 2 | tail -n 8 | while read SAMPLE_ID; do
htseq-count -r pos -f bam -i gene_id -s reverse $IN_DIR/"$SAMPLE_ID"Aligned.sortedByCoord.out.bam $REF_GTF > "$OUT_DIR"/"$SAMPLE_ID".counts
done
```