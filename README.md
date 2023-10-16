# ESACT Metabolic and Bioprocess Modelling for Animal Cells Course 2023

## Introduction

---
This respository enables the reproduction of the differential expression analysis of the RNA-seq data described in: 

Tzani *et al.* **Sub physiological temperature induces pervasive alternative splicing in Chinese hamster ovary cells** *Biotechnology and Bioengineering 2020* [https://doi.org/10.1002/bit.27365]( https://doi.org/10.1002/bit.27365)

## Preparation

---

### Dependancies

All the software must be in the PATH to run this workflow

- [FastQC v0.12.0](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [multiqc v1.16](https://multiqc.info)
- [TrimGalore](https://github.com/FelixKrueger/TrimGalore)
- [STAR-2.7.11a](https://github.com/alexdobin/STAR)
- [HTSeq v2.0.4](https://htseq.readthedocs.io/en/latest/)

### Data availability

Information on the deposited data can be found here:

[NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA593052/)  
[European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/view/PRJNA593052)

### Download the data from ENA

This is a straightforward way to dowload from ENA, for higher speed download use the Aspera client.

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

## RNASeq data preprocesssing

---

### Initial QC of reads

FASTQC serves as an initial quality assessment tool for individual RNA-seq samples, while MULTIQC enhances the efficiency of data quality assessment by aggregating and visualizing results from multiple samples or datasets. Together, these tools play a crucial role in ensuring the reliability and accuracy of RNA-seq data analysis.

```bash
# sample by sample evaluation 
fastqc data/ena/* -o data/quality_test/before/

# aggregate the results
multiqc data/quality_test/ --filename raw_data_qc
```

### Filter and trim reads

The TrimGalore tool is used to remove any adapter sequences and trim low quality bases

```bash

mkdir data/trim_galore

IN_DIR=data/ena
OUT_DIR=data/trim_galore

# loop over samples
cat sample_info.txt | cut -f 2 | tail -n 8 | while read SAMPLE_ID; do
    trim_galore --paired $IN_DIR/"$SAMPLE_ID"_1.fastq.gz $IN_DIR/"$SAMPLE_ID"_2.fastq.gz \
    --illumina -o $OUT_DIR -q 30 \
    -j 8 --fastqc --fastqc_args "-o data/quality_test/after/"
done 
```

### Final quality assesement

Compare the post-trimming quality of the data

```bash
# sample by sample evaluation 
fastqc data/preprocessed/paired/* -o data/quality_test/after/

# aggregate the results 
multiqc data/quality_test/after --filename preprocessed_qc
```

## Mapping reads to reference genome

In this tutorial we use the latest assembly of the Chinese hamster PICRH genome availiable through ENSEMBL. The corresponding GTF annotation file is downloaded. 

---

### Download CHO cell PICR reference genome

```bash
mkdir reference_genome

genome_ftp=https://ftp.ensembl.org/pub/release-110/fasta/cricetulus_griseus_picr/
wget $genome_ftp/dna/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel.fa.gz -P reference_genome

gtf_ftp=https://ftp.ensembl.org/pub/release-110/gtf/cricetulus_griseus_picr
wget $gtf_ftp/Cricetulus_griseus_picr.CriGri-PICRH-1.0.110.gtf.gz -P reference_genome

gunzip reference_genome/*
```

### Create STAR genome index

In this tutorial, we will use the STAR to align our reads for each sample to the genome. First an index must be constructed using the FASTA and GTF files downloaded from ENSEMBL

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

### Alignment 

Map the reads to the genome index

```bash

mkdir -p data/mapped

OUT_DIR=data/mapped
IN_DIR=data/trim_galore

# loop over samples
cat sample_info.txt | cut -f 2 | tail -n 8 | while read SAMPLE_ID; do
    STAR \
    --runThreadN 16 \
    --readFilesIn $IN_DIR/"$SAMPLE_ID"_1_val_1.fq.gz $IN_DIR/"$SAMPLE_ID"_2_val_2.fq.gz \
    --genomeDir reference_genome/star_index \
    --readFilesCommand gunzip -c \
    --outFileNamePrefix $OUT_DIR/"$SAMPLE_ID" \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic
done
```

### Mapping statistics

MultiQC can also aggregate the mapping results

```bash
multiqc data/mapping/*Log.final.out --filename mapping_stats
```

## Counting reads aligned to features

---

### HTSeq count

```bash

mkdir -p data/counts

OUT_DIR=data/counts
IN_DIR=data/mapped

GTF=reference_genome/Cricetulus_griseus_picr.CriGri-PICRH-1.0.110.gtf

# loop over samples
cat sample_info.txt | cut -f 2 | tail -n 8 | while read SAMPLE_ID; do
    # Index the bam file
    samtools index $IN_DIR/"$SAMPLE_ID"Aligned.sortedByCoord.out.bam
    # count 
    htseq-count \
    -r pos -f bam -i gene_id -s reverse \
    $IN_DIR/"$SAMPLE_ID"Aligned.sortedByCoord.out.bam \
    $GTF > \
    "$OUT_DIR"/"$SAMPLE_ID".counts

done
```

### Count statistics

Finally multiqc is used to evaluate the counting process

```bash
multiqc data/counts/* --filename count_stats
```
