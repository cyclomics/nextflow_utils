# nextflow_utils

This repo contains resuable processes and workflows for other nextflow pipelines within Cyclomics.

## Folders

### Consensus
This folder contains the consensus generators.

#### Input patern
For each fastq, containing at least 1 sequence to extract reads from we expect the following data in  this exact order:
1. Sample ID
1. Fastq ID
1. Fastq file object
1. Synthtic DNA subunit file object e.g. backbone
1. Reference genome file object

#### Output patern
1. Sample ID
1. Fastq ID
1. Fastq file object

### Primary
Here we analyze the sequence itself, eg. we look for the presence of barcode sequences, adapters or specific read lengths.
The input and output of these analyses are symetric, allowing easy chaining of different steps.

#### Input patern
For each fastq we expect the following data in  this exact order:
1. Sample ID
1. Fastq ID
1. Fastq file object
Another tuple or value, specific to the analysis, eg a value or a file

#### Output patern
1. Sample ID
1. Fastq ID
1. Fastq file object

### Secondary
Abc

#### Input patern
For each fastq we expect the following data in  this exact order:
1. Sample ID
1. Fastq ID
1. Fastq file object
Another tuple or value, specific to the analysis, eg a value or a file

#### Output patern
1. Sample ID
1. Fastq ID
1. New format file, e.g. bam file.

### Tertiary
Abc

### Quaternary
Abc
