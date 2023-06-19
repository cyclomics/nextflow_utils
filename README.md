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

### Sequence analysis
Here we analyze the sequence itself, eg. we look for the presence of barcode sequences, adapters or specific read lengths.
The input and output of these analyses are symetric, allowing easy chaining of different steps.

#### Naming convention
tbd

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

### parse_covert
Workflows and modules to convert data structures from one to another. A good example of this is mapping against a reference genome. 

#### Naming convention
tbd

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

### Post Hoc
Code to perform additional actions, wheter its enriching or filtering, on the input data. eg Checking cosmic id's in a vcf using an API. 

### Reporting
Code for visualization and reporting of the results.

