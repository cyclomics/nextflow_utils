# nextflow_utils

This repo contains resuable processes and workflows for other nextflow pipelines within Cyclomics.

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://code.askimed.com/nf-test)
![Main CI status](https://github.com/cyclomics/nextflow_utils/actions/workflows/check_input_data.yml/badge.svg)
![development CI status](https://github.com/cyclomics/nextflow_utils/actions/workflows/check_input_data.yml/badge.svg)

## Folders

### Consensus
This folder contains the consensus generators.

#### Input patern
For each fastq, containing at least 1 sequence to extract reads from we expect the following data in  this exact order:
1. Sample ID
1. Fastq ID
1. Fastq file object
And optionally the following as seperate inputs:
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


# Testing lessons learned with nf-test:

1. Dont compare to large objects, as this makes the stdout hard to read in case of a test failure:
    ```json
    def first_line = lines.get(0)
    assert first_line == "@64fee029-4..."
    ```
    is much better than:
    ```json
    assert lines.get(0) == "@64fee029-4..."
    ```
    as the later will cause the entire lines object to be dumped into stdout, while the prior will only put the first line into your stdout.

2. Regex is kinda weird and undocumented: this fails:
    ```json
    assert output_file == ".*/empty_filtered.fastq"
    ```
    Buth this passes:
    ```json
    assert output_file ==~ ".*/empty_filtered.fastq"
    ```

