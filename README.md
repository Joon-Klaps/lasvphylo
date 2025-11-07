# LASV-phylo V2.0

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A525.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![test status](https://github.com/Joon-Klaps/lasvphylo/actions/workflows/test.yml/badge.svg)](https://github.com/Joon-Klaps/lasvphylo/actions/workflows/test.yml)

## Introduction

**nf-core/lasvphylo** is a bioinformatics "best-practice" analysis pipeline for creating a high quality alignment of LASV segements and a maximum likelihood tree. I use this pipeline for a fast phylogentic analysis of newly identified LASV cases/outbreaks to determine the rise of new clades or interesting mutations.

## Pipeline summary

1. Orient and isolate the genes of LASV ([`MAFFT`](https://mafft.cbrc.jp/alignment/software/))
2. Remove the reference sequences from the alignment ([`SeqTk`](https://github.com/lh3/seqtk))
3. Concatenate the genes to each other ([`SeqKit`](https://bioinf.shenwei.me/seqkit/))
4. Perform a **constrained** tree search using a previous tree and new sequences ([`IQ-TREE`](http://www.iqtree.org/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=25.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)) or [`Conda`](https://conda.io/miniconda.html).

3. Download the pipeline (git clone) and give in the correct variables :

   ```bash
   nextflow run Joon-Klaps/lasvphylo -profile <docker|conda|singularity> -c <CONFIG> --outdir <OUTDIR>
   ```

Where `<CONFIG>` is a configuration file (see `conf/` for examples) and `<OUTDIR>` is the desired output directory.

Example configuration file:

```groovy
params {
    // Maximum resource limits
    max_cpus   = 6
    max_memory = '8.GB'
    max_time   = '1.h'

    // Necessary inputs
    input_L      = '/path/to/input/sequences.L.fasta'
    input_S      = '/path/to/input/sequences.S.fasta'

    // Optional Inputs - if none provided, everything will be run from scratch
    alignment_L = '/path/to/existing/alignment.L.fasta'
    alignment_S = '/path/to/existing/alignment.S.fasta'
    tree_L      = '/path/to/existing/tree.L.treefile'
    tree_S      = '/path/to/existing/tree.S.treefile'

    // Output directory
    outdir      = 'testing'
}
```
