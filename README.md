# LASV-phylo V1.0

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**nf-core/lasvphylo** is a bioinformatics "best-practice" analysis pipeline for creating a high quality alignment of LASV segements and a maximum likelihood tree. I use this pipeline for a fast phylogentic analysis of newly identified LASV cases/outbreaks to determine the rise of new clades or interesting mutations.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

<!-- TODO create a graphical overview -->

1. Orient and isolate the genes of LASV ([`MAFFT`](https://mafft.cbrc.jp/alignment/software/))
2. Remove the reference sequences from the alignment ([`SeqTk`](https://github.com/lh3/seqtk))
3. Concatenate the genes to each other ([`SeqKit`](https://bioinf.shenwei.me/seqkit/))
4. Align the concatenated genes to an existing alignment ([`MUSCLE`](https://www.drive5.com/muscle))
5. Perform a __constrained__ tree search using a previous tree and new sequences ([`IQ-TREE`](http://www.iqtree.org/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.
    > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
    > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

3. Download the pipeline (got clone) and give in the correct variables :

   ```bash
   nextflow run lasvphylo/main.nf -profile YOURPROFILE -c <CONFIG> --outdir <OUTDIR>
   ```
   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. You can chain multiple config profiles in a comma-separated string. In this config file you'll have to specify the following variables:

   - input_S : S segment of the new sequence(s)
   - input_L : L segment of the new sequence(s)
   - input_id_S : id for S segment files
   - input_id_L : id for L segment files
   - alignment_S : Previous alignment of S segment
   - alignment_L : Previous alignment of L segment
   - tree_S : Previous tree of S segment
   - tree_L : Previous tree of L segment





