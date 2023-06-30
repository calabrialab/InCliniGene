# InCliniGene
Graph Database solution for Clonal tracking with vector integration sites (Integrated Clinical Genomics DB)

## Introduction

InCliniGene is the first database of comprehensive archive of clones from gene therapy datasets using integrating viral vectors.
This repository contains data and database structure for clonal tracking in gene therapy applications. In particular, the project InCliniGene aims at collecting the output of **gamma-TRIS** (Calabria et al, Bioinformatics 2020, [DOI 10.1093/bioinformatics/btz747](https://doi.org/10.1093/bioinformatics/btz747)) which is able to identify integation sites also in repetitive or low complexity genomic regions, and building a graph database of all the identified clones for tracking and data mining purposes by only data querying.

## Exploring the repository

The folder *"gtris_output/clonal_expansion"* contains the output files of gamma-TRIS for our validation assay, that is an experimental in vitro dataset of vector integration sites (IS), retrieved by the protocol Sonication Linker Merdiate PCR (described in Cesana et al, Nature Medicine 2021) and sequenced by Illumina paired-ends. We then run gamma-TRIS to identify the IS after trimming the vector sequences with VISPA2 (Spinozzi et al, BMC Bioinformatics 2017 [DOI 10.1186/s12859-017-1937-9](https://doi.org/10.1186/s12859-017-1937-9)).

The folder *"import_data_graphdb"* contains the scripts for data import into the graph database. The database export is here stored in the folder *"graph_db_export"*.

The query output to build the tracking matrix is reported in the folder *"query_graphdb"*.

We compared the results of gamma-TRIS in terms of output tracking matrix with the output matrix exported from **InCliniGene** and we reported the results in the filder *"compare_dataset"*. Moreover, we compared the IS landing in a repetitive region between the quantification of gamma-TRIS and VISPA2 to demonstrate the reliability of gamma-TRIS and the importance of using its results in all clonal tracking analyses.

Further details will be available in the publication, once released.

## Tutorials and Demo

We released a [tutorial page](TUTORIAL_explore.md) to guide all users through the database configuration, import and query.

We also added a demo video ([external link via Zenodo](https://doi.org/10.5281/zenodo.8102448)) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8102448.svg)](https://doi.org/10.5281/zenodo.8102448)
