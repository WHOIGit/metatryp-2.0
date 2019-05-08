Saito Lab: METATRYP 2.0 Proteomics Toolkit
==============
NOTE: This documentation is currently in progress and is incomplete.

This repository contains a set of computational tools for analyzing proteomics data.

It is revised and expanded from METATRYP 1.0 found here: https://github.com/saitomics/metatryp

METARYP 2.0 now supports the ingestion and search of Genomes, Metagenomes and Metagenome Assembled Genomes (MAGs).  This reposititory is designed both to be used as a stand alone application and in support of the Ocean Protein Portal (https://proteinportal.whoi.edu).

## General System Description

These tools are intended to help researchers answer questions about proteomics data. 

Some specific questions that these tools help provide answers for include:

- Is tryptic peptide XXX found in Proteome A and Proteome B?
- How many tryptic peptides do Proteome A and Proteome B have in common?
- Are tryptic peptides XXX and YYY both found in Proteome A?
- How many tryptic peptides are there in Proteome A that differ from peptide XXX by n positions?

The analysis tools consist of:

1): a database, stored in a PostgreSQL database.  
2): python libraries, for processing, ingesting, and analyzing data 
3): command-line scripts, to act as interfaces to the python libraries and the database.

## Installation

These instructions assume that you are running within a unix environment that has python3.X .

Further instructions to follow.

## Quick Usage Guide

### 1: Digest and ingest genomic data
1. Run the script bin/digest_and_ingest.sh with FASTA proteome files you wish to digest and ingest. e.g.:
````
bin/digest_and_ingest.sh file1.fasta file2.fasta ...
````

This script reads the FASTA files, and runs digestions on their sequences. You should see a fair amount of output as these files are processed.


### 2.  Query sequences for taxon matches:
You can query the database for exact and fuzzy matches. e.g.:
````
bin/query_by_sequence.sh --sequence MGFPCNR --max-distance 1
````

### 3.(optional, expected to occur rarely): Clear data for a given set of taxa.
If you wish to **delete** data for a given set of taxa in the db, run a command like this:
````
bin/clear_taxon_data.sh --taxon-ids croc5801
````


