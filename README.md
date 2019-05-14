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

These instructions assume that you are running within a unix environment that has python3.X, and that Anaconda or Miniconda is installed on the system.

1. Download this Repository: https://github.com/WHOIGit/metatryp-2.0/archive/master.zip
2. Unpack the Repository into a folder: Unzip it.
3. Create a virtual environment: Navigate into the unzipped directory you created above. It should contain this README.md file. Then run this command: conda env create -f environment.yml . This will create a new conda virtual environment directory named 'metatryp2_env'.
4. Activate the environment: conda activate metatryp2_env

Further instructions to follow.

## Quick Usage Guide

### 1: Digest and ingest genomic data
1. Run the script bin/digest_and_ingest.sh with FASTA proteome files you wish to digest and ingest. e.g.:
````
bin/digest_and_ingest.sh file1.fasta file2.fasta ...
````

This script reads the FASTA files, and runs digestions on their sequences. You should see a fair amount of output as these files are processed.


### 2.  Query sequences for taxon matches:
You can query the database for exact and fuzzy matches*. e.g.:
````
bin/query_by_sequence.sh --sequence MGFPCNR --max-distance 1
````
*fuzzy matches is currently not working

You can query by other types such as metagenomes or specialized assemblies by using the --type parameter

For example:
````
bin/query_by_sequence.sh --sequence MGFPCNR --type sa
````
will query for Specialized Assemblies such as MAGs.

````
bin/query_by_sequence.sh --sequence MGFPCNR --type all
````
will query across genomes, metagenomes and specialized assemblies

--type parameters (default is genomes)

g - genomes

m - meta

sa - specialized assemblies

all - all types

### 3: Generate redundancy tables
1.: See available taxon ids by querying DB: e.g. 
````
bin/list_taxon_ids.sh
````
2.: Generate redundancy tables for groups of taxons e.g.
````
bin/generate_redundancy_tables.sh --taxon-ids syn8102 syn7502 syn7503 --output-dir exampleRedundancyTables
````

Note that you can also specify a file that contains a list of taxon IDs, e.g

````
bin/generate_redundancy_tables.sh --taxon-id-file taxon_id_list.txt --output-dir exampleRedundancyTables
````

3.: View resulting files in exampleRedundancyTables
    - counts.csv contains counts of redundant peptides
    - percents.csv contains the values in counts.csv, divided by the number of unique peptides in the *union* of digestions of a taxa pair.

### 3.(optional, expected to occur rarely): Clear data for a given set of taxa.
If you wish to **delete** data for a given set of taxa in the db, run a command like this:
````
bin/clear_taxon_data.sh --taxon-ids croc5801
````


