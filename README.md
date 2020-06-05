# METATRYP v. 2.0 Documentation
These tools are intended to help researchers answer questions about proteomics data. Written in Python 3.6 on a Unix operating system (linux/macOS).

Some specific questions that these tools help provide answers for include:

* Is tryptic peptide XXX found in Proteome A and Proteome B?
* How many tryptic peptides do Proteome A and Proteome B have in common?
* Are tryptic peptides XXX and YYY both found in Proteome A?
* How many tryptic peptides are there in Proteome A that differ from peptide XXX by n positions?
* What is the Lowest Common Ancestor of taxa containing peptide XXX?

METATRYP v. 2.0 has faster performance and higher throughput than [METATRYP v. 1.0](https://github.com/saitomics/metatryp) which it is based upon. A web server instance of METATRYP v. 2.0 using a marine microbe database can be found at [https://metatryp.whoi.edu/](https://metatryp.whoi.edu/), and a web instance using a coronavirus focused database can be found at [https://metatryp-coronavirus.whoi.edu/](https://metatryp-coronavirus.whoi.edu/).

METARYP 2.0 now supports the ingestion and search of Genomes, Metagenomes and Metagenome Assembled Genomes (MAGs) and calculation of the Least Common Ancestor (LCA) for each peptide. This repository is designed both to be used as a stand-alone application and in support of the Ocean Protein Portal (https://proteinportal.whoi.edu).

The METATRYP software consists of:
1. a database, stored in a standalone postgreSQL database.
2. python libraries, for processing, ingesting, and analyzing data. 
3. command-line scripts, to act as interfaces to the python libraries and the database.
***
* [Installation Instructions](https://github.com/WHOIGit/metatryp-2.0/wiki/Installation-Instructions)
    * [Download](https://github.com/WHOIGit/metatryp-2.0/wiki/Installation-Instructions#1-download-this-repository)
    * [Create virtual environment](https://github.com/WHOIGit/metatryp-2.0/wiki/Installation-Instructions#2-create-a-virtual-environment)
     [Setup PostgreSQL database](https://github.com/WHOIGit/metatryp-2.0/wiki/Installation-Instructions#3-setup-postgresql-to-build-database-from-scratch)
    * [Add METATRYP Schema to empty database](https://github.com/WHOIGit/metatryp-2.0/wiki/Installation-Instructions#4-add-metatryp--schema-to-the-empty-database)
    * [Run tests](https://github.com/WHOIGit/metatryp-2.0/wiki/Installation-Instructions#5-run-tests)
* [Quick Usage Guide](https://github.com/WHOIGit/metatryp-2.0/wiki/Quick-Usage-Guide)
    * [Digest and ingest data](https://github.com/WHOIGit/metatryp-2.0/wiki/Quick-Usage-Guide#1-digest-and-ingest-data)
    * [Query a data category for a peptide](https://github.com/WHOIGit/metatryp-2.0/wiki/Quick-Usage-Guide#2-query-a-data-category-for-a-peptide)
    * [List taxa ids](https://github.com/WHOIGit/metatryp-2.0/wiki/Quick-Usage-Guide#3-list-the-taxon-ids)
    * [Generate redundancy tables](https://github.com/WHOIGit/metatryp-2.0/wiki/Quick-Usage-Guide#4-generate-redundancy-tables)
    * [Remove taxa from database](https://github.com/WHOIGit/metatryp-2.0/wiki/Quick-Usage-Guide#5-remove-taxa-from-the-database)
* [Tutorial](https://github.com/WHOIGit/metatryp-2.0/wiki/Tutorial)
    * [Taxa in Tutorial Database](https://github.com/WHOIGit/metatryp-2.0/wiki/Tutorial/_edit#taxa-in-the-pre-built-tutorial-database)
***
When using METATRYP v. 2.0 please cite the following:

[Saunders, J. K., Gaylord, D., Held, N., Symmonds, N., et al., METATRYP v 2.0: Metaproteomic Least Common Ancestor Analysis for Taxonomic Inference Using Specialized Sequence Assemblies - Standalone Software and Web Servers for Marine Microorganisms and Coronaviruses. bioRxiv 2020, 2020.2005.2020.107490.](https://www.biorxiv.org/content/10.1101/2020.05.20.107490v1)
